SUBROUTINE P_Change_Decomposition(Decmp_name,nato,vpx,vpy,vpz,nstart&
     &,nend,nstart1,nend1,node,nprocs)

!!$***********************************************************************
!!$   Time-stamp: <2009-03-12 15:53:12 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 19 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
  REAL(8) :: vpx(*),vpy(*),vpz(*)
  INTEGER :: nstart,nend,nstart1,nend1
  INTEGER :: node,nprocs,nato
  CHARACTER(*) :: Decmp_name

!!$----------------------- VARIABLES IN COMMON --------------------------*

#if defined PARALLEL   
  INCLUDE 'mpif.h'
  include 'mpi_size.h'

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: exchg,i
  INTEGER, DIMENSION (:), ALLOCATABLE :: isize
  INTEGER, DIMENSION (:,:), ALLOCATABLE :: partners
  INTEGER, DIMENSION (:), ALLOCATABLE :: starts,ends,starts1,ends1&
       &,sfirst
  INTEGER :: istatus(MPI_STATUS_SIZE),source,dest,k&
       &,first,last,Buffer_Dim

  REAL(8), DIMENSION (:), ALLOCATABLE :: buf
  CHARACTER, DIMENSION(:), ALLOCATABLE :: Buffer_MPI
  INTEGER :: Buffer_p

  INTEGER :: ierr,Buffer_pp,Buffer_pa
  INTEGER :: count,j,nnst,nned
  INTEGER, DIMENSION (:), ALLOCATABLE :: nstart_o,nend_o

  TYPE NODES
     INTEGER :: exchg_n,nstart_n,nend_n,nstart1_n,nend1_n
     INTEGER, DIMENSION (:,:), Allocatable :: partner_n
     INTEGER, DIMENSION (:), Allocatable :: sfirst_n,isize_n
     CHARACTER(80) :: label
     TYPE (NODES), POINTER :: next
  END TYPE NODES
  TYPE(NODES), POINTER, SAVE :: list
  TYPE(NODES), POINTER, SAVE :: nfirst

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  CALL Inquire_Decmp(list,nfirst)
  ALLOCATE(nstart_o(exchg),nend_o(exchg))

  ALLOCATE(buf(3*nato))
  Buffer_Dim=SUM(isize(1:exchg))  
  CALL MPI_PACK_SIZE(Buffer_Dim,MPI_REAL8,MPI_COMM_WORLD,Buffer_p&
       &,ierr)
  Buffer_p=Buffer_p+MPI_BSEND_OVERHEAD
  ALLOCATE(Buffer_MPI(Buffer_p))

  Buffer_pa=Buffer_p

  count=0
  DO k=1,exchg
     dest=partners(1,k)
     source=partners(2,k)
     first=sfirst(k)
     last=first+isize(k)-1
     IF(node == dest) THEN
        count=count+1
        CALL MPI_RECV(buf(first),isize(k),MPI_REAL8,source,0&
             &,MPI_COMM_WORLD, istatus ,ierr)
        nstart_o(count)=(first-1)/3+1
        nend_o(count)=isize(k)/3-1+nstart_o(count)
     ELSE IF(node == source) THEN
        DO i=nstart,nend
           buf((i-1)*3+1)=vpx(i)
           buf((i-1)*3+2)=vpy(i)
           buf((i-1)*3+3)=vpz(i)
        END DO
        CALL MPI_BUFFER_ATTACH(Buffer_MPI,Buffer_pa,ierr)
        CALL MPI_BSEND(buf(first),isize(k),MPI_REAL8,dest,0,&
             & MPI_COMM_WORLD,ierr)
        CALL MPI_BUFFER_DETACH(Buffer_MPI,Buffer_pa,ierr)
     END IF
  END DO
  DO j=1,count
     DO i=nstart_o(j),nend_o(j)
        vpx(i)=buf((i-1)*3+1)
        vpy(i)=buf((i-1)*3+2)
        vpz(i)=buf((i-1)*3+3)
     END DO
  END DO
  DEALLOCATE(isize,partners,sfirst,nstart_o,nend_o)

  DEALLOCATE(Buffer_MPI,buf)
        
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
CONTAINS
  SUBROUTINE Get_Partners(starts,ends,starts1,ends1,sfirst,partners,isize,ex)
    IMPLICIT none
    INTEGER :: ex
    INTEGER :: partners(:,:),isize(:),starts(:),starts1(:),ends(:)&
         &,ends1(:),sfirst(:)

    INTEGER :: nsti,nstj,neni,nenj,i,j
    INTEGER :: nstart,nend

    ex=0
    DO i=1,nprocs
       nsti=starts(i)
       neni=ends(i)
       DO j=1,nprocs
          IF(j .NE. i) THEN
             nstj=starts1(j)
             nenj=ends1(j)
             IF(nstj .LE. neni .AND. nenj .GE. nsti) THEN
                nend=nenj
                IF(nend .GT. neni) nend=neni
                nstart=nstj
                IF(nstart .LT. nsti) nstart=nsti
                ex=ex+1
                partners(2,ex)=i-1
                partners(1,ex)=j-1
                isize(ex)=3*(nend-nstart+1)
                sfirst(ex)=3*(nstart-1)+1
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE Get_Partners
  SUBROUTINE Inquire_Decmp(list,nfirst)
    IMPLICIT none
    TYPE(NODES), POINTER:: list
    TYPE(NODES), POINTER:: nfirst
    LOGICAL :: there

    
    there=.FALSE.
    IF(.NOT. ASSOCIATED(list)) THEN
       ALLOCATE(list)
       NULLIFY(list%next)
       nfirst=>list
    ELSE
       count=0
       list=>nfirst
       DO WHILE (ASSOCIATED(list%next))
          IF(Decmp_Name == list%label) THEN
             nstart=list%nstart_n
             nend=list%nend_n
             nstart1=list%nstart1_n
             nend1=list%nend1_n
             exchg=list%exchg_n
             ALLOCATE(isize(exchg))
             ALLOCATE(sfirst(exchg))
             ALLOCATE(partners(2,exchg))
             partners(1:2,1:exchg)=list%partner_n(1:2,1:exchg)
             sfirst(1:exchg)=list%sfirst_n(1:exchg)
             isize(1:exchg)=list%isize_n(1:exchg)
             there=.TRUE.
          END IF
          list => list%next
          count=count+1
       END DO
    END IF
    IF(.NOT. there) THEN
       ALLOCATE (starts(nprocs),ends(nprocs))
       ALLOCATE (starts1(nprocs),ends1(nprocs))
       ALLOCATE(isize(nprocs*nprocs))
       ALLOCATE(sfirst(nprocs*nprocs))
       ALLOCATE(partners(2,nprocs*nprocs))
       CALL MPI_ALLGATHER(nstart,1,MPI_INTEGER4,starts,1,MPI_INTEGER4&
            &,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHER(nend,1,MPI_INTEGER4,ends,1,MPI_INTEGER4&
            &,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHER(nstart1,1,MPI_INTEGER4,starts1,1,MPI_INTEGER4&
            &,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHER(nend1,1,MPI_INTEGER4,ends1,1,MPI_INTEGER4&
            &,MPI_COMM_WORLD,ierr)
       
       
       CALL Get_Partners(starts,ends,starts1,ends1,sfirst,partners,isize,exchg)
       
       list%nstart_n=nstart
       list%nend_n=nend
       list%nstart1_n=nstart1
       list%nend1_n=nend1
       list%exchg_n=exchg
       list%label=Decmp_Name
       
       ALLOCATE(list%partner_n(2,exchg),list%sfirst_n(exchg),list&
            &%isize_n(exchg))
       list%partner_n(1:2,1:exchg)=partners(1:2,1:exchg)
       list%sfirst_n(1:exchg)=sfirst(1:exchg)
       list%isize_n(1:exchg)=isize(1:exchg)
       DEALLOCATE(starts,ends,starts1,ends1)
       
       ALLOCATE(list%next)
       list=>list%next
       NULLIFY(list%next)
    END IF
  END SUBROUTINE Inquire_Decmp

#endif

END SUBROUTINE P_Change_Decomposition
