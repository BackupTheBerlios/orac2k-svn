SUBROUTINE P_Comm_Intra(nstart,nend,node,nprocs,x,y,z,xc,yc,zc,nato&
     &,atomp,lstrtch,lstretch,lbnd_x,lbndg,lbend,lbndg_x,ltor,ltors&
     &,ltor_x,litr,litor,litr_x,int14,int14p,int14_x,int13,int13p&
     &,int13_x,ingrp,ingrpp,ingrp_x)

!!$***********************************************************************
!!$   Time-stamp: <2009-03-09 12:45:06 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Dec 11 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: nstart,nend,node,nprocs
  REAL(8) :: x(*),y(*),z(*),xc(*),yc(*),zc(*)
  INTEGER, OPTIONAL :: nato,atomp(*),lstrtch(2,*),lstretch,lbnd_x(*)&
       &,lbndg(3,*)&
       &,lbend,lbndg_x(*),ltor(4,*),ltors,ltor_x(*),litr(4,*),litor&
       &,litr_x(*),int14(2,*),int14p,int14_x(*),int13(2,*),int13p&
       &,int13_x(*),ingrp(2,*),ingrpp,ingrp_x(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

#if defined PARALLEL
  INCLUDE 'mpif.h'
  include 'mpi_size.h'

!!$------------------------- LOCAL VARIABLES ----------------------------*

  CHARACTER(1), DIMENSION(:), ALLOCATABLE :: Buffer_MPI
  INTEGER, SAVE :: Buffer_p,natoa
  INTEGER, DIMENSION (:), ALLOCATABLE :: starts,ends
  TYPE NODES
     INTEGER :: natoms
     INTEGER, DIMENSION (:), POINTER :: atoms
     TYPE (NODES), POINTER :: next
  END TYPE NODES
  TYPE LISTA
     INTEGER :: natoms
     INTEGER, DIMENSION (:), POINTER :: atoms
  END TYPE LISTA
  TYPE(LISTA), DIMENSION (:), ALLOCATABLE, SAVE :: Takes,Takes_m&
       &,Gives,Gives_m
  TYPE(NODES), POINTER, SAVE :: list,nfirst
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: ipartners
  INTEGER, SAVE :: exchg

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  IF(PRESENT(nato)) THEN
     natoa=nato
     CALL Setup_Comm
     RETURN
  END IF

  CALL Comm_Particles(x,y,z,Takes,Gives)
  CALL Comm_Particles(xc,yc,zc,Takes_m,Gives_m)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

CONTAINS
  SUBROUTINE Comm_Particles(x,y,z,Takes,Gives)
!!$***********************************************************************
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Dec 11 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8), DIMENSION (*) :: x,y,z
    TYPE(LISTA), DIMENSION (:) :: Takes,Gives

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8), DIMENSION(:), ALLOCATABLE :: Rbuff_r,Rbuff_s

    INTEGER :: Dims_Ibuff,count,i,j,jj,m,mp,mp3,Takes_Dim,Gives_Dim
    INTEGER :: ierr,k,dest,source,istatus(MPI_STATUS_SIZE),ii

!!$------------------------- LOCAL VARIABLES ----------------------------*
    
    Takes_Dim=SUM(Takes(1:nprocs)%natoms)*3
    Gives_Dim=SUM(Gives(1:nprocs)%natoms)*3

    ALLOCATE(Rbuff_r(Takes_Dim),Rbuff_s(Gives_Dim))
    mp=Takes_Dim+Gives_Dim

    CALL MPI_PACK_SIZE(mp,MPI_REAL8,MPI_COMM_WORLD,Buffer_p&
         &,ierr)
    Buffer_p=Buffer_p+MPI_BSEND_OVERHEAD
    ALLOCATE(Buffer_MPI(Buffer_p))

    DO k=1,exchg
       dest=ipartners(1,k)
       source=ipartners(2,k)
       IF(node == dest) THEN
          mp=Takes(source+1)%natoms
          mp3=mp*3
          CALL MPI_RECV(Rbuff_r,mp3,MPI_REAL8,source,0,MPI_COMM_WORLD&
               &, istatus ,ierr)
          DO ii=1,mp
             i=Takes(source+1)%atoms(ii)
             x(i)=Rbuff_r((ii-1)*3+1)
             y(i)=Rbuff_r((ii-1)*3+2)
             z(i)=Rbuff_r((ii-1)*3+3)
          END DO
       ELSE IF(node == source) THEN
          mp=Gives(dest+1)%natoms
          mp3=mp*3
          DO ii=1,mp
             i=Gives(dest+1)%atoms(ii)
             Rbuff_s((ii-1)*3+1)=x(i)
             Rbuff_s((ii-1)*3+2)=y(i)
             Rbuff_s((ii-1)*3+3)=z(i)
          END DO
          CALL Buffer_Attach
          CALL MPI_BSEND(Rbuff_s,mp3,MPI_REAL8,dest,0,MPI_COMM_WORLD,ierr)
          CALL Buffer_Detach
       END IF
    END DO
    
    DEALLOCATE(Buffer_MPI)

    DEALLOCATE(Rbuff_r,Rbuff_s)

  END SUBROUTINE Comm_Particles

!!$------------------------------------------------------------------------

  SUBROUTINE Setup_Comm
    IMPLICIT none 
    INTEGER :: i,j,count,nst,ned,nst_m,ned_m,inode,jj,m,n,k,source,dest
    INTEGER, DIMENSION(:), ALLOCATABLE :: Atom_node,counts
    LOGICAL, DIMENSION (:), ALLOCATABLE :: tags
    INTEGER :: ierr


    ALLOCATE(Tags(nato))
    Tags=.FALSE.
    ALLOCATE(Atom_node(nato))

    ALLOCATE(starts(nprocs),ends(nprocs))

    CALL MPI_ALLGATHER(nstart,1,MPI_INTEGER4,starts,1,MPI_INTEGER4&
         &,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLGATHER(nend,1,MPI_INTEGER4,ends,1,MPI_INTEGER4&
         &,MPI_COMM_WORLD,ierr)
    DO i=1,nprocs
       nst=starts(i)
       ned=ends(i)
       nst_m=atomp(nst)
       ned_m=atomp(ned)
       Atom_node(nst:ned)=i
    END DO


    ALLOCATE(list)
    NULLIFY(list%next)
    nfirst=>list

    CALL Get_Extra_Atoms(nstart,nend,lstrtch,2,lbnd_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,lbndg,3,lbndg_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,ltor,4,ltor_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,litr,4,litr_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,int14,2,int14_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,int13,2,int13_x,Tags,list,nfirst)
    CALL Get_Extra_Atoms(nstart,nend,ingrp,2,ingrp_x,Tags,list,nfirst)

 
    ALLOCATE(counts(nprocs))
    count=0
    counts=0
    list=>nfirst
    DO WHILE (ASSOCIATED(list%next))
       count=count+list%natoms
       DO i=1,list%natoms
          inode=Atom_node(list%atoms(i))
          counts(inode)=counts(inode)+1
       END DO
       list => list%next
    END DO

    ALLOCATE(Takes(nprocs),Gives(nprocs))
    DO i=1,nprocs
       IF(counts(i) /= 0) ALLOCATE(Takes(i)%atoms(counts(i)))
    END DO


    Takes(:)%natoms=0
    Gives(:)%natoms=0
    list=>nfirst
    DO WHILE (ASSOCIATED(list%next))
       DO i=1,list%natoms
          inode=Atom_node(list%atoms(i))
          Takes(inode)%natoms=Takes(inode)%natoms+1
          count=Takes(inode)%natoms
          Takes(inode)%atoms(count)=list%atoms(i)
       END DO
       list => list%next
    END DO

    CALL Get_partners

    CALL Pass_Takes_Gives

    CALL Make_Molecular(Tags)

    DEALLOCATE(counts,Tags,list,starts,ends,Atom_node)
  END SUBROUTINE Setup_Comm
  SUBROUTINE Get_Extra_Atoms(nstart,nend,intra,p1,intra_x,Tags,list,nfirst)
    IMPLICIT none 
    INTEGER :: p1
    INTEGER, DIMENSION(p1,*) :: intra
    INTEGER, DIMENSION(*) :: intra_x
    INTEGER :: nstart,nend
    TYPE(NODES), POINTER:: list,nfirst
    LOGICAL, DIMENSION (:) :: Tags

    INTEGER :: i,ii,n,count,ip1,k
    INTEGER, DIMENSION (:), ALLOCATABLE :: ats

    ALLOCATE(ats(p1))

    n=intra_x(1)
    DO ii=1,n
       i=IABS(intra_x(1+ii))
       count=0
       DO k=1,p1
          ip1=intra(k,i)
          IF((ip1 .LT. nstart .OR. ip1 .GT. nend) .AND. (.NOT.&
               & Tags(ip1))) THEN
             count=count+1
             ats(count)=ip1
             Tags(ip1)=.TRUE.
          END IF
       END DO
       IF(count .GT. 0) THEN
          list%natoms=count
          ALLOCATE(list%atoms(count))
          list%atoms(1:count)=ats(1:count)
          ALLOCATE(list%next)
          list=>list%next
          NULLIFY(list%next)
       END IF
    END DO
    DEALLOCATE(ats)
  END SUBROUTINE Get_Extra_Atoms
  SUBROUTINE Pass_Takes_Gives
    IMPLICIT none 
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ibuff_r,Ibuff_s
    INTEGER :: Dims_Ibuff,count,i,j,jj,m
    INTEGER :: ierr,mp,k,dest,source,max_dims,istatus(MPI_STATUS_SIZE)
    INTEGER, DIMENSION (:), ALLOCATABLE :: Dims

!!$
!!$------------------------------------------------------------------
!!$

    ALLOCATE(Dims(nprocs))
    Dims_Ibuff=SUM(Takes(1:nprocs)%natoms)+nprocs
    CALL MPI_ALLGATHER(Dims_Ibuff,1,MPI_INTEGER4,Dims,1,MPI_INTEGER4&
         &,MPI_COMM_WORLD,ierr)
    max_dims=0
    DO i=1,nprocs
       IF(max_dims < Dims(i)) max_dims=Dims(i)
    END DO
    Dims_Ibuff=max_dims
    ALLOCATE(ibuff_s(Dims_Ibuff))
    ALLOCATE(ibuff_r(Dims_Ibuff))

    CALL MPI_PACK_SIZE(Dims_Ibuff,MPI_REAL8,MPI_COMM_WORLD,Buffer_p&
         &,ierr)
    Buffer_p=Buffer_p+MPI_BSEND_OVERHEAD
    ALLOCATE(Buffer_MPI(Buffer_p))

    DO k=1,exchg
       dest=ipartners(2,k)
       source=ipartners(1,k)
       mp=0
       IF(node == dest) THEN
          CALL MPI_RECV(mp,1,MPI_INTEGER4,source,0,MPI_COMM_WORLD&
               &, istatus ,ierr)
          CALL MPI_RECV(Ibuff_r,mp,MPI_INTEGER4,source,1,MPI_COMM_WORLD&
               &, istatus ,ierr)
          Gives(source+1)%natoms=mp
          ALLOCATE(Gives(source+1)%atoms(mp))
          Gives(source+1)%atoms(1:mp)=Ibuff_r(1:mp)
       ELSE IF(node == source) THEN
          mp=Takes(dest+1)%natoms
          DO j=1,mp
             Ibuff_s(j)=Takes(dest+1)%atoms(j)
          END DO
          CALL Buffer_Attach
          CALL MPI_BSEND(mp,1,MPI_INTEGER4,dest,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BSEND(Ibuff_s,mp,MPI_INTEGER4,dest,1,MPI_COMM_WORLD,ierr)
          CALL Buffer_Detach
       END IF
    END DO

    DEALLOCATE(Buffer_MPI)
    DEALLOCATE(ibuff_s,Ibuff_r)
  END SUBROUTINE Pass_Takes_Gives
  SUBROUTINE Get_Partners
    IMPLICIT none 

    INTEGER :: i,j,ierr
    INTEGER, DIMENSION (:), ALLOCATABLE ::  pair
    INTEGER, DIMENSION (:,:), ALLOCATABLE ::  pairs
!!$
!!$------------------------------------------------------------------
!!$

    ALLOCATE(pair(nprocs),pairs(nprocs,nprocs))
    pair=0
    DO i=1,nprocs
       IF(Takes(i)%natoms /= 0) pair(i)=1
    END DO
    CALL MPI_ALLGATHER(pair,nprocs,MPI_INTEGER4,pairs,nprocs,MPI_INTEGER4&
         &,MPI_COMM_WORLD,ierr)
    
    exchg=0
    DO i=1,nprocs
       DO j=1,nprocs
          IF(pairs(j,i) == 1) THEN
             exchg=exchg+1
          END IF
       END DO
    END DO
    ALLOCATE(ipartners(2,exchg))
    exchg=0
    DO i=1,nprocs
       DO j=1,nprocs
          IF(pairs(j,i) == 1) THEN
             exchg=exchg+1
             ipartners(1,exchg)=i-1
             ipartners(2,exchg)=j-1
          END IF
       END DO
    END DO
    DEALLOCATE(pair,pairs)
  END SUBROUTINE Get_Partners
  SUBROUTINE Make_Molecular(Tags)
    IMPLICIT none 
    INTEGER :: count,m,i,j,nn,jj
    INTEGER, DIMENSION(:), ALLOCATABLE :: keep
    LOGICAL, DIMENSION (:) :: Tags

    ALLOCATE(Gives_m(nprocs))
    Gives_m(:)%natoms=0
    DO i=1,nprocs
       m=Gives(i)%natoms
       IF( m /= 0) THEN
          Tags=.FALSE.
          ALLOCATE(keep(m))
          count=0
          DO j=1,m
             nn=atomp(Gives(i)%atoms(j))
             IF(.NOT. Tags(nn)) THEN
                Tags(nn)=.TRUE.
                count=count+1
                keep(count)=nn
             END IF
          END DO
          Gives_m(i)%natoms=count
          IF(count /= 0) THEN
             ALLOCATE(Gives_m(i)%atoms(count))
             Gives_m(i)%atoms(1:count)=keep(1:count)
          END IF
          DEALLOCATE(keep)
       END IF
    END DO

    ALLOCATE(Takes_m(nprocs))
    Takes_m(:)%natoms=0
    DO i=1,nprocs
       m=Takes(i)%natoms
       IF( m /= 0) THEN
          Tags=.FALSE.
          ALLOCATE(keep(m))
          count=0
          DO j=1,m
             nn=atomp(Takes(i)%atoms(j))
             IF(.NOT. Tags(nn)) THEN
                Tags(nn)=.TRUE.
                count=count+1
                keep(count)=nn
             END IF
          END DO
          Takes_m(i)%natoms=count
          IF(count /= 0) THEN
             ALLOCATE(Takes_m(i)%atoms(count))
             Takes_m(i)%atoms(1:count)=keep(1:count)
          END IF
          DEALLOCATE(keep)
       END IF
    END DO
  END SUBROUTINE Make_Molecular
  SUBROUTINE Buffer_Attach
    IMPLICIT none 
    INTEGER :: ierr
    CALL MPI_BUFFER_ATTACH(Buffer_MPI,Buffer_p,ierr)
  END SUBROUTINE Buffer_Attach
  SUBROUTINE Buffer_Detach
    IMPLICIT none 
    INTEGER :: ierr
    CALL MPI_BUFFER_DETACH(Buffer_MPI,Buffer_p,ierr)
  END SUBROUTINE Buffer_Detach

#endif

END SUBROUTINE P_Comm_Intra
