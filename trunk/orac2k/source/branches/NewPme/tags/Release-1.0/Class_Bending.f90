!!$***********************************************************************
!!$   Time-stamp: <2005-03-03 21:36:53 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 18 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
MODULE Class_Bending
  USE Class_Connect, Connect_Init=> Init, Connect_Destroy=>Destroy,&
       & Connect_Print=> Print 
  TYPE Bending
     INTEGER, DIMENSION(:,:), POINTER :: atm=>NULL()
  END TYPE Bending
CONTAINS
  FUNCTION Init(table) RESULT (out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    TYPE(Connect), DIMENSION(:), POINTER :: table
    TYPE(Bending) :: out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: dum(2)
    TYPE NODE
       INTEGER, DIMENSION(3) :: bnd
       TYPE(NODE), POINTER :: next
    END TYPE NODE
    TYPE(NODE), POINTER :: list,lista,nfirst
    INTEGER :: nato,count,i,m,n,a,b,ia,ib,ijk
    INTEGER, DIMENSION (:), ALLOCATABLE :: ind,ind1
    INTEGER, DIMENSION (:,:), ALLOCATABLE :: bnd_a

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    dum=Connect_Print()
    nato=dum(1)
    ALLOCATE(list)
    NULLIFY(list%next)
    nfirst=>list
    count=0
    DO i=1,nato
       m=table(i) % m
       DO a=1,m
          ia=table(i) % cnt(a)
          n=table(ia) % m
          DO b=1,n
             ib=table(ia) % cnt(b)
             IF((i /= ib .AND. ib /= ia) .AND. i < ib) THEN
                count=count+1
                list % bnd(1)=i
                list % bnd(2)=ia
                list % bnd(3)=ib
                ALLOCATE(list%next)
                list=>list % next
                NULLIFY(list % next)
             END IF
          END DO
       END DO
    END DO
    WRITE(*,*) 'count =',count
    ALLOCATE(bnd_a(3,count),ind(count),ind1(count))
    list => nfirst
    count=0
    DO WHILE (ASSOCIATED(list%next))
       count=count+1
       ind(count)=list % bnd(1)
       ind1(count)=count
       bnd_a(1,count)=list % bnd(1)
       bnd_a(2,count)=list % bnd(2)
       bnd_a(3,count)=list % bnd(3)
       lista=>list
       list => list%next
       DEALLOCATE(lista)
    END DO
    WRITE(*,*) 'count =',count
    CALL iindexx(count,ind1,ind)
    m=count
    count=0
    bnd_a(1,1)=bnd_a(1,ind1(1))
    bnd_a(2,1)=bnd_a(2,ind1(1))
    bnd_a(3,1)=bnd_a(3,ind1(1))
    DO i=2,m
       IF(ind(ind1(i)) == ind(ind1(i-1)) ) THEN

          count=count+1
          bnd_a(1,count)=bnd_a(1,ind1(i))
          bnd_a(2,count)=bnd_a(2,ind1(i))
          bnd_a(3,count)=bnd_a(3,ind1(i))
       END IF
    END DO
    ALLOCATE(out % atm (3,count))
    out % atm(1,:)=bnd_a(1,1:count)
    out % atm(2,:)=bnd_a(2,1:count)
    out % atm(3,:)=bnd_a(3,1:count)
    DEALLOCATE(bnd_a,ind,ind1)
    WRITE(*,*) 'Exit '

  END FUNCTION Init
  
END MODULE Class_Bending
