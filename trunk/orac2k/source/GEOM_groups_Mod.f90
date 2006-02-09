MODULE GEOM_groups_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-02-09 14:32:36 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Feb  9 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
  
  USE Xerror_Mod
  USE INPUT_Mod, ONLY: Read_String, Parser, Parse_Numbers,err_open&
       &,err_end,err_unr,err_fnf,err_args
  USE GROUPS_Mod, ONLY: GR_groups=>groups, group, molecs

  LOGICAL, SAVE :: Geom_groups=.FALSE.
  INTEGER, SAVE :: n_write=1,k_write=0,natom
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: Properties
  TYPE Geom
     CHARACTER(20) :: type
     CHARACTER(20), DIMENSION(:), POINTER :: data
  END TYPE Geom
  TYPE(Geom), DIMENSION(:), ALLOCATABLE, SAVE :: Geoms
CONTAINS

  SUBROUTINE Init(ntap)

!!$======================== DECLARATIONS ================================*


    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: ntap

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: m,i,k
    CHARACTER(80) :: errmsg
    LOGICAL :: ok

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(.NOT. GR_groups) THEN
       errmsg='Cannot run GEOMS(&PROPERTIES) without GROUPS defined.'&
            &//' Abort/.'
       CALL abort_now(errmsg)
    ELSE
       natom=ntap
       ALLOCATE(Properties(SIZE(Geoms)))
       m=SIZE(Geoms)
       DO i=1,m
          ok=.FALSE.
          DO k=1,SIZE(molecs)
             IF(molecs(k)%type == Geoms(i)%data(1)) THEN
                ok=.TRUE.
             END IF
          END DO
          IF(.NOT. ok) THEN
             errmsg='First group type not found. Abort'
             CALL abort_now(errmsg)
          END IF
          ok=.FALSE.
          DO k=1,SIZE(molecs)
             IF(molecs(k)%type == Geoms(i)%data(2)) THEN
                ok=.TRUE.
             END IF
          END DO
          IF(.NOT. ok) THEN
             errmsg='Second group type not found. Abort'
             CALL abort_now(errmsg)
          END IF
       END DO
    END IF
  END SUBROUTINE Init
  SUBROUTINE Compute(xpa,ypa,zpa,mass,co)

!!$======================== DECLARATIONS ================================*

    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xpa(:),ypa(:),zpa(:),mass(:),co(3,3)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: m,i,ia,ib,k
    CHARACTER(80) :: errmsg
    INTEGER :: n_properties

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    m=SIZE(Geoms)
    DO i=1,m
       n_properties=i
       SELECT CASE(Geoms(i)%type)
       CASE('dist')
          ia=0
          ib=0
          DO k=1,SIZE(molecs)
             IF(molecs(k)%type == Geoms(i)%data(1)) THEN
                ia=k
             END IF
             IF(molecs(k)%type == Geoms(i)%data(2)) THEN
                ib=k
             END IF
          END DO
          CALL Get_Dist('N')
       CASE('dist_mass')
          ia=0
          ib=0
          DO k=1,SIZE(molecs)
             IF(molecs(k)%type == Geoms(i)%data(1)) THEN
                ia=k
             END IF
             IF(molecs(k)%type == Geoms(i)%data(2)) THEN
                ib=k
             END IF
          END DO
          CALL Get_Dist('Y')
       END SELECT
    END DO
  CONTAINS
    SUBROUTINE Get_dist(label)
      IMPLICIT NONE 
      CHARACTER(1) :: label
      REAL(8) :: xcma,ycma,zcma,xcmb,ycmb,zcmb,mass_b
      REAL(8) :: xa,ya,za,xc,yc,zc,rsp,PBC
      INTEGER :: ka,ja,kb,jb
      
      xcma=0.0D0
      ycma=0.0D0
      zcma=0.0D0
      DO ka=1,molecs(ia)%n
         ja=molecs(ia)%index(ka)
         mass_b=1.0D0
         IF(label == 'Y') THEN
            mass_b=1.0D0/mass(ja)
         END IF
         xcma=xcma+mass_b*xpa(ja)
         ycma=ycma+mass_b*ypa(ja)
         zcma=zcma+mass_b*zpa(ja)
      END DO
      xcma=xcma/DBLE(molecs(ia)%n)
      ycma=ycma/DBLE(molecs(ia)%n)
      zcma=zcma/DBLE(molecs(ia)%n)
      
      xcmb=0.0D0
      ycmb=0.0D0
      zcmb=0.0D0
      DO kb=1,molecs(ib)%n
         jb=molecs(ib)%index(kb)
         mass_b=1.0D0
         IF(label == 'Y') THEN
            mass_b=1.0D0/mass(jb)
         END IF
         xcmb=xcmb+mass_b*xpa(jb)
         ycmb=ycmb+mass_b*ypa(jb)
         zcmb=zcmb+mass_b*zpa(jb)
      END DO
      xcmb=xcmb/DBLE(molecs(ib)%n)
      ycmb=ycmb/DBLE(molecs(ib)%n)
      zcmb=zcmb/DBLE(molecs(ib)%n)

      xa=xcmb-xcma
      ya=ycmb-ycma
      za=zcmb-zcma
      xa=xa-2.0D0*PBC(xa)
      ya=ya-2.0D0*PBC(ya)
      za=za-2.0D0*PBC(za)
      xc=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
      yc=           co(2,2)*ya+co(2,3)*za
      zc=                      co(3,3)*za
      rsp=SQRT(xc*xc+yc*yc+zc*zc)
      Properties(n_properties)=rsp
    END SUBROUTINE Get_dist
  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)

!!$======================== DECLARATIONS ================================*

    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: fstep
    INTEGER :: m
    CHARACTER(80) :: errmsg
    INTEGER, SAVE :: Time_of_calls=0

!!$------------------------- LOCAL VARIABLES ----------------------------*

    m=SIZE(Properties)
    IF(m > 15) THEN
       errmsg='Cannot compute more than 15 group geometries at once.'&
            &//' Abort.'
       CALL abort_now(errmsg)
    END IF
    IF(Time_of_calls == 0) THEN
       WRITE(k_write,100) 
    END IF
    Time_of_calls=Time_of_calls+1
    WRITE(k_write,'(f12.3,2x,15f15.6)') fstep,Properties(1:m)
    
100 FORMAT('**************************************************************'&
         &/'*                                                            *'&
         &/'*        Writing group properties for the systems            *'&
         &/'*         Format: Timestep  prop1 prop2 prop3 etc.           *'&
         &/'*                                                            *'&
         &/'**************************************************************'/)

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  END SUBROUTINE Write_it
  REAL(8) FUNCTION PBC(x)
    REAL(8) :: x
    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  INCLUDE 'GEOM_groups_Readit.f90'
END MODULE GEOM_groups_Mod
