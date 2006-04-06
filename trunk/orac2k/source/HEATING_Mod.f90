MODULE HEATING_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-04-03 12:06:46 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Mar 31 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  INTEGER, SAVE :: ndgree_Tot,ndgree_slv,ndgree_slt,ntap,nstart,nend&
       &,node
  REAL(8), SAVE :: efact,T_slt,gascon=8.3143d0,t,dtemp
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ss_index

!!$======================== DECLARATIONS ================================*

CONTAINS
  SUBROUTINE Init(Ta,dtempa,natom,nato_slv,nmol,cnstpp,cnstpp_slv&
       &,ss_indexa,efacta,nstarta,nenda,nodea)
    IMPLICIT NONE 
    INTEGER :: natom,nato_slv,nmol,cnstpp,cnstpp_slv,ss_indexa(*)&
         &,nstarta,nenda,nodea
    REAL(8) :: efacta,Ta,dtempa
    
    T=Ta
    dtemp=dtempa
    nstart=nstarta
    nend=nenda
    node=nodea
    ntap=natom
    ALLOCATE(ss_index(ntap))

    ss_index=ss_indexa(1:ntap)

    ndgree_Tot=3*ntap-cnstpp
    ndgree_slv=3*nato_slv*nmol-cnstpp_slv
    ndgree_slt=ndgree_Tot-ndgree_slv

  END SUBROUTINE Init
  SUBROUTINE Set_slv(ucek,vpx,vpy,vpz)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: vpx(*),vpy(*),vpz(*),ucek

!!$------------------------- LOCAL VARIABLES ----------------------------*
    
    INTEGER :: i,type
    REAL(8) :: temp_slv,scale

!!$----------------------------------------------------------------------*

    temp_slv=2.0D0*ucek/(gascon*ndgree_slv)
    IF(DABS(T-temp_slv) > dtemp) THEN
       IF(node == 0) WRITE(*,'(''Solvent Temperature at ---->'',f10.3&
         &,'' Reset to ----> '',f10.3)') temp_slv,T
       scale=SQRT(T/temp_slv)
       DO i=nstart,nend
          type=ss_index(i)
          IF(type == 2) THEN
             vpx(i)=vpx(i)*scale
             vpy(i)=vpy(i)*scale
             vpz(i)=vpz(i)*scale
          END IF
       END DO
    END IF
  END SUBROUTINE Set_slv
  SUBROUTINE Set_slt(T_slt,pucek,vpx,vpy,vpz)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: T_slt,vpx(*),vpy(*),vpz(*),pucek

!!$------------------------- LOCAL VARIABLES ----------------------------*
    
    INTEGER :: i,type
    REAL(8) :: temp_slt,scale

!!$----------------------------------------------------------------------*

    temp_slt=2.0D0*pucek/(gascon*ndgree_slt)
    IF(node == 0) WRITE(*,'(''Solute Temperature at ---->'',f10.3&
         &,'' Set Temperature to ---->'',f10.3)') temp_slt,T_slt
    scale=SQRT(T_slt/temp_slt)
    DO i=nstart,nend
       type=ss_index(i)
       IF(type == 1) THEN
          vpx(i)=vpx(i)*scale
          vpy(i)=vpy(i)*scale
          vpz(i)=vpz(i)*scale
       END IF
    END DO
  END SUBROUTINE Set_slt

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE HEATING_Mod
