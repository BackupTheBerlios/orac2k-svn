SUBROUTINE GetKineticTest(nstart_2,nend_2,nstart_cm,nend_cm,UTotKin,vpx &
     ,vpy,vpz,vcax,vcay,vcaz,massp,mass,co)

!!$***********************************************************************
!!$   Time-stamp: <01/03/07 12:31:19 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Mar  7 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program orac ----*


!!$======================== DECLARATIONS ================================*

 IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

 INTEGER(4) :: nstart_2,nend_2,nstart_cm,nend_cm
 REAL(8) :: UTotKin,massp(*),co(3,3),vpx(*),vpy(*),vpz(*),vcax(*) &
      ,vcay(*),vcaz(*),mass(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

 INTEGER :: i
 REAL(8) :: TotKin,TotKinCM,vxyz(3)

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


 TotKin=0.0D0
 TotKinCM=0.0D0
 DO i=nstart_2,nend_2
    TotKin=TotKin+(0.50D+0)*mass(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
 END DO
 
 DO i=nstart_cm,nend_cm
    vxyz(1)=co(1,1)*vcax(i)+co(1,2)*vcay(i)+co(1,3)*vcaz(i)
    vxyz(2)=co(2,1)*vcax(i)+co(2,2)*vcay(i)+co(2,3)*vcaz(i)
    vxyz(3)=co(3,1)*vcax(i)+co(3,2)*vcay(i)+co(3,3)*vcaz(i)
    TotKinCM=TotKinCM+0.5D0*massp(i)*(vxyz(1)**2+vxyz(2)**2+vxyz(3)**2)
 END DO

 UTotKin=(TotKinCM+TotKin)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE GetKineticTest
