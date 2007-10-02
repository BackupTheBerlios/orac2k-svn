SUBROUTINE Get_vp_with_vc(vpx,vpy,vpz,co,nstart,nend,nprot,protl&
     &,vpcmx,vpcmy,vpcmz)

!!$***********************************************************************
!!$   Time-stamp: <04/11/19 15:15:07 marchi>                           *
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

  REAL(8) ::  vpx(*),vpy(*),vpz(*),vpcmx(*),vpcmy(*),vpcmz(*),co(3,3)
  INTEGER :: nprot,protl(*),nstart,nend

!!$----------------------- VARIABLES IN COMMON --------------------------*

  
!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: i,j,count,m,i1
  REAL(8) :: xc,yc,zc

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  count=0
  DO j=1,nstart-1
     m=protl(1+count)
     count=count+m+1
  END DO
  DO j=nstart,nend
     xc=co(1,1)*vpcmx(j)+co(1,2)*vpcmy(j)+co(1,3)*vpcmz(j)
     yc=co(2,1)*vpcmx(j)+co(2,2)*vpcmy(j)+co(2,3)*vpcmz(j)
     zc=co(3,1)*vpcmx(j)+co(3,2)*vpcmy(j)+co(3,3)*vpcmz(j)
     m=protl(1+count)
     DO i=1,m
        i1=protl(count+1+i)
        vpx(i1)=vpx(i1)+xc
        vpy(i1)=vpy(i1)+yc
        vpz(i1)=vpz(i1)+zc
     END DO
     count=count+m+1
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Get_vp_with_vc
