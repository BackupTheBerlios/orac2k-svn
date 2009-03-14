SUBROUTINE Get_vp_without_vc(vpx,vpy,vpz,oc,nstart,nend,nprot,protl&
     &,massa,massb,vpcmx,vpcmy,vpcmz)

!!$***********************************************************************
!!$   Time-stamp: <04/11/25 11:38:50 marchi>                           *
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

  REAL(8) ::  vpx(*),vpy(*),vpz(*),vpcmx(*),vpcmy(*),vpcmz(*),oc(3,3)&
       &,massa(*),massb(*)
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
     vpcmx(j)=0.0D0
     vpcmy(j)=0.0D0
     vpcmz(j)=0.0D0
     m=protl(1+count)
     DO i=1,m
        i1=protl(count+1+i)
        vpcmx(j)=vpcmx(j)+massa(i1)*vpx(i1)
        vpcmy(j)=vpcmy(j)+massa(i1)*vpy(i1)
        vpcmz(j)=vpcmz(j)+massa(i1)*vpz(i1)
     END DO
     vpcmx(j)=vpcmx(j)/massb(j)
     vpcmy(j)=vpcmy(j)/massb(j)
     vpcmz(j)=vpcmz(j)/massb(j)
     DO i=1,m
        i1=protl(count+1+i)
        vpx(i1)=vpx(i1)-vpcmx(j)
        vpy(i1)=vpy(i1)-vpcmy(j)
        vpz(i1)=vpz(i1)-vpcmz(j)
     END DO
     
     xc=oc(1,1)*vpcmx(j)+oc(1,2)*vpcmy(j)+oc(1,3)*vpcmz(j)
     yc=oc(2,1)*vpcmx(j)+oc(2,2)*vpcmy(j)+oc(2,3)*vpcmz(j)
     zc=oc(3,1)*vpcmx(j)+oc(3,2)*vpcmy(j)+oc(3,3)*vpcmz(j)
     vpcmx(j)=xc
     vpcmy(j)=yc
     vpcmz(j)=zc
     count=count+m+1
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Get_vp_without_vc
