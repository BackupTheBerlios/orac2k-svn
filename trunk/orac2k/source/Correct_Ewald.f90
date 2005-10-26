SUBROUTINE Correct_Ewald(alphal,charge,dipole,xp0,yp0,zp0,list,nlist&
     &,sp,Utotal,Edx,Edy,Edz,fpx,fpy,fpz)
  
!!$======================== DECLARATIONS ================================*
  
  IMPLICIT NONE 
  
!!$----------------------------- ARGUMENTS ------------------------------*
  
  INTEGER :: nlist,list(2,*),sp(*)
  REAL(8) :: alphal,fudge,Energy,charge(*),dipole(3,*),xp0(*)&
       &,yp0(*),zp0(*),Edx(*),Edy(*),Edz(*),fpx(*),fpy(*),fpz(*)&
       &,Utotal 
  
!!$----------------------- VARIABLES IN COMMON --------------------------*
  
  
!!$------------------------- LOCAL VARIABLES ----------------------------*
  
  INTEGER :: i,i1,j,mma,ii,iam
  REAL(8) :: cgi,cgj,cgij,mui_x,mui_y,mui_z,muj_x,muj_y,muj_z&
       &,dotir,dotjr,dotij,xpi,ypi,zpi,xpj,ypj,zpj,Udd,Ucc&
       &,termi,termj,xc,yc,zc,rsq,rsp,rsqi,aux2,facteur,rspi,ucoula
  REAL(8) :: qt,expcst,erfcst,rspqi,alphar,furpar,twrtpi,fac,switch&
       &,d_switch,aux1,arg,fact,cexp,B0,B1,B2,B3,Gij0,Gij1&
       &,Gij2,B0_0,B1_0,B2_0,pi,B3_0,dphij_dx,dphij_dy,dphij_dz&
       &,dphii_dx,dphii_dy,dphii_dz,term
  REAL(8), SAVE :: a1=0.2548296d0,a2=-0.28449674d0,a3=1.4214137d0&
       &,a4=-1.453152d0,a5=1.0614054d0,qp=0.3275911d0
  
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
  
  Udd=0.0D0
  Ucc=0.0D0
  Utotal=0.0D0
  pi=4.0D0*DATAN(1.0D0)
  twrtpi=1.0d0/DSQRT(pi)/alphal
  fac=2.0D0*alphal*alphal
  mma=sp(1)
  DO ii=1,mma
     iam=sp(ii+1)
     i=IABS(iam)
     facteur=FACTOR(iam)

     i1=list(1,i)
     j=list(2,i)
     cgi=charge(i1)
     cgj=charge(j)
     cgij=cgi*cgj

     mui_x=dipole(1,i1)
     mui_y=dipole(2,i1)
     mui_z=dipole(3,i1)
     muj_x=dipole(1,j)
     muj_y=dipole(2,j)
     muj_z=dipole(3,j)

     xpi=xp0(i1)
     ypi=yp0(i1)
     zpi=zp0(i1)
     xpj=xp0(j)
     ypj=yp0(j)
     zpj=zp0(j)
     
     
     xc=xpi-xpj
     yc=ypi-ypj
     zc=zpi-zpj
     rsq=xc**2+yc**2+zc**2

     dotir=mui_x*xc+mui_y*yc+mui_z*zc
     dotjr=muj_x*xc+muj_y*yc+muj_z*zc
     dotij=mui_x*muj_x+mui_y*muj_y+mui_z*muj_z
               
     rsqi=1.0d0/rsq
     rsp=DSQRT(rsq)
     rspi=1.0D0/rsp
     alphar=alphal*rsp
                        
     qt=1.0d0/(1.0d0+qp*alphar)
     expcst=DEXP(-alphar*alphar)
     erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
     
     switch=erfcst
     d_switch=twrtpi*expcst
     
     B0=switch*rspi
     fact=d_switch*fac
     B1=(B0+fact)*rsqi
     fact=fact*fac
     B2=(3.0D0*B1+fact)*rsqi
     fact=fact*fac
     B3=(5.0D0*B2+fact)*rsqi
     
     B0_0=rspi
     B1_0=B0_0*rsqi
     B2_0=3.0D0*B1_0*rsqi
     B3_0=5.0D0*B2_0*rsqi

     B0=B0-B0_0
     B1=B1-B1_0
     B2=B2-B2_0
     B3=B3-B3_0
     
     Gij0=cgi*cgj
     Gij1=dotij+cgi*dotjr-cgj*dotir
     Gij2=-dotjr*dotir

     term=facteur*(Gij1*B2+Gij2*B3)
     
     Utotal=Utotal+facteur*(Gij1*B1+Gij2*B2)
     
     termj=-dotjr*B2-cgj*B1
     termi=-dotir*B2+cgi*B1
     dphii_dx = termj*xc + B1*muj_x
     dphii_dy = termj*yc + B1*muj_y
     dphii_dz = termj*zc + B1*muj_z
     dphij_dx = termi*xc + B1*mui_x
     dphij_dy = termi*yc + B1*mui_y
     dphij_dz = termi*zc + B1*mui_z
     
     Edx(i1)=Edx(i1)-dphii_dx
     Edy(i1)=Edy(i1)-dphii_dy
     Edz(i1)=Edz(i1)-dphii_dz
     Edx(j)=Edx(j)-dphij_dx
     Edy(j)=Edy(j)-dphij_dy
     Edz(j)=Edz(j)-dphij_dz
     
     dphii_dx = term*xc-termj*mui_x-termi*muj_x
     dphii_dy = term*yc-termj*mui_y-termi*muj_y
     dphii_dz = term*zc-termj*mui_z-termi*muj_z
     
     fpx(i1)=fpx(i1)+dphii_dx
     fpy(i1)=fpy(i1)+dphii_dy
     fpz(i1)=fpz(i1)+dphii_dz
     fpx(j)=fpx(j)-dphii_dx
     fpy(j)=fpy(j)-dphii_dy
     fpz(j)=fpz(j)-dphii_dz

  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
!!$       B3=(5.0D0*B2+fact)*rsqi
  CONTAINS
    FUNCTION FACTOR(i)
      IMPLICIT NONE 
      REAL(8) :: FACTOR
      INTEGER :: i
      REAL(8), SAVE :: half=0.5D0
      FACTOR=(1.0D0+DBLE(i/IABS(i)))*half
    END FUNCTION FACTOR
END SUBROUTINE Correct_Ewald
