
               muj_x=dipole(1,j)
               muj_y=dipole(2,j)
               muj_z=dipole(3,j)
               cgj=charge(j)
               cgij=cgi*cgj
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
               
               Gij0=cgi*cgj
               Gij1=dotij+cgi*dotjr-cgj*dotir
               Gij2=-dotjr*dotir
               term=Gij0*B1+Gij1*B2+Gij2*B3
               
               echarge=echarge+B0*Gij0
               epol=epol+Gij1*B1+Gij2*B2
               ucoula=ucoula+B0*Gij0+Gij1*B1+Gij2*B2
               
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

c$$$========================================================================
c$$$--- Do Thole correction 
c$$$========================================================================
               
               arg = rsp/plrzbij(lij)
               fact = 1.0d0 + arg + 0.50d0*arg**2
               cexp = exp(-arg)
               D1=fact*cexp
               D2=(fact+arg**3/6.0d0)*cexp
               
               BB1=B1
               BB2=B2
               BB3=B3

               B1=B1*D1
               B2=B2*D2

               Gij11 =  dotij

               
               termj=-dotjr*B2
               termi=-dotir*B2

               dphii_dx = termj*xc + B1*muj_x
               dphii_dy = termj*yc + B1*muj_y
               dphii_dz = termj*zc + B1*muj_z
               dphij_dx = termi*xc + B1*mui_x
               dphij_dy = termi*yc + B1*mui_y
               dphij_dz = termi*zc + B1*mui_z

               Edx(i1)=Edx(i1)+dphii_dx
               Edy(i1)=Edy(i1)+dphii_dy
               Edz(i1)=Edz(i1)+dphii_dz
               Edx(j)=Edx(j)+dphij_dx
               Edy(j)=Edy(j)+dphij_dy
               Edz(j)=Edz(j)+dphij_dz



               Uddd=Uddd-Gij11*B1-Gij2*B2

               fact=arg+arg**2
               d_D1=-(fact*cexp-arg*D1)*rsqi
               fact=fact+0.5D0*arg**3
               d_D2=-(fact*cexp-arg*D2)*rsqi

               d_B3=BB3*D2+BB2*d_D2
               d_B2=BB2*D1+BB1*d_D1


               term=Gij11*d_B2+Gij2*d_B3


               dphii_dx = term*xc-termj*mui_x-termi*muj_x
               dphii_dy = term*yc-termj*mui_y-termi*muj_y
               dphii_dz = term*zc-termj*mui_z-termi*muj_z
               
               fpx(i1)=fpx(i1)-dphii_dx
               fpy(i1)=fpy(i1)-dphii_dy
               fpz(i1)=fpz(i1)-dphii_dz
               fpx(j)=fpx(j)+dphii_dx
               fpy(j)=fpy(j)+dphii_dy
               fpz(j)=fpz(j)+dphii_dz


c$$$               IF(rsp .LT. 3.0D0) THEN
c$$$             
c$$$                  Gij11 =  dotij
c$$$
c$$$                  arg = rsp/plrzbij(lij)
c$$$                  fact = 1.0d0 + arg + 0.50d0*arg**2
c$$$                  cexp = exp(-arg)
c$$$                  
c$$$                  B1_0=(1.0D0/rsp)*rsqi
c$$$                  B2_0=3.0D0*B1_0*rsqi
c$$$                  
c$$$                  B1_0=B1_0*fact*cexp
c$$$                  B2_0=B2_0*(fact+arg**3/6.0d0)*cexp
c$$$
c$$$                  termj=-dotjr*B2_0
c$$$                  termi=-dotir*B2_0
c$$$                  dphii_dx = termj*xc + B1_0*muj_x
c$$$                  dphii_dy = termj*yc + B1_0*muj_y
c$$$                  dphii_dz = termj*zc + B1_0*muj_z
c$$$                  dphij_dx = termi*xc + B1_0*mui_x
c$$$                  dphij_dy = termi*yc + B1_0*mui_y
c$$$                  dphij_dz = termi*zc + B1_0*mui_z
c$$$                  Edx(i1)=Edx(i1)+dphii_dx
c$$$                  Edy(i1)=Edy(i1)+dphii_dy
c$$$                  Edz(i1)=Edz(i1)+dphii_dz
c$$$                  Edx(j)=Edx(j)+dphij_dx
c$$$                  Edy(j)=Edy(j)+dphij_dy
c$$$                  Edz(j)=Edz(j)+dphij_dz
c$$$                  Uddd=Uddd-Gij11*B1_0-Gij2*B2_0
c$$$
c$$$               END IF
