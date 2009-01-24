!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
SUBROUTINE Imph
  INTEGER, SAVE :: MyCalls=0
  REAL(8), SAVE :: quasi_zero=1.0D-12

  INTEGER :: i,l1,l2,l3,l4,ntph,tmp,type,ntors,nWeight
  REAL(8) ::  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32&
       &,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32&
       &,rsp43,pota,potb,potc
  REAL(8) ::  cb1,cb2,cb3,sb1,sb2,sb3,aux,aux1,aux2,auxa,auxb,auxc&
       &,dbx11,dbx31,dbx12,dbx22,dbx32,dbx13,dbx23,dbx33,dbx24&
       &,dbx34,dby11,dby31,dby12,dby22,dby32,dby13,dby23,dby33&
       &,dby24,dby34,dbz11,dbz31,dbz12,dbz22,dbz32,dbz13,dbz23&
       &,dbz33,dbz24,dbz34,uux1,uux2,uux3,uux4,uuy1,uuy2,uuy3&
       &,uuy4,uuz1,uuz2,uuz3,uuz4,aux3,bb,soa
  REAL(8) ::  coa,enrg,qforce,sig,utors(2),Weight



  IF(Calls == 0) THEN
     Calls=Calls+1
     Conv_Fact=1000.0D0*4.184/(unite*avogad)
  END IF

  ntors=SIZE(Indx_Imph,2)


  utors=0.0_8
  DO i=1,ntors
     l1=Indx_Imph(1,i)
     l2=Indx_Imph(2,i)
     l3=Indx_Imph(3,i)
     l4=Indx_Imph(4,i)
     pota=Param_Imph(i) % pot(1)*Conv_Fact
     potc=Param_Imph(i) % pot(3)*Degree_To_Rad
     nWeight=SIZE(Param_Imph(i) % pot)
     Weight=Param_Imph(i) % pot(nWeight)

     type=Slv(l1)


     xr1=xp0(l1)
     yr1=yp0(l1)
     zr1=zp0(l1)
     xr2=xp0(l2)
     yr2=yp0(l2)
     zr2=zp0(l2)
     xr3=xp0(l3)
     yr3=yp0(l3)
     zr3=zp0(l3)
     xr4=xp0(l4)
     yr4=yp0(l4)
     zr4=zp0(l4)
     x21=xr2-xr1
     y21=yr2-yr1
     z21=zr2-zr1
     x32=xr3-xr2
     y32=yr3-yr2
     z32=zr3-zr2
     x43=xr4-xr3
     y43=yr4-yr3
     z43=zr4-zr3
     rsq21=x21**2+y21**2+z21**2
     rsq32=x32**2+y32**2+z32**2
     rsq43=x43**2+y43**2+z43**2
     rsp21=DSQRT(rsq21)
     rsp32=DSQRT(rsq32)
     rsp43=DSQRT(rsq43)
     cb1=(x21*x32+y21*y32+z21*z32)/(rsp21*rsp32)
     cb2=(x43*x32+y43*y32+z43*z32)/(rsp43*rsp32)
     cb3=(x21*x43+y21*y43+z21*z43)/(rsp21*rsp43)
     sb1=DSQRT(DABS(1.0d0-cb1**2))
     sb2=DSQRT(DABS(1.0d0-cb2**2))
     sb3=DSQRT(DABS(1.0d0-cb3**2))
     IF(sb1 .EQ. 0.0D0) sb1=quasi_zero
     IF(sb2 .EQ. 0.0D0) sb2=quasi_zero
     IF(sb3 .EQ. 0.0D0) sb3=quasi_zero
     aux=sb1*sb2
     coa=(cb1*cb2-cb3)/aux
     IF(coa.GT.1.AND.coa.LT.1.000001) coa=1.0d0
     IF(coa.LT.-1.0d0.AND.coa.GT.-1.000001) coa=-1.0d0
     aux1=cb1*sb2
     aux2=sb1*cb2
     auxa=-(aux1*coa+aux2)/aux
     auxb=-(aux2*coa+aux1)/aux
     auxc=sb3/aux
     aux1=rsp21*rsp32*sb1
     aux2=rsp32*rsp43*sb2
     aux3=rsp21*rsp43*sb3
     dbx11= x32/aux1-cb1*x21/(rsq21*sb1)
     dbx13=-x21/aux1+cb1*x32/(rsq32*sb1)
     dbx12=-dbx13-dbx11
     dbx22= x43/aux2-cb2*x32/(rsq32*sb2)
     dbx24=-x32/aux2+cb2*x43/(rsq43*sb2)
     dbx23=-dbx22-dbx24
     dbx31=x43/aux3-cb3*x21/(rsq21*sb3)
     dbx33=x21/aux3-cb3*x43/(rsq43*sb3)
     dbx32=-dbx31
     dbx34=-dbx33
     dby11= y32/aux1-cb1*y21/(rsq21*sb1)
     dby13=-y21/aux1+cb1*y32/(rsq32*sb1)
     dby12=-dby13-dby11
     dby22= y43/aux2-cb2*y32/(rsq32*sb2)
     dby24=-y32/aux2+cb2*y43/(rsq43*sb2)
     dby23=-dby22-dby24
     dby31=y43/aux3-cb3*y21/(rsq21*sb3)
     dby33=y21/aux3-cb3*y43/(rsq43*sb3)
     dby32=-dby31
     dby34=-dby33
     dbz11= z32/aux1-cb1*z21/(rsq21*sb1)
     dbz13=-z21/aux1+cb1*z32/(rsq32*sb1)
     dbz12=-dbz13-dbz11
     dbz22= z43/aux2-cb2*z32/(rsq32*sb2)
     dbz24=-z32/aux2+cb2*z43/(rsq43*sb2)
     dbz23=-dbz22-dbz24
     dbz31=z43/aux3-cb3*z21/(rsq21*sb3)
     dbz33=z21/aux3-cb3*z43/(rsq43*sb3)
     dbz32=-dbz31
     dbz34=-dbz33
     
     bb=DACOS(coa)
     IF(bb.LT.1.0D-06) bb=1.0D-06
     soa=SIN(bb)
     qforce=(2.0d0*pota*(bb-potc))/soa
     utors(type)=utors(type)+Weight*pota*(bb-potc)**2
     uux1=auxa*dbx11+auxc*dbx31
     uux2=auxa*dbx12+auxb*dbx22+auxc*dbx32
     uux3=auxa*dbx13+auxb*dbx23+auxc*dbx33
     uux4=auxb*dbx24+auxc*dbx34
     uuy1=auxa*dby11+auxc*dby31
     uuy2=auxa*dby12+auxb*dby22+auxc*dby32
     uuy3=auxa*dby13+auxb*dby23+auxc*dby33
     uuy4=auxb*dby24+auxc*dby34
     uuz1=auxa*dbz11+auxc*dbz31
     uuz2=auxa*dbz12+auxb*dbz22+auxc*dbz32
     uuz3=auxa*dbz13+auxb*dbz23+auxc*dbz33
     uuz4=auxb*dbz24+auxc*dbz34
     fpx(l1)=fpx(l1)+qforce*uux1
     fpx(l2)=fpx(l2)+qforce*uux2
     fpx(l3)=fpx(l3)+qforce*uux3
     fpx(l4)=fpx(l4)+qforce*uux4
     fpy(l1)=fpy(l1)+qforce*uuy1
     fpy(l2)=fpy(l2)+qforce*uuy2
     fpy(l3)=fpy(l3)+qforce*uuy3
     fpy(l4)=fpy(l4)+qforce*uuy4
     fpz(l1)=fpz(l1)+qforce*uuz1
     fpz(l2)=fpz(l2)+qforce*uuz2
     fpz(l3)=fpz(l3)+qforce*uuz3
     fpz(l4)=fpz(l4)+qforce*uuz4
  END DO
  Uitors_slt=utors(1)
  Uitors_slv=utors(2)
  CALL EN_Imph_(Uitors_Slv,Uitors_Slt,0.0D0)

END SUBROUTINE Imph
