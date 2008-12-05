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
SUBROUTINE Angles(fpx,fpy,fpz)
  REAL(8) :: fpx(:),fpy(:),fpz(:)
  REAL(8), POINTER :: xp0(:),yp0(:),zp0(:)
  REAL(8) ::  xr1,xr2,xr3,yr1,yr2,yr3,zr1,zr2,zr3,x12,x32,y12,y32&
       &,z12,z32,rs12,rs32,uux1,uux2,uux3,uuy1,uuy2,uuy3,uuz1,uuz2&
       &,uuz3,xr31,yr31,zr31,rs31,rsp31,x31,y31,z31,rsp12,rsp32,k12&
       &,r1,r2 
  REAL(8) :: dcc2,cb,sb,bb,qforce,pforce,pi,ubend(2),pota,potb&
       &,potc,potd
  INTEGER i,la,lb,lc,type,nbend
  
  
  
  IF(Calls == 0) THEN
     Calls=Calls+1
     Conv_Fact=1000.0D0*4.184/(unite*avogad)
  END IF
  
  
  xp0=>Atoms(:) % x
  yp0=>Atoms(:) % y
  zp0=>Atoms(:) % z
  
  nbend=SIZE(Indx_Angles,2)
  ubend=0.0_8
  DO i=1,nbend
     la=Indx_Angles(1,i)
     lb=Indx_Angles(2,i)
     lc=Indx_Angles(3,i)
     pota=Param_Angles(i) % pot(2)*Degree_To_Rad
     potb=Param_Angles(i) % pot(1)*Conv_Fact
     type=Atoms(la) % Id_Slv
     
     xr1=xp0(la)
     yr1=yp0(la)
     zr1=zp0(la)
     xr2=xp0(lb)
     yr2=yp0(lb)
     zr2=zp0(lb)
     xr3=xp0(lc)
     yr3=yp0(lc)
     zr3=zp0(lc)
     x12=xr1-xr2
     y12=yr1-yr2
     z12=zr1-zr2
     x32=xr3-xr2
     y32=yr3-yr2
     z32=zr3-zr2
     rs12=x12**2+y12**2+z12**2
     rs32=x32**2+y32**2+z32**2
     
     dcc2=DSQRT(rs12*rs32)
     cb=(x12*x32+y12*y32+z12*z32)/dcc2
     sb=DSQRT(1.0d0-cb**2)
     bb=DACOS(cb)
!!$c--      switch to singularity free potential for 
!!$c--      linear bending V = 2*K*(cos(theta) + 1)  
     IF(ABS(pota-pi) < 0.01) THEN 
        qforce=-2.0d0*potb
     ELSE
        qforce=2.0d0*potb*(bb-pota)/sb
     END IF
     uux1=x32/dcc2-cb*x12/rs12
     uux3=x12/dcc2-cb*x32/rs32
     uux2=-uux1-uux3
     uuy1=y32/dcc2-cb*y12/rs12
     uuy3=y12/dcc2-cb*y32/rs32
     uuy2=-uuy1-uuy3
     uuz1=z32/dcc2-cb*z12/rs12
     uuz3=z12/dcc2-cb*z32/rs32
     uuz2=-uuz1-uuz3
     fpx(la)=fpx(la)+qforce*uux1
     fpx(lc)=fpx(lc)+qforce*uux3
     fpx(lb)=fpx(lb)+qforce*uux2
     fpy(la)=fpy(la)+qforce*uuy1
     fpy(lc)=fpy(lc)+qforce*uuy3
     fpy(lb)=fpy(lb)+qforce*uuy2
     fpz(la)=fpz(la)+qforce*uuz1
     fpz(lc)=fpz(lc)+qforce*uuz3
     fpz(lb)=fpz(lb)+qforce*uuz2
     
!!$c--      switch to singularity free potential for 
!!$c--      linear bending V = 2*K(cos(theta) + 1)  
     IF(ABS(pota-pi) < 0.01) THEN 
        ubend(type)=ubend(type)+2.d0*potb*(cb + 1.d0) 
     ELSE
        ubend(type)=ubend(type)+potb*(bb-pota)**2
     END IF
     IF(SIZE(Param_Angles(i) % pot) /= 4) CYCLE
     
!!$*=======================================================================
!!$*---  Compute the Urey-Bradley term ------------------------------------
!!$*=======================================================================
     
     potc=Param_Angles(i) % pot(4)
     potd=Param_Angles(i) % pot(3)*Conv_Fact
     x31=xr3-xr1
     y31=yr3-yr1
     z31=zr3-zr1
     rs31=x31**2+y31**2+z31**2
     rsp31=DSQRT(rs31)
     
     qforce=-2.0D0*potd*(potc-rsp31)
     uux1=x31/rsp31
     uuy1=y31/rsp31
     uuz1=z31/rsp31
     uux2=-uux1
     uuy2=-uuy1
     uuz2=-uuz1
     fpx(la)=fpx(la)+qforce*uux1
     fpy(la)=fpy(la)+qforce*uuy1
     fpz(la)=fpz(la)+qforce*uuz1
     fpx(lc)=fpx(lc)+qforce*uux2
     fpy(lc)=fpy(lc)+qforce*uuy2
     fpz(lc)=fpz(lc)+qforce*uuz2
     ubend(type)=ubend(type)+potd*(rsp31-potc)**2 
  END DO
  Ubend_slt=ubend(1)
  Ubend_slv=ubend(2)
END SUBROUTINE Angles
