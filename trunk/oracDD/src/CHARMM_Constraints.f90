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
SUBROUTINE Constraints(fpx,fpy,fpz)
  REAL(8) :: fpx(:),fpy(:),fpz(:)
  REAL(8), POINTER :: xp0(:),yp0(:),zp0(:)
  INTEGER :: i,la,lb,type,nbond
  REAL(8) ::  xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21,uux1,uux2&
       &,uuy1,uuy2,uuz1,uuz2,ubond(2),pota,potb
  REAL(8) ::  qforce,dsq
  
  IF(Calls == 0) THEN
     Calls=Calls+1
     Conv_Fact=1000.0D0*4.184/(unite*avogad)
  END IF
  
  xp0=>Atoms(:) % x
  yp0=>Atoms(:) % y
  zp0=>Atoms(:) % z
  
  nbond=SIZE(Indx_Bonds,2)
  ubond=0.0_8
  DO i=1,nbond
     la=Indx_Bonds(1,i)
     lb=Indx_Bonds(2,i)
     pota=Param_Bonds(i) % pot(2)
     
     type=Atoms(la) % Id_Slv
     
     xr1=xp0(la)
     yr1=yp0(la)
     zr1=zp0(la)
     xr2=xp0(lb)
     yr2=yp0(lb)
     zr2=zp0(lb)
     
     x21=xr2-xr1
     y21=yr2-yr1
     z21=zr2-zr1
     rs21=DSQRT(x21**2+y21**2+z21**2)
     
     qforce=-2.0D0*(pota-rs21)
     uux1=x21/rs21
     uuy1=y21/rs21
     uuz1=z21/rs21
     uux2=-uux1
     uuy2=-uuy1
     uuz2=-uuz1
     fpx(la)=fpx(la)+qforce*uux1
     fpy(la)=fpy(la)+qforce*uuy1
     fpz(la)=fpz(la)+qforce*uuz1
     
     fpx(lb)=fpx(lb)+qforce*uux2
     fpy(lb)=fpy(lb)+qforce*uuy2
     fpz(lb)=fpz(lb)+qforce*uuz2
     ubond(type)=ubond(type)+(rs21-pota)**2
  END DO
  nbond=SIZE(Indx_Constr,2)
  DO i=1,nbond
     la=Indx_Constr(1,i)
     lb=Indx_Constr(2,i)
     pota=Param_Constr(i) % pot(1)
     
     type=Atoms(la) % Id_Slv
     
     
     
     xr1=xp0(la)
     yr1=yp0(la)
     zr1=zp0(la)
     xr2=xp0(lb)
     yr2=yp0(lb)
     zr2=zp0(lb)
     
     x21=xr2-xr1
     y21=yr2-yr1
     z21=zr2-zr1
     rs21=DSQRT(x21**2+y21**2+z21**2)
     
     qforce=-2.0D0*(pota-rs21)
     uux1=x21/rs21
     uuy1=y21/rs21
     uuz1=z21/rs21
     uux2=-uux1
     uuy2=-uuy1
     uuz2=-uuz1
     fpx(la)=fpx(la)+qforce*uux1
     fpy(la)=fpy(la)+qforce*uuy1
     fpz(la)=fpz(la)+qforce*uuz1
     
     fpx(lb)=fpx(lb)+qforce*uux2
     fpy(lb)=fpy(lb)+qforce*uuy2
     fpz(lb)=fpz(lb)+qforce*uuz2
     ubond(type)=ubond(type)+(rs21-pota)**2
  END DO

  UConstr_slt=ubond(1)
  UConstr_slv=ubond(2)
END SUBROUTINE Constraints
