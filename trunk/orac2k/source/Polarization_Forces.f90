!==============================================================================
!  Program to calculatr self-consistently dipoles on a molecule
!  Matteo Ceccarelli, CECAM-ENSL, November 99
!==============================================================================
SUBROUTINE Polarization_Forces(fscnstr_slt,fscnstr_slv,ntap,chargeb&
     &,ss_index,plrzbij,kpol_inp&
     &,lbohr,lpnbd,type,label,nbtype,protl,anprot,annpro,anpoint&
     &,polar,pol,Ext,xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma&
     &,ypcma,zpcma,co,oc,pmass,mapnl,fpx,fpy,fpz,worka&
     &,cpu,ncpu,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,node&
     &,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs,ncube,fstep&
     &,tag_bndg,ngrp,grppt,ingrpp,ingrp,ingrp_x,errmsg,iret,nstep&
     &,unitc,efact,pi,alphal,U_conf,U_ele,Udirect,Urecip,Uind&
     &,uself_dip,uself,Ugrp,U_Thole,ucoul_slt,ucoul_slv,ucoul_ss&
     &,uconf_slt,uconf_slv,uconf_ss,U_solv,Old_Dipoles,Fixt_Dipoles&
     &,Mesos,Mesos_Rho,volume,Polarization_Models)

  USE Gauss_Overlap_Mod
  implicit none

! Argument
  INTEGER :: ngrp,grppt(2,*),ingrpp,ingrp(2,*),iret,nstep,op
  INTEGER :: annpro,anpoint(2,*),protl(*),ss_index(*)
  CHARACTER(len=80) :: errmsg
  INTEGER :: ntap,kpol_inp,lpnbd,nbtype(*)
  REAL(8) chargeb(*),plrzbij(*),lbohr,Fixt_Dipoles(3,*)&
       &,fscnstr_slt,fscnstr_slv
  LOGICAL :: anprot,polar,print_dp,print_ef,clewld,Old_Dipoles&
       &,Solvation,Mesos
  CHARACTER (len=7) :: label(*),type(*)
  INTEGER :: node,nprocs,ncube,mapnl(*)
  INTEGER :: ingrp_x(*),tag_bndg(*)
  INTEGER :: nodex,nodey,nodez,ictxt,npy,npz,descQ(*)
  REAL(8) :: pol(*),Ext(3),co(3,3),oc(3,3),pmass(*)
  REAL(8) ::  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)
  REAL(8) :: xpga(*),ypga(*),zpga(*),xp0(*),yp0(*),zp0(*)
  REAL(8) :: fpx(*),fpy(*),fpz(*),cpu,fstep,unitc,efact,pi&
       &,alphal
  INTEGER :: nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,worka(*)
  INTEGER :: npp_m,ncpu
  REAL(8) :: U_conf,U_ele,Udirect,Urecip,Uind,Uself_dip,uself,Ugrp&
       &,U_Thole,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
       &,uconf_ss,U_solv,Mesos_Rho,volume
  CHARACTER(80) :: Polarization_Models(2)
  TYPE(Overlap) :: c

  INTERFACE
     SUBROUTINE Do_Electric_fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
           &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
           &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
           &,Ey_rec,Edx,Edy,Edz,Ez_rec,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
           &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
           &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
           &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Utotal&
           &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
           &,uconf_slt,uconf_slv,uconf_ss,charges,dipole,ntap,What_To_do)
      INTEGER(4) :: node,nprocs,ncube,mapnl(*),ingrpp_ef,ingrp_ef(2,*)
      INTEGER(4) :: nodex,nodey,nodez,ictxt,npy,npz,descQ(*),ingrp_ef_x(*)
      INTEGER(4) :: tag_bndg(*)
      REAL(8) ::  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)&
           &,xpga(*),ypga(*),zpga(*),xp0(*),yp0(*),zp0(*)
      REAL(8) ::  Edx(*),Edy(*),Edz(*),dipole(3,*),charges(*)
      REAL(8) ::  Ex_rec(*),Ey_rec(*),Ez_rec(*)
      REAL(8) ::  Ex_Cor(*),Ey_Cor(*),Ez_Cor(*)
      REAL(8) ::  fpx(*),fpy(*),fpz(*),cpu,Ugrp,Utotal,eer,uself&
           &,ene_ferrf,uself_dip,U_Thole,ELJ,flj_x(*),flj_y(*),flj_z(&
           &*),ucoul_slt,ucoul_slv,ucoul_ss& 
           &,uconf_slt,uconf_slv,uconf_ss
      INTEGER(4) :: nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,worka(*)
      INTEGER(4) :: npp_m,ncpu,ntap
      CHARACTER(80) :: What_To_Do
    END SUBROUTINE Do_Electric_fields
 END INTERFACE

! Local Variables

  INTEGER ii,i,j,ij,k,count,nbti,idx,npol,Count_Max
  REAL(8) :: Cell_sta(3),Cell_ind(3),Atom_sta(3),Atom_ind(3),diff_ene
  REAL(8) :: U_cgd,Ued,Utot,diff,ddot,U_charges
  REAL(8) :: au_to_aa,dip_sta_nor,dip_ind_nor,norma,diff_dip
  REAL(8), SAVE :: apara,gamma,maxdiff,Utot_old,maxdiff0,gamma0
  REAL(8) :: sum_sta,sum_ind,auxe,auxd&
       &,ene_ferrf,ELJ,sumx,sumy,sumz
  INTEGER, allocatable, save :: mapnl2(:)
  LOGICAL :: static,self,first_step=.TRUE.,ok
  CHARACTER(7), DIMENSION (:) , ALLOCATABLE, SAVE :: type_aux
  REAL(8), DIMENSION (:) , ALLOCATABLE, SAVE :: pol_au,pol_au0&
       &,pol_atom
  CHARACTER(7) :: dummy_c
  REAL(8) :: dummy_r,chi,PARAS,paras0,pol_max,U_total
  REAL(8) :: xc,yc,zc,xv,yv,zv,xd,yd,zd,sum_cha,Urecip_old
  REAL(8), allocatable, save :: Etotal(:,:),Etotal0(:,:),Etotal_rec0(:,:)
  REAL(8), DIMENSION(:,:), allocatable, save :: dip_new,dip_old,Dipoles,dip_zero
  REAL(8), dimension(:), allocatable, save :: Edx,Edy,Edz,charge
  REAL(8), dimension(:), allocatable, save :: Grad_x,Grad_y,Grad_z
  REAL(8), dimension(:), allocatable, save :: Ex_rec,Ey_rec,Ez_rec
  REAL(8), dimension(:), allocatable, save :: flj_x,flj_y,flj_z
  REAL(8), dimension(:), allocatable, save :: fd0x,fd0y,fd0z
  REAL(8), dimension(:), allocatable, save :: fd1x,fd1y,fd1z
  REAL(8), dimension(:), allocatable, save :: Ex_Cor,Ey_Cor,Ez_Cor
  REAL(8), dimension(:), allocatable, save :: Ed00x,Ed00y,Ed00z
  REAL(8), dimension(:), allocatable, save :: Ed0x,Ed0y,Ed0z
  REAL(8), dimension(:), allocatable, save :: ene,enedip,ene_rec
  REAL(8), dimension(:), allocatable, save :: Ecx,Ecy,Ecz,qc,xpg&
       &,ypg,zpg
  REAL(8), dimension(:,:), allocatable, save :: Energies
  REAL(8) :: deriv(3),energ(3),num_d,bin,xpgg

  REAL(8) rot(3,3),error,qt(4),sumcha,factor
  REAL(8), allocatable, save :: xyz(:,:),xyz0(:,:),xyzfit(:,:)&
       &,wca2(:),work(:)
  REAL(8), SAVE :: gg,Utot_0,fac_dip,Udirect0,Urecip0
  REAL(8) :: fac_dip2,Tot_Dip_x,Tot_Dip_y,Tot_Dip_z,U_Tot_Dip,Polzed&
       &,aux_pol,Emod
  REAL(8), SAVE :: dip_sat,sat_dip
  INTEGER :: No_Polzed
  INTEGER, SAVE :: No_of_Calls=0,Call_GR=0
  INTEGER :: nword
  CHARACTER(80) :: What_To_Do
  LOGICAL :: do_LJ,do_All,do_Gauss,do_CC,do_CD,do_DD
  LOGICAL, SAVE :: Saturation=.TRUE.
  INTEGER :: mapnl0(10000),ingrpp0,ingrp0(2,30000),ingrp0_x(30000)
  CHARACTER(80) :: line,strngs(40)
  CHARACTER(8) :: fmt
  CHARACTER(1) :: sep(2)=(/' ',','/),comm(2)=(/'(',')'/)

  REAL(8), DIMENSION(:),ALLOCATABLE, SAVE :: hist_E0,hist_E,Hist_GR,Hist_Dip
  INTEGER, DIMENSION(:),ALLOCATABLE, SAVE :: No_hist_E0,No_hist_E&
       &,No_Hist_GR,No_Hist_Dip
  REAL(8), SAVE :: Bins=0.05,Rmax=20.0D0,Max_Dip,Shell=2.5D0
  REAL(8), SAVE :: Debey_to_elang=0.208194346224757D0
  INTEGER, SAVE :: Max_Bin,Max_Shell
  INTEGER, SAVE :: Times_of_Calls=0


! initialize
  Times_of_Calls=Times_of_Calls+1
  if( first_step ) then
!!$     Max_Shell=IDINT(0.5D0*Shell/Bins)+1
     Max_Shell=0
     fac_dip=(4.0/3.0)*alphal**3/dsqrt(pi)
     maxdiff = 1.0d-6


     au_to_aa=lbohr**3
     read(kpol_inp,*) apara
     read(kpol_inp,*) gamma0
     read(kpol_inp,*) dip_sat

     read(kpol_inp,*) maxdiff
     npol=0
     DO 
        read(kpol_inp,'(a7,f6.4)',END=1200,ERR=1300) dummy_c,dummy_r
        npol=npol+1
     END DO
1300 CONTINUE
     WRITE(*,*) 'Error in input'
     STOP
1200 CONTINUE
     ALLOCATE(type_aux(npol),pol_au(npol))
     ALLOCATE(pol_atom(ntap))

     REWIND(kpol_inp)
     read(kpol_inp,*) apara
     read(kpol_inp,*) gamma0
     read(kpol_inp,*) dip_sat
     read(kpol_inp,*) maxdiff
     
     IF(dip_sat > 0.0D0) THEN
        Saturation=.TRUE.
        sat_dip=dip_sat
        sat_dip=sat_dip*Debey_to_Elang/DSQRT(unitc)
     ELSE
        Saturation=.FALSE.
     END IF
  
     pol_max=-1.0d10
     DO i=1,npol
        READ(kpol_inp,'(a78)') line(1:78)
        CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
        type_aux(i)=strngs(1)(1:7)
        CALL fndfmt(2,strngs(2),fmt)
        READ(strngs(2),fmt) pol_au(i)
        IF(mesos) THEN
           pol_au(i)=pol_au(i)/Mesos_rho
           WRITE(*,'(''Polarization in Ang. '',a7,2x,f15.9)') &
                & type_aux(i),pol_au(i)
        ELSE
           pol_au(i)=pol_au(i)*au_to_aa
        END IF
        IF(pol_au(i) .GT. pol_max) pol_max=pol_au(i)
     END DO
     
     DO i=1,ntap
        nbti = nbtype(i)
        ok=.FALSE.
        DO j=1,npol
           if(type_aux(j) .EQ. type(nbti)) then
              pol(nbti) = pol_au(j)
              ok=.TRUE.
           end if
        END DO
        IF(.NOT. ok) THEN
           DO j=1,npol
              IF(type_aux(j)(2:2) == '=') THEN
                 if(type_aux(j)(1:1) .EQ. type(nbti)(1:1)) then
                    pol(nbti) = pol_au(j)
                    ok=.TRUE.
                 END if
              END IF
           END DO
        END IF
        IF(.NOT. ok) THEN
           write(6,*)
           write(6,*) 'TYPE MISMATCH'
           write(6,*) 'Type ', type(nbti),'does not match polarization types '
           write(6,*)
           STOP
        END IF
     END DO
     IF(Saturation) THEN
        WRITE(*,'(''Dipoles *with* saturation threshold ='',f12.4)') dip_sat
     ELSE
        WRITE(*,'(''Dipoles *without* saturation threshold'')')
     END IF

     DO i=1,ntap
        nbti = nbtype(i)
        pol_atom(i)=pol(nbti)
     END DO
        

     do i=1,lpnbd
        do j=i,lpnbd
           ij = j*(j-1)/2 + i
           plrzbij(ij) = (pol(i)*pol(j))**(1.d0/6.d0)
           plrzbij(ij) = plrzbij(ij) / apara
        end do
     end do
     write(6,*)
     write(6,*) 'Polarization: Initialization OK'
     write(6,*)
     
     allocate(dip_old(1:3,1:ntap))
     allocate(Dipoles(1:3,1:ntap))
     allocate(Dip_zero(1:3,1:ntap))
     Dip_zero=0.0D0
     allocate(dip_new(1:3,1:ntap))
     allocate(Etotal(1:3,1:ntap))
     allocate(Etotal0(1:3,1:ntap))
     allocate(Etotal_rec0(1:3,1:ntap))
     
     allocate(Ex_rec(1:ntap))
     allocate(Ey_rec(1:ntap))
     allocate(Ez_rec(1:ntap))
     allocate(Ex_Cor(1:ntap))
     allocate(Ey_Cor(1:ntap))
     allocate(Ez_Cor(1:ntap))
     allocate(ene(1:ntap))
     allocate(enedip(1:ntap))
     allocate(ene_rec(1:ntap))
     allocate(Edx(1:ntap),Edy(1:ntap),Edz(1:ntap))
     allocate(Ed0x(1:ntap),Ed0y(1:ntap),Ed0z(1:ntap))
     allocate(Ed00x(1:ntap),Ed00y(1:ntap),Ed00z(1:ntap))
     allocate(xpg(1:ngrp))
     allocate(ypg(1:ngrp))
     allocate(zpg(1:ngrp))
     allocate(qc(1:ntap))
     allocate(Energies(20,300))
     allocate(Grad_x(1:ntap))
     allocate(Grad_y(1:ntap))
     allocate(Grad_z(1:ntap))
     allocate(flj_x(1:ntap))
     allocate(flj_y(1:ntap))
     allocate(flj_z(1:ntap))
     allocate(charge(1:ntap))
     allocate(fd0x(1:ntap))
     allocate(fd0y(1:ntap))
     allocate(fd0z(1:ntap))
     allocate(fd1x(1:ntap),fd1y(1:ntap),fd1z(1:ntap))
     DO i=1,ntap
        charge(i)=chargeb(i)
        IF(charge(i) == 0.0D0) charge(i)=1.0D-15
     END DO
     gg=efact/1000.0D0
     IF(mesos) THEN
        Max_Bin=INT(Rmax/Bins)+1
        ALLOCATE(hist_E0(1:Max_Bin),hist_E(1:Max_Bin)&
             &,hist_GR(1:Max_Bin),hist_Dip(1:Max_Bin))
        ALLOCATE(No_hist_E0(1:Max_Bin),No_hist_E(1:Max_Bin)&
             &,No_hist_GR(1:Max_Bin),No_hist_Dip(1:Max_Bin))
        hist_E=0.0D0
        hist_E0=0.0D0
        hist_GR=0.0D0
        hist_Dip=0.0D0
        No_hist_E=0
        No_hist_E0=0
        No_hist_GR=0
        No_hist_Dip=0
     END IF
  end if
  dip_old = 0.0d0
!!$  IF(Polarization_Models(2)(1:5) == 'Gauss') THEN
!!$     CALL Gauss_init(ntap, ngrp, grppt, ss_index, nnlpp0, xpa, ypa,&
!!$          & zpa, co, oc, Dip_old, pol_atom, alphal)
!!$  END IF
!   Now calculate the external electric field that dipoles interact with
!   mapnl not zero and reciprocal part

  
!!$     CALL Get_Solvent_Out(charge,qc,ss_index,ntap)
!!$     qc=charge

!!$  No_of_Calls=No_of_Calls+1

!!$========================================================================
!!$---- Compute contribution from charges only
!!$========================================================================
  
  What_To_DO='Charges Direct and LJ'
  CALL Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
       &,xpcma,ypcma,zpcma,mapnl,ingrpp,ingrp,ingrp_x,fd0x,fd0y,fd0z&
       &,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Ed00x,Ed00y,Ed00z&
       &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu&
       &,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex&
       &,nodey,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg&
       &,Ugrp,Urecip0,uself,uself_dip,Udirect0,ene_ferrf,U_Thole,ELJ&
       &,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
       &,uconf_ss,charge,dip_old,ntap,What_To_Do)

  U_conf=ELJ

  Utot_0 = Udirect0 + Ugrp 

  do i=1,ntap
     Etotal0(1,i) = Ed00x(i) 
     Etotal0(2,i) = Ed00y(i) 
     Etotal0(3,i) = Ed00z(i) 
     Ed0x(i)=Ed00x(i)
     Ed0y(i)=Ed00y(i)
     Ed0z(i)=Ed00z(i)
  END do
  What_To_DO='Charges Reciprocal no LJ'
  CALL Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
       &,xpcma,ypcma,zpcma,mapnl,ingrpp,ingrp,ingrp_x,fd1x,fd1y,fd1z&
       &,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Ed00x,Ed00y,Ed00z&
       &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu&
       &,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex&
       &,nodey,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg&
       &,Ugrp,Urecip0,uself,uself_dip,Udirect0,ene_ferrf,U_Thole,ELJ&
       &,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
       &,uconf_ss,charge,dip_old,ntap,What_To_Do)

  Urecip0 = Urecip0 + fscnstr_slt+fscnstr_slv
  Utot_0 = Utot_0 + uself + Urecip0

  DO i=1,ntap
     Etotal0(1,i) = Etotal0(1,i) + Ex_rec(i) + Ex_Cor(i)
     Etotal0(2,i) = Etotal0(2,i) + Ey_rec(i) + Ey_Cor(i)
     Etotal0(3,i) = Etotal0(3,i) + Ez_rec(i) + Ez_Cor(i)
     Etotal_rec0(1,i) = Ex_rec(i) + Ex_Cor(i)
     Etotal_rec0(2,i) = Ey_rec(i) + Ey_Cor(i)
     Etotal_rec0(3,i) = Ez_rec(i) + Ez_Cor(i)
  END DO


  IF(Old_Dipoles .OR. (.NOT. first_step)) THEN
     Dip_old=Dipoles
  ELSE
     do i=1,ntap
        ddot = 0.0d0
        nbti = nbtype(i)
        aux_pol=pol(nbti)
        IF(Saturation) THEN
           Emod=DSQRT(Etotal0(1,i)**2+Etotal0(2,i)**2+Etotal0(3,i)**2)
           aux_pol=(sat_dip/Emod)*langevin(3.0D0*pol(nbti)*Emod&
                &/sat_dip)
        END IF
        do k=1,3
           dip_old(k,i) = aux_pol * (Etotal0(k,i)+Ext(k))
        end do
     END DO
     U_conf=ELJ
  END IF

  count = 0
  diff  = 1.0d5
  ! begin to iterate
  Utot =0.0D0
  count=0
  Count_Max=400000
!!$  IF(First_Step) Count_Max=10000
!!$  IF(.NOT. First_Step) maxdiff=1.0D-30
  do while(diff .gt. maxdiff .AND. count .LT. Count_Max)
     count = count + 1
     What_To_Do=Polarization_Models(2)
     CALL Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
          &,xpcma,ypcma,zpcma,mapnl,ingrpp,ingrp,ingrp_x,fpx,fpy,fpz&
          &,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Edx,Edy,Edz&
          &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu,nstart&
          &,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex,nodey&
          &,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg,Ugrp&
          &,Urecip,uself,uself_dip,Udirect,ene_ferrf,U_Thole,ELJ&
          &,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
          &,uconf_ss,charge,dip_old,ntap,What_To_Do)

     IF(Polarization_Models(2)(1:5) == 'Gauss') THEN
        WRITE(*,*) 9.0D0
        DO i=1,1000
           c=Gauss_Overlap(8.0D0)
           WRITE(*,*) c%x(2)
        END DO
        c=Gauss_Overlap(5.0D0)
        WRITE(*,*) c%x(1:10)
        STOP
     END IF
     Ued = 0.0d0
     Uind = 0.0d0
     U_total=0.0D0

     Tot_dip_x = 0.0d0
     Tot_dip_y = 0.0d0
     Tot_dip_z = 0.0d0
!!$     DO i=1,ntap
!!$        Tot_dip_x=Tot_dip_x+dip_old(1,i)
!!$        Tot_dip_y=Tot_dip_y+dip_old(2,i)
!!$        Tot_dip_z=Tot_dip_z+dip_old(3,i)
!!$     END DO
!!$     fac_dip2=-(4.0D0*pi/3.0D0)/volume

     fac_dip2=0.0D0

     DO i=1,ntap
        ddot = 0.0d0
        nbti = nbtype(i)
        Etotal(1,i) = Edx(i) + Ex_rec(i) + Ex_Cor(i) +&
             & fac_dip*dip_old(1,i) + fac_dip2*Tot_dip_x + Ext(1) 
        Etotal(2,i) = Edy(i) + Ey_rec(i) + Ey_Cor(i) +&
             & fac_dip*dip_old(2,i) + fac_dip2*Tot_dip_y + Ext(2)
        Etotal(3,i) = Edz(i) + Ez_rec(i) + Ez_Cor(i) +&
             & fac_dip*dip_old(3,i) + fac_dip2*Tot_dip_z + Ext(3)
        Ued=Ued-Ext(1)*dip_old(1,i)-Ext(2)*dip_old(2,i)-Ext(3)*dip_old(3,i)
        gamma=gamma0

        aux_pol=pol(nbti)
        IF(Saturation) THEN
           Emod=DSQRT(Etotal(1,i)**2+Etotal(2,i)**2+Etotal(3,i)**2)
           aux_pol=(sat_dip/Emod)*langevin(3.0D0*pol(nbti)*Emod&
                &/sat_dip) 
        END IF
        DO k=1,3
           ddot = ddot + dip_old(k,i)*dip_old(k,i)
           dip_new(k,i) = gamma * aux_pol * Etotal(k,i)
           dip_new(k,i) = dip_new(k,i) + (1.0d0-gamma) * dip_old(k,i)
        END DO
        Uind = Uind + 0.50d0*ddot/pol(nbti)
        Grad_x(i)=dip_old(1,i)/pol(nbti)-Etotal(1,i)
        Grad_y(i)=dip_old(2,i)/pol(nbti)-Etotal(2,i)
        Grad_z(i)=dip_old(3,i)/pol(nbti)-Etotal(3,i)
     END DO
     U_Tot_Dip=0.0D0
!!$     U_Tot_Dip=0.5D0*fac_dip2*(Tot_Dip_x**2+Tot_Dip_y**2+Tot_Dip_z**2)

     Urecip=Urecip + fscnstr_slt + fscnstr_slv
     Uind=Uind+U_tot_dip
     Udirect=Udirect+Udirect0
     Utot_old = Udirect + Urecip + Uind + uself_dip + uself + Ugrp +&
          & U_Thole + Ued
     diff_ene = dabs(Utot_old-Utot)*efact/1000.0d0
     diff_dip = diff_dipole(ntap,dip_new,dip_old)*dsqrt(unitc)*4.805d0
     diff = diff_dip
     Utot = Utot_old

     IF(diff .GT. maxdiff) THEN
        dip_old = dip_new
     END IF
     if(node.eq.0) then
        write(*,'(i4,2x,f15.4,f19.5,f16.7,2x,e16.7,f16.7,e16.7)')  count,diff_ene&
             &,Utot*efact/1000.0D0,diff_dip,diff_dip,U_Thole*efact/1000.0D0&
             &,Grad_x(2)
     end if
  end do ! while

  U_ele=Utot

  
  Dipoles=Dip_old
!!$  WRITE(41,'(''Fields  '',i8,3e18.8)') nstep,Ext(1),Ext(2),Ext(3)
!!$  DO i=1,ntap
!!$     Etotal(:,i)=Etotal(:,i)-Ext(:)
!!$     WRITE(41,'(i8,3(e18.8,2x))') i,Etotal(1,i)-Ext(1),Etotal(2,i)&
!!$          &-Ext(2),Etotal(3,i)-Ext(3) 
!!$     WRITE(42,'(i8,3(e18.8,2x))') i,Etotal(1,i),Etotal(2,i),Etotal(3,i) 
!!$  END DO
!!$  STOP
  IF(Polarization_Models(2) .EQ. 'Gauss no LJ') THEN
     Etotal0(1,:) = Etotal0(1,:) - Ed0x(:)
     Etotal0(2,:) = Etotal0(2,:) - Ed0y(:)
     Etotal0(3,:) = Etotal0(3,:) - Ed0z(:)
  ELSE IF(Polarization_Models(2) .EQ. 'Full no LJ') THEN
     fd0x=0.0D0
     fd0y=0.0D0
     fd0z=0.0D0
  END IF
  Utot=-0.5D0*SUM(Etotal0(1,:)*dipoles(1,:)+Etotal0(2,:)*dipoles(2,:)&
       &+Etotal0(3,:)*dipoles(3,:)) + Utot_0
  WRITE(*,*) 'Rule ',Utot*efact/1000.0D0
  fpx(1:ntap)=fpx(1:ntap)+fd0x(1:ntap)+flj_x(1:ntap)
  fpy(1:ntap)=fpy(1:ntap)+fd0y(1:ntap)+flj_y(1:ntap)
  fpz(1:ntap)=fpz(1:ntap)+fd0z(1:ntap)+flj_z(1:ntap)

  Tot_Dip_x=Tot_Dip_x/DBLE(ntap-1)
  Tot_Dip_y=Tot_Dip_y/DBLE(ntap-1)
  Tot_Dip_z=Tot_Dip_z/DBLE(ntap-1)
  Polzed=DSQRT(Tot_Dip_x**2+Tot_Dip_y**2+Tot_Dip_z**2)

  CALL Get_Max_Dip(Dipoles,Max_Dip)
  WRITE(*,'(''Maximum Dipole = '',f14.5,e16.9)') Max_Dip*DSQRT(unitc)&
       &/Debey_to_elang,Max_dip
  WRITE(*,*) 'Sails ',U_conf*efact/1000.0D0,(U_ele+U_conf)*efact/1000.0D0
  first_step=.FALSE.
  IF(MOD(No_of_Calls,10) == 0) THEN
     U_solv=Utot-Utot_0
  ELSE
     U_solv=0.0D0
  END IF
  CALL Do_Field_Histo0(Etotal0,hist_E0,No_hist_E0)
  CALL Do_Field_Histo0(Etotal,hist_E,No_hist_E)
  CALL Do_Field_Histo0(Dipoles,hist_Dip,No_hist_Dip)
  CALL Do_GR_Histo(hist_GR,No_hist_GR)

  IF(MOD(Times_of_Calls,1) .EQ. 0) THEN
     OPEN(31,file='Hist_E0.dat')
     OPEN(32,file='Hist_E.dat')
     OPEN(33,file='Hist_GR.dat')
     OPEN(34,file='Hist_Dip.dat')
     CALL Write_Field_Histo(31,hist_E0,No_hist_E0)
     CALL Write_Field_Histo(32,hist_E,No_hist_E)
     CALL Write_GR_Histo(33,hist_GR,No_hist_GR)
     CALL Write_Field_Histo(34,hist_Dip,No_hist_Dip)
     CLOSE(31)
     CLOSE(32)
     CLOSE(33)
     CLOSE(34)
  END IF

CONTAINS
!!$=======================================================================
  REAL(8) function langevin(x)
    implicit none
    REAL(8) :: x
    langevin=(1.0D0/tanh(x))-1.0D0/x
  end function langevin


  REAL(8) function diff_dipole(ntap,dip_new,dip_old)
    implicit none
    
    REAL(8) :: dip_new(3,*),dip_old(3,*)
    INTEGER :: ntap
    
    REAL(8) :: norma,diff,sum2,sum
    INTEGER :: i,k
    
    sum = 0.0d0
    sum2 = 0.0d0
    do i=1,ntap
       norma = 0.0d0
       diff=0.0D0
       do k=1,3
          diff=diff + (dip_new(k,i)-dip_old(k,i))**2
          norma = norma + dip_old(k,i)**2
       end do
       sum = sum + DSQRT(norma)
       sum2 = sum2 + diff
    end do
    sum=sum/DBLE(ntap)
    sum2=DSQRT(sum2/DBLE(ntap))
    diff_dipole = 100.0D0*sum2/sum
    diff_dipole = sum2
    return
  end function diff_dipole
!=======================================================================
  SUBROUTINE Get_Solvent_Out(charge,chargeb,ss_index,nato)
    IMPLICIT NONE 
    REAL(8) :: charge(*),chargeb(*)
    INTEGER :: ss_index(*),nato

    INTEGER :: i
    
    DO i=1,nato
       chargeb(i)=charge(i)
       IF(ss_index(i) == 2) THEN
          chargeb(i)=0.0D0
       END IF
    END DO
  END SUBROUTINE Get_Solvent_Out

  SUBROUTINE Do_Field_Histo0(E,hist,no_hist)
    IMPLICIT NONE 
    REAL(8), DIMENSION(:,:) :: E
    REAL(8), DIMENSION(:) :: hist
    INTEGER, DIMENSION(:) :: no_hist
    REAL(8) :: xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,xg,yg,zg
    
    INTEGER :: i_dex,j
    REAL(8) :: EE,Ex,Ey,Ez,ver(3)
    
    xpi=xpa(1)
    ypi=ypa(1)
    zpi=zpa(1)

    DO j=2,ntap
       xg=xpi-xpa(j)
       yg=ypi-ypa(j)
       zg=zpi-zpa(j)
       xg=xg-2.0D0*PBC(xg)
       yg=yg-2.0D0*PBC(yg)
       zg=zg-2.0D0*PBC(zg)
       xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
       yc=           co(2,2)*yg+co(2,3)*zg
       zc=                      co(3,3)*zg
       rsq=xc**2+yc**2+zc**2
       rsp=DSQRT(rsq)
       ver(1)=xc/rsp
       ver(2)=yc/rsp
       ver(3)=zc/rsp
       IF(rsp .LE. Rmax) THEN
          i_dex=INT(rsp/Bins)+1
          Ex=E(1,j)
          Ey=E(2,j)
          Ez=E(3,j)
          EE=DSQRT(Ex**2+Ey**2+Ez**2)
          hist(i_dex)=hist(i_dex)+EE
          no_hist(i_dex)=no_hist(i_dex)+1
       END IF
    END DO

  END SUBROUTINE Do_Field_Histo0
  SUBROUTINE Do_Field_Histo(E,hist,no_hist)
    IMPLICIT NONE 
    REAL(8), DIMENSION(:,:) :: E
    REAL(8), DIMENSION(:) :: hist
    INTEGER, DIMENSION(:) :: no_hist
    REAL(8) :: xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,xg,yg,zg
    
    INTEGER :: i_dex,j,ii,k,Start,End
    REAL(8) :: EE,Ex,Ey,Ez,ver(3)
    
    xpi=xpa(1)
    ypi=ypa(1)
    zpi=zpa(1)

    DO j=2,ntap
       xg=xpi-xpa(j)
       yg=ypi-ypa(j)
       zg=zpi-zpa(j)
       xg=xg-2.0D0*PBC(xg)
       yg=yg-2.0D0*PBC(yg)
       zg=zg-2.0D0*PBC(zg)
       xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
       yc=           co(2,2)*yg+co(2,3)*zg
       zc=                      co(3,3)*zg
       rsq=xc**2+yc**2+zc**2
       rsp=DSQRT(rsq)
       ver(1)=xc/rsp
       ver(2)=yc/rsp
       ver(3)=zc/rsp
       IF(rsp .LE. Rmax) THEN
          i_dex=INT(rsp/Bins)+1
          Ex=E(1,j)
          Ey=E(2,j)
          Ez=E(3,j)
          EE=DSQRT(Ex**2+Ey**2+Ez**2)
          Start=-Max_Shell+1
          End=Max_shell+1
          IF(i_dex+Start .LE. 0) Start=-i_dex+1
          IF(i_dex+End .GT. Max_Bin) End=Max_Bin-i_dex
          DO k=Start,End
             ii=i_dex+k
             hist(ii)=hist(ii)+EE
             no_hist(ii)=no_hist(ii)+1
          END DO
       END IF
    END DO

  END SUBROUTINE Do_Field_Histo
  SUBROUTINE Do_GR_Histo(hist,no_hist)
    IMPLICIT NONE 
    REAL(8), DIMENSION(:) :: hist
    INTEGER, DIMENSION(:) :: no_hist
    REAL(8) :: xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,xg,yg,zg
    
    INTEGER :: i_dex,j
    REAL(8) :: EE,xmap,ymap,zmap
    
    xpi=xpa(1)
    ypi=ypa(1)
    zpi=zpa(1)
    DO j=2,ntap
       xg=xpi-xpa(j)
       yg=ypi-ypa(j)
       zg=zpi-zpa(j)
       xmap=-2.0D0*PBC(xg)
       ymap=-2.0D0*PBC(yg)
       zmap=-2.0D0*PBC(zg)
       xg=xg+xmap
       yg=yg+ymap
       zg=zg+zmap
       xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
       yc=           co(2,2)*yg+co(2,3)*zg
       zc=                      co(3,3)*zg
       rsq=xc**2+yc**2+zc**2
       rsp=DSQRT(rsq)
       IF(rsp .LE. Rmax) THEN
          i_dex=INT(rsp/Bins)+1
          EE=1.0D0
          hist(i_dex)=hist(i_dex)+EE
       END IF
    END DO

  END SUBROUTINE Do_GR_Histo
  SUBROUTINE Write_Field_Histo(Unit,hist,no_hist)
    IMPLICIT NONE 
    REAL(8), DIMENSION(:) :: hist
    INTEGER, DIMENSION(:) :: no_hist
    INTEGER :: Unit

    REAL(8) :: EE,rx
    INTEGER, SAVE :: No_Calls=0
    
    No_Calls=No_Calls+1
    
    rx=0.0D0
    DO i=1,Max_Bin
       rx=Bins*(i-1)
       IF(no_hist(i) .NE. 0) THEN
          EE=hist(i)/DBLE(no_hist(i))
          WRITE(Unit,'(f12.4,e18.8)') rx,EE
       END IF
    END DO
  END SUBROUTINE Write_Field_Histo
  SUBROUTINE Write_GR_Histo(Unit,hist,no_hist)
    IMPLICIT NONE 
    REAL(8), DIMENSION(:) :: hist
    INTEGER, DIMENSION(:) :: no_hist
    INTEGER :: Unit

    REAL(8) :: EE,rx,dv,facrdf
    INTEGER, SAVE :: No_Calls=0
    
    facrdf=1.0D0/(4.0d0*pi*Bins)
    
    DO i=1,Max_Bin
       rx=Bins*(i-1)
       dv=4.0D0*pi*rx**2
       IF(hist(i) .NE. 0.0D0) THEN
          EE=facrdf*hist(i)/DBLE(Times_of_Calls)/(rx*rx)
          WRITE(Unit,'(f12.4,e18.8)') rx,EE
       END IF
    END DO
  END SUBROUTINE Write_GR_Histo
  SUBROUTINE Get_Max_Dip(Dipoles,Max_Dip)
    IMPLICIT NONE 
    REAL(8), DIMENSION (:,:) :: Dipoles
    REAL(8) :: Max_Dip

    INTEGER :: i 
    REAL(8) :: Dipx,Dipy,Dipz,Dips
    Max_dip=-1.0D90
    DO i=1,ntap
       Dipx=Dipoles(1,i)
       Dipy=Dipoles(2,i)
       Dipz=Dipoles(3,i)
       Dips=DSQRT(Dipx**2+Dipy**2+Dipz**2)
       IF(Max_Dip .LT. Dips) Max_Dip=Dips
    END DO
  END SUBROUTINE Get_Max_Dip

  REAL(8) FUNCTION PBC(x)
    REAL(8) :: x
    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
  SUBROUTINE Do_Finite_Differences_Field

    xpgg=dip_old(2,3)
    bin=0.000001D0
    count=0
    DO i=-1,1
       count=count+1
       dip_old(2,3)=xpgg+bin*i
       CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
       What_To_Do=Polarization_Models(2)
       CALL Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
            &,xpcma,ypcma,zpcma,mapnl,ingrpp,ingrp,ingrp_x,fpx,fpy,fpz&
            &,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Edx,Edy,Edz&
            &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu,nstart&
            &,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex,nodey&
            &,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg,Ugrp&
            &,Urecip,uself,uself_dip,Udirect,ene_ferrf,U_Thole,ELJ&
            &,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
            &,uconf_ss,charge,dip_old,ntap,What_To_Do)
       
       Ued = 0.0d0
       Uind = 0.0d0
       U_total=0.0D0
       do j=1,ntap
          ddot = 0.0d0
          nbti = nbtype(j)
          Ex_rec(j) = Ex_rec(j) + (4.0/3.0)*alphal**3/dsqrt(pi)*dip_old(1,j)
          Ey_rec(j) = Ey_rec(j) + (4.0/3.0)*alphal**3/dsqrt(pi)*dip_old(2,j)
          Ez_rec(j) = Ez_rec(j) + (4.0/3.0)*alphal**3/dsqrt(pi)*dip_old(3,j)
          Etotal(1,j) = Edx(j) + Ex_rec(j) + Ex_Cor(j)
          Etotal(2,j) = Edy(j) + Ey_rec(j) + Ey_Cor(j)
          Etotal(3,j) = Edz(j) + Ez_rec(j) + Ez_Cor(j)
          do k=1,3
             ddot = ddot + dip_old(k,j)*dip_old(k,j)
          end do
          Uind = Uind + 0.50d0*ddot/pol(nbti)
       end do
       
       Uind = 0.0d0
       
       Utot = Udirect + Urecip + Uind + uself_dip + uself + Ugrp + U_Thole
       energ(count)=Utot
       deriv(count)=Etotal(2,3)
    END DO
    WRITE(*,*) 'Analytical ',deriv(2)
    WRITE(*,*) 'Numerical  ',-0.5D0*(energ(3)-energ(1))/bin
    
    STOP
    
  END SUBROUTINE Do_Finite_Differences_Field
  SUBROUTINE Do_Finite_Differences_Force
    IMPLICIT none
    REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: dip

    ALLOCATE(dip(3,ntap))

    dip=dip_old
    xpgg=xp0(6)
    bin=0.0001D0
    count=0
    
    What_To_Do=Polarization_Models(2)
    WRITE(*,*) What_To_Do
    DO i=-1,1
       count=count+1
       xp0(6)=xpgg+bin*i
       CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
       
       CALL Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
            &,xpcma,ypcma,zpcma,mapnl,ingrpp,ingrp,ingrp_x,fpx,fpy,fpz&
            &,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Edx,Edy,Edz&
            &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu,nstart&
            &,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex,nodey&
            &,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg,Ugrp&
            &,Urecip,uself,uself_dip,Udirect,ene_ferrf,U_Thole,ELJ&
            &,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
            &,uconf_ss,charge,dip,ntap,What_To_Do)
       U_conf=0.0D0
       Ued = 0.0d0
       Uind = 0.0d0
       U_total=0.0D0
       do j=1,ntap
          ddot = 0.0d0
          nbti = nbtype(j)
          do k=1,3
             ddot = ddot + dip(k,j)*dip(k,j)
          end do
          Uind = Uind + 0.50d0*ddot/pol(nbti)
       end do
       
       Utot = Udirect + Ugrp + Urecip + Uind + uself_dip + uself + U_Thole + U_conf
       energ(count)=Utot
       deriv(count)=fpx(6)
       WRITE(*,*) Utot
    END DO
    WRITE(*,*) 'Analytical ',deriv(2)
    WRITE(*,*) 'Numerical  ',-0.5D0*(energ(3)-energ(1))/bin
    
    DEALLOCATE(dip)
    STOP
  END SUBROUTINE Do_Finite_Differences_Force
END SUBROUTINE Polarization_Forces
