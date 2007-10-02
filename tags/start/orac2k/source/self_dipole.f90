!==============================================================================
!  Program to calculate self-consistently dipoles on a molecule
!  Matteo Ceccarelli, CECAM-ENSL, November 99
!==============================================================================
   SUBROUTINE self_dipole(ntap,charge,plrzbij,kpol_inp,lbohr,lpnbd,type,label &
                         ,nbtype,protl,polar,pol,Ext                          &
                         ,xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma  &
                         ,zpcma,co,mapnl,fpx,fpy,fpz,nnlpp0,npp_m,worka,cpu   &
                         ,ncpu,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a    &
                         ,node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs   &
                         ,ncube,fstep,xpt,ypt,zpt,wca,tag_bndg,ngrp,grppt     &
                         ,ingrpp,ingrp,ingrp_x,mapnl0,ingrpp0,ingrp0,ingrp0_x &
                         ,errmsg,iret,nstep,unitc,print_ef,kout_dp,print_dp   &
                         ,nfreq_dp,kout_dp_sta,kout_dp_tot,polar_counter,efact&
                         ,clewld)

   implicit none

! Argument
  integer :: ngrp,grppt(2,*),ingrpp,ingrp(2,*),iret,nstep
  integer :: protl(*),kout_dp,kout_dp_sta,kout_dp_tot,nfreq_dp
  character(len=80) :: errmsg
  integer :: ntap,kpol_inp,lpnbd,nbtype(*)
  real(kind=8) charge(*),plrzbij(*),lbohr
  logical :: polar,print_dp,print_ef,clewld
  character (len=7) :: label(*),type(*)
  integer :: node,nprocs,ncube,mapnl(*),mapnl0(*),ingrpp0,ingrp0(2,*)
  integer :: ingrp_x,ingrp0_x,tag_bndg(*)
  integer :: nodex,nodey,nodez,ictxt,npy,npz,descQ(*)
  integer(kind=4) :: nnlpp0(*)
  real(kind=8) :: pol(*),Ext(3),co(3,3)
  real(kind=8) ::  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)
  real(kind=8) :: xpga(*),ypga(*),zpga(*),xp0(*),yp0(*),zp0(*)
  real(kind=8) :: xpt(*),ypt(*),zpt(*),wca(*)
  real(kind=8) :: fpx(*),fpy(*),fpz(*),cpu,fstep,unitc,polar_counter,efact
  integer :: nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,worka(*)
  integer :: npp_m,ncpu

! Local Variables

  integer ii,i,j,ij,k,count,nbti,idx
  real(kind=8) :: Cell_sta(3),Cell_ind(3),Atom_sta(3),Atom_ind(3),diff_ene
  real(kind=8) :: ucoul,Ustatic,Udd,Ued,Uind,Utot,diff,diff_dipole,ddot,pol_au
  real(kind=8) :: au_to_aa,dip_sta_nor,dip_ind_nor,norma,Urecip,uself,eer,diff_dip
  real(kind=8) :: apara,gamma,maxdiff,npol,sum_sta,sum_ind,auxe,auxd,Utot_old,norm
  integer, allocatable, save :: mapnl2(:)
  logical :: static,self,first_step=.TRUE.
  character (len=7) :: type_aux
  real(kind=8) :: xc,yc,zc,xv,yv,zv,xd,yd,zd
  real(kind=8), allocatable, save :: Etotal(:,:),Estatic(:,:)
  real(kind=8), allocatable, save :: dip_new(:,:),dip_old(:,:)
  real(kind=8), dimension(:), allocatable, save :: Edx,Edy,Edz,dipx,dipy,dipz
  real(kind=8), dimension(:), allocatable, save :: fx_rec,fy_rec,fz_rec
  real(kind=8), dimension(:), allocatable, save :: ene,enedip,ene_rec

! initialize
  
  if( first_step ) then
      
     maxdiff = 1.0d-6
  
     if(polar) then

        polar_counter = 0.0d0
        au_to_aa=lbohr**3
        read(kpol_inp,*) apara
        read(kpol_inp,*) gamma
        read(kpol_inp,*) maxdiff
        read(kpol_inp,*) npol
        INPUT: do i=1,npol
           read(kpol_inp,'(a7,f6.4)') type_aux,pol_au
           do j=1,lpnbd
              if(type_aux.eq.type(j)) then
                 pol(j) = pol_au*au_to_aa
                 CYCLE INPUT
              end if
           end do
           write(6,*)
           write(6,*) 'WARNING: file kpol_inp'
           write(6,*) 'line:', i, '  error, type does not match'
           write(6,*)
           STOP
        end do INPUT

        do i=1,lpnbd
           do j=i,lpnbd
              ij = j*(j-1)/2 + i
                  plrzbij(ij) = (pol(i)*pol(j))**(1.d0/6.d0)
                  plrzbij(ij) = plrzbij(ij) / apara
           end do
        end do
        write(6,*)
        write(6,*) 'Self_dipole: Initialization OK'
        write(6,*)

        allocate(dip_old(1:3,1:ntap))
        allocate(dip_new(1:3,1:ntap))
        allocate(Etotal(1:3,1:ntap))
        allocate(Estatic(1:3,1:ntap))

     end if ! polar

     allocate(ene(1:ntap))
     allocate(enedip(1:ntap))
     allocate(ene_rec(1:ntap))
     allocate(fx_rec(1:ntap))
     allocate(fy_rec(1:ntap))
     allocate(fz_rec(1:ntap))
     allocate(Edx(1:ntap))
     allocate(Edy(1:ntap))
     allocate(Edz(1:ntap))
     allocate(dipx(1:ntap))
     allocate(dipy(1:ntap))
     allocate(dipz(1:ntap))
     dipx=0.0d0
     dipy=0.0d0
     dipz=0.0d0

  end if

  if(polar) then
!   Now calculate the external electric field that dipoles interact with
!   mapnl not zero and reciprocal part

  self =   .TRUE.
  static = .TRUE.
  polar =  .FALSE.

  CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma  &
                         ,zpcma,mapnl,ingrpp,ingrp,ingrp_x                    &
                         ,fpx,fpy,fpz,fx_rec,fy_rec,fz_rec,nnlpp0,npp_m,worka &
                         ,cpu,ncpu,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a&
                         ,node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs   &
                         ,ncube,tag_bndg,ene,enedip,ene_rec&
                         ,ucoul,eer,uself,Udd,Edx,Edy,Edz,dipx,dipy,dipz,self,static)

     Estatic(1,1:ntap) = (fpx(1:ntap)+fx_rec(1:ntap)) / charge(1:ntap)
     Estatic(2,1:ntap) = (fpy(1:ntap)+fy_rec(1:ntap)) / charge(1:ntap)
     Estatic(3,1:ntap) = (fpz(1:ntap)+fz_rec(1:ntap)) / charge(1:ntap)
     Urecip  = eer+uself
     Ustatic = ucoul+Urecip

     if(first_step) then
        do i=1,ntap
           nbti = nbtype(i)
           dipx(i) = pol(nbti) * (Estatic(1,i)+Ext(1))
           dipy(i) = pol(nbti) * (Estatic(2,i)+Ext(2))
           dipz(i) = pol(nbti) * (Estatic(3,i)+Ext(3))
        end do
     end if

  count = 0
  diff  = 1.0d5

! begin to iterate

  static = .FALSE.
  polar  = .TRUE.
  self   = .FALSE.

  do while(diff .gt. maxdiff)

     count = count + 1
     Ued = 0.0d0
     Uind = 0.0d0
                    
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma  &
                            ,zpcma,mapnl0,ingrpp0,ingrp0,ingrp0_x                &
                            ,fpx,fpy,fpz,fx_rec,fy_rec,fz_rec,nnlpp0,npp_m,worka &
                            ,cpu,ncpu,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a&
                            ,node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs   &
                            ,ncube,tag_bndg,ene,enedip,ene_rec&
                            ,ucoul,eer,uself,Udd,Edx,Edy,Edz,dipx,dipy,dipz,self,static) 

        norm = 0.0d0
        do i=1,ntap
           ddot = 0.0d0
           nbti = nbtype(i)
           Etotal(1,i) = Edx(i) + Estatic(1,i) + Ext(1)
           Etotal(2,i) = Edy(i) + Estatic(2,i) + Ext(2)
           Etotal(3,i) = Edz(i) + Estatic(3,i) + Ext(3)
           dip_old(1,i) = dipx(i)
           dip_old(2,i) = dipy(i)
           dip_old(3,i) = dipz(i)
           do k=1,3
              Ued = Ued - (Estatic(k,i)+Ext(k))*dip_old(k,i)
              ddot = ddot + dip_old(k,i)*dip_old(k,i)
              dip_new(k,i) = gamma * pol(nbti) * Etotal(k,i)
              dip_new(k,i) = dip_new(k,i) + (1.0d0-gamma) * dip_old(k,i)
           end do
           Uind = Uind + 0.50d0*ddot/pol(nbti)
            norm = norm + dsqrt(ddot)
        end do
        norm = norm / ntap

        Utot_old = Ustatic + Ued + Udd + Uind

        diff_ene = dabs(Utot_old-Utot)*efact/1000.0d0
        diff_dip = diff_dipole(ntap,dip_new,dip_old)*dsqrt(unitc)*4.805d0
        diff = diff_dip

        if(node.eq.0) then
        write(100,'(i4,2x,f15.4,2(f15.9))')  count,diff_ene,diff_dip,norm
        end if
                   
        Utot = Utot_old
        dipx(1:ntap) = dip_new(1,1:ntap)
        dipy(1:ntap) = dip_new(2,1:ntap)
        dipz(1:ntap) = dip_new(3,1:ntap)

     end do ! while

     if(print_dp .AND. (node .eq. 0)) then
        if(first_step) then
           write(kout_dp,*) '   Dipole Moments (Debye) '
           write(kout_dp,*) '   Time       Static    Total'
        end if
        IF(MOD(nstep,nfreq_dp) .EQ. 0) THEN
           CALL write_mol_dipole(ntap,fstep,kout_dp,kout_dp_sta,kout_dp_tot,unitc &
                                ,protl,xp0,yp0,zp0,charge,dip_new)
        END IF
     end if

  end if ! polar

! recalculate dipole-dipole interactions: only forpp fnbgrp
! if only efield requested then calculate also clewld (true)
!
     static = .TRUE.
     if(.NOT. polar) then
        polar= .FALSE.
     else
        polar=.TRUE.
     end if
     self   = .FALSE.
     if(polar) clewld = .FALSE.

     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma  &
                            ,zpcma,mapnl0,ingrpp0,ingrp0,ingrp0_x                &
                            ,fpx,fpy,fpz,fx_rec,fy_rec,fz_rec,nnlpp0,npp_m,worka &
                            ,cpu,ncpu,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a&
                            ,node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs   &
                            ,ncube,tag_bndg,ene,enedip,ene_rec&
                            ,ucoul,eer,uself,Udd,Edx,Edy,Edz,dipx,dipy,dipz,self,static)
   
     if(polar) then
        Ued  = 0.0d0
        Uind = 0.0d0
        do i=1,ntap
           ddot = 0.0d0
           nbti = nbtype(i)
           do k=1,3
              Ued = Ued - (Estatic(k,i)+Ext(k))*dip_new(k,i)
              ddot = ddot + dip_new(k,i)*dip_new(k,i)
           end do
           Uind = Uind + 0.50d0*ddot/pol(nbti)
        end do

        Utot = Ustatic + Udd + Ued + Uind
        diff_ene = dabs(Utot_old-Utot)*efact/1000.0d0
   
        polar_counter = polar_counter + count

        write(6,*)
        write(6,200) count,diff_ene,diff_dip
200  FORMAT(' OK Iteration:',i4,2x,'d_ene =',f6.2,2x,'d_dip =',e12.4)
        write(6,*)
        write(6,*) 'Utot   =', Utot*efact/1000.0d0
        write(6,*) 'Ucc    =', Ustatic*efact/1000.0d0
        write(6,*) 'Ued    =', Ued*efact/1000.0d0
        write(6,*) 'Udd    =', Udd*efact/1000.0d0
        write(6,*) 'Uind   =', Uind*efact/1000.0d0
        write(6,*) 'test0 %=', (0.50d0*Ued + Udd + Uind)/Utot
        write(6,*) 

     end if !polar

! Add the reciprocal term and write results
!
     if(clewld) Urecip  = eer+uself
     ucoul = ucoul + Urecip
     ene(1:ntap) = ene(1:ntap) + ene_rec(1:ntap)
     fpx(1:ntap) = fpx(1:ntap) + fx_rec(1:ntap)
     fpy(1:ntap) = fpy(1:ntap) + fy_rec(1:ntap)
     fpz(1:ntap) = fpz(1:ntap) + fz_rec(1:ntap)

     if(polar) clewld = .TRUE.
 
  if(print_ef) then

     CALL write_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpt,ypt   &
                ,zpt,wca,fstep,fpx,fpy,fpz,Edx,Edy,Edz,dipx      &
                ,dipy,dipz,ene,enedip,ucoul,Udd,Ued,Uind,node)
  end if

  first_step=.FALSE.
  
  RETURN
  END

!=======================================================================
  real(kind=8) function diff_dipole(ntap,dip_new,dip_old)
  implicit none

  real(kind=8) :: dip_new(3,*),dip_old(3,*)
  integer :: ntap

  real(kind=8) :: norma_new,norma_old,sum
  integer :: i,k

  do i=1,ntap
     norma_new = 0.0d0
     norma_old = 0.0d0
     do k=1,3
        norma_new = norma_new + dip_new(k,i)**2
        norma_old = norma_old + dip_old(k,i)**2
     end do
     norma_new = dsqrt(norma_new)
     norma_old = dsqrt(norma_old)
     sum = sum + (norma_new-norma_old)**2
  end do
  sum = sum / ntap

  diff_dipole = dsqrt(sum/ntap)

  return
  end
!=======================================================================
