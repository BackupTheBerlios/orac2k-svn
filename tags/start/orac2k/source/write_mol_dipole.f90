  subroutine write_mol_dipole(ntap,fstep,kout_dp,kout_dp_sta,kout_dp_tot,unitc&
                             ,protl,xp0,yp0,zp0,charge,dip)
  implicit none

! Argument
  real(8) :: fstep,dip(3,*),unitc,charge(*),xp0(*),yp0(*),zp0(*)
  integer :: ntap,protl(*),kout_dp,kout_dp_sta,kout_dp_tot

! Local
  integer n,m,i,k,ncount
  real(kind=8) :: Cell_sta(3),Cell_ind(3),Atom_sta(3),Atom_ind(3)  
  real(8) :: sum_sta,sum_ind,norma,dip_sta_nor,dip_ind_nor

! Initialization

  sum_sta = 0.0d0
  sum_ind = 0.0d0
  norma = 0.0d0
  Cell_sta = 0.0d0
  Cell_ind = 0.0d0
    
! Calculate

  n=0
  ncount=1

  do while(n .LT. ntap)
     m=protl(ncount)   
     Atom_sta = 0.0d0
     Atom_ind = 0.0d0

     do i=1,m
        do k=1,3
           Atom_ind(k) = Atom_ind(k) + dip(k,i)
        end do
        Atom_sta(1) = Atom_sta(1) + xp0(i)*charge(i)
        Atom_sta(2) = Atom_sta(2) + yp0(i)*charge(i)
        Atom_sta(3) = Atom_sta(3) + zp0(i)*charge(i)
     end do

     dip_sta_nor = 0.0d0
     dip_ind_nor = 0.0d0
     do k=1,3
        dip_sta_nor = dip_sta_nor + Atom_sta(k)**2
        dip_ind_nor = dip_ind_nor + Atom_ind(k)**2
     end do

     dip_sta_nor = dsqrt(dip_sta_nor)*dsqrt(unitc)*4.805d0
     dip_ind_nor = dsqrt(dip_ind_nor)*dsqrt(unitc)*4.805d0

     sum_sta = sum_sta + dip_sta_nor
     sum_ind = sum_ind + dip_ind_nor
     norma = norma + 1.0d0

     Cell_sta = Cell_sta + Atom_sta
     Cell_ind = Cell_ind + Atom_ind

     ncount=ncount+m+1
     n=n+m
  end do !while

  sum_sta = sum_sta / norma
  sum_ind = sum_ind / norma

  write(kout_dp,'(f12.2,3(1x,f9.4))') fstep, sum_sta,sum_sta+sum_ind
  write(kout_dp_sta,'(f12.2,3(1x,f12.3))') fstep, (Cell_sta(k),k=1,3)
  write(kout_dp_tot,'(f12.2,3(1x,f12.3))') fstep, (Cell_sta(k)+Cell_ind(k),k=1,3)

  return
  end 
