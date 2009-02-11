/* Energies.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _ENERGIES_H
#define _ENERGIES_H
  CHARACTER(9), SAVE :: ctime='; Time = '&
       &                ,ctot='; Tot  = '&
       &               ,ctpot='; TPot = '&
       &                , clj='; Lj   = '&
       &               ,ccdir='; CDir = '&
       &               ,ccrec='; CRec = '&
       &               ,cstre='; Stre = '&
       &              ,cangle=';Angle = '&
       &               ,cimph='; Imph = '&
       &               ,c14lj='; 14LJ = '&
       &             ,c14coul='; 14Co = '&
       &              ,cdihed=';Dihed = '&
       &               ,cbond=';Bnded = '&
       &                ,ckin='; Kin  = '&
       &               ,ctemp='; Temp = '
  CHARACTER(12), SAVE ::blank12='            '
  CHARACTER(9), SAVE ::blank9=';        '
  CHARACTER(1), SAVE :: pipe='|'
  Type :: Energy_Type
     REAL(8) :: tot=0.0_8,Slv=0.0_8,Slt=0.0_8,Mix=0.0_8
     INTEGER :: n0=0
  END type Energy_Type
#define _fact Dble(PI_Nprocs)
#define __avg(a,b)				\
    a%Tot=a%Tot+b%Tot*_fact;				\
    a%Slv=a%Slv+b%Slv*_fact;				\
    a%Slt=a%Slt+b%Slt*_fact;				\
    a%Mix=a%Mix+b%Mix*_fact;				\
    a%n0=a%n0+1;
#endif
