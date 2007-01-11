MODULE UNITS
  REAL(8), PARAMETER :: pi=3.1415926535897931D0, twopi=2.0D0*pi, avogad=6.0225D23&
       &, boltz=1.38054D-23,gascon=8.3143D0,planck=6.6256D-34,elechg=1.602D-19&
       &,epso=8.854D-12,boxl=2.0D0,unitm=1.0D0/(avogaD*1000.0D0),unitl=1.0D-10&
       &,unitt=1.D-15,unite=unitm*(unitl/unitt)**2&
       &,unitc=4.0D0*pi*epso*unitl*unite/(elechg*elechg)&
       &,unitp=(unite/unitl**3)/1.0D6,efact=unite*avogad&
       &,hartree = 4.35981 * 1.0D-18 &
       &,lbohr = 0.52917706D0
  REAL(8), SAVE :: LJ_Fact,unitepot,unitefield
CONTAINS
  SUBROUTINE Units_
    LJ_fact=1.0D0/(2.0D0**(1.0D0/6.0D0))
    unitepot=unite/unitc**(0.5D0)/hartree
    unitefield=unitepot*lbohr 
  END SUBROUTINE Units_
END MODULE UNITS
