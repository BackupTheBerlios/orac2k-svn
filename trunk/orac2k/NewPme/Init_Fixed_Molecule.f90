SUBROUTINE Init_Fixed_Molecule(mass,nato,vpx,vpy,vpz)
  IMPLICIT none
  INTEGER :: nato
  REAL(8) :: mass(*),vpx(*),vpy(*),vpz(*)
  
  INTEGER :: i

  DO i=1,nato
     IF(mass(i) > 1.0D5) THEN
        vpx(i)=0.0D0
        vpy(i)=0.0D0
        vpz(i)=0.0D0
     END IF
  END DO
END SUBROUTINE Init_Fixed_Molecule
