      SUBROUTINE change_virtual_potential(ntap,pnbd2,betb,nbtype,r_min
     &     ,Virtual_type,types_virt,virtual_atoms)

      IMPLICIT none

      REAL*8 pnbd2(*),r_min(*)
      CHARACTER*7 betb(*),types_virt(*)
      INTEGER nbtype(*),Virtual_Type,ntap,virtual_atoms(*)

      INTEGER i,j,m

      virtual_atoms(1)=0
      DO i=1,Virtual_type
         DO j=1,ntap
            IF(betb(j) .EQ. types_virt(i)) THEN
               pnbd2(nbtype(j))=r_min(i)
               virtual_atoms(1)=virtual_atoms(1)+1
               m=virtual_atoms(1)
               virtual_atoms(1+m)=j
            END IF
         END DO
      END DO

      RETURN
      END
