SUBROUTINE Transform_CHARMM(Residue)
  IMPLICIT NONE 
  TYPE(Resid) :: Residue
  INTEGER :: ikeys,n
  INTEGER, SAVE :: First_time=0
    
  CALL Group_CHARMM(Residue)
  IF(First_time == 0) THEN
     CALL Init_hash_Charmm
     First_Time=1
  END IF
  ikeys=IHash('bond',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%bonds,ikeys)
  ikeys=IHash('imph',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%imph,ikeys)
  ikeys=IHash('acc',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%acc,ikeys)
  ikeys=IHash('acc_',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%acc_,ikeys)
  ikeys=IHash('don',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%don,ikeys)
  ikeys=IHash('don_',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%don_,ikeys)
  ikeys=IHash('dele',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%dele,ikeys)
  ikeys=IHash('term',Hash_Charmm)
  CALL Topol_Charmm(Residue,Resi%ends,ikeys)

  IF(.NOT. ALLOCATED(Resi%ends)) THEN
     CALL Ends_Charmm(Residue,Resi)
  END IF
END SUBROUTINE Transform_CHARMM

SUBROUTINE Topol_CHARMM(Residue,Resi,ikeys)
  IMPLICIT NONE 
  TYPE(Resid) :: Residue
  CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: Resi
  INTEGER :: ikeys
  
  INTEGER :: n,m,o,om
  INTEGER :: count,nword,nstop
  CHARACTER(len=max_char) :: line
  LOGICAL :: ok
  
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     nword=SIZE(strngs)
     ok=.FALSE.
     DO m=1,SIZE(Hash_Charmm(ikeys)%keys)
        line=TRIM(Hash_Charmm(ikeys)%keys(m))
        IF(MY_FAM(line,strngs(1))) ok=.TRUE.
     END DO
     
     IF(ok) THEN
        IF(Hash_Charmm(ikeys)%Method == 0) THEN
           nstop=INT((nword-1)/Hash_Charmm(ikeys)%n)
           count=count+nstop
        ELSE
           count=count+1
        END IF
     END IF
  END DO
  IF(count == 0) RETURN
  ALLOCATE(Resi(Hash_Charmm(ikeys)%n,count))
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     nword=SIZE(strngs)
     ok=.FALSE.
     DO m=1,SIZE(Hash_Charmm(ikeys)%keys)
        line=TRIM(Hash_Charmm(ikeys)%keys(m))
        IF(MY_FAM(line,strngs(1))) THEN
           ok=.TRUE.
        END IF
     END DO
     IF(ok) THEN
        nstop=INT((nword-1)/Hash_Charmm(ikeys)%n)
        IF(Hash_Charmm(ikeys)%Method == 0) THEN
           DO m=1,nstop
              count=count+1
              DO o=1,Hash_Charmm(ikeys)%n
                 om=(m-1)*Hash_Charmm(ikeys)%n+o+1
                 Resi(o,count)=strngs(om)
              END DO
           END DO
        ELSE
           count=count+1
           Resi(1,count)=ST_Concat(2,strngs)
        END IF
     END IF
  END DO
END SUBROUTINE Topol_CHARMM
SUBROUTINE Group_CHARMM(Residue)
  IMPLICIT NONE 
  TYPE(Resid) :: Residue
  INTEGER :: n,m
  INTEGER :: count
  INTEGER, PARAMETER :: max_group=500
  TYPE(list), DIMENSION(:), ALLOCATABLE :: group
  INTEGER, DIMENSION(max_group) :: ig=0

  ig=0
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     IF(MY_FAM('grou',strngs(1))) THEN 
        count=count+1
     END IF
     IF(MY_FAM('atom',strngs(1))) THEN 
        ig(count)=ig(count)+1
     END IF
  END DO
  
  ALLOCATE(Resi%group(count))
  
  DO n=1,count
     ALLOCATE(Resi%group(n)%g(ig(n)))
     ig(n)=0
  END DO
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     IF(MY_FAM('grou',strngs(1))) THEN 
        count=count+1
     END IF
     IF(MY_FAM('atom',strngs(1))) THEN 
        ig(count)=ig(count)+1
        Resi%group(count)%g(ig(count))=ST_Concat(2,strngs)
     END IF
  END DO
END SUBROUTINE Group_CHARMM
SUBROUTINE Ends_Charmm(Residue,Resi)
  TYPE(Resid) :: Residue
  TYPE(Unit_Char) :: Resi
  INTEGER :: count,n,n_extra
  CHARACTER(len=max_char) :: line1,line2
  CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: new_bonds
  
  IF(.NOT. ALLOCATED(Resi%bonds)) RETURN
  ALLOCATE(new_bonds(2,SIZE(Resi%bonds,2)))
  n_extra=-1
  count=0
  DO n=1,SIZE(Resi%bonds,2)
     line1=Resi%bonds(1,n)
     line2=Resi%bonds(2,n)
     IF(MY_FAM('+',line1) .OR. MY_FAM('+',line2)) THEN
        n_extra=n
     ELSE
        count=count+1
        new_bonds(1,count)=line1
        new_bonds(2,count)=line2
     END IF
  END DO
  IF(SIZE(Resi%bonds,2)-count > 1) STOP
  
  IF(n_extra > 0) THEN
     IF(.NOT. ALLOCATED(Resi%ends)) THEN
        n=n_extra
        line1=TRIM(Resi%bonds(1,n))
        line2=TRIM(Resi%bonds(2,n))
        IF(MY_FAM('+',line1)) THEN 
           CALL CHR_Delete('+',line1)
           line1=TRIM(line1)
        END IF
        IF(MY_FAM('+',line2)) THEN 
           CALL CHR_Delete('+',line2)
           line2=TRIM(line2)
        END IF
        ALLOCATE(Resi%ends(2,1))
        Resi%ends(1,1)=line2
        Resi%ends(2,1)=line1
     END IF
     DEALLOCATE(Resi%bonds)
     ALLOCATE(Resi%bonds(2,count))
     Resi%bonds=new_bonds
  END IF

  DEALLOCATE(new_bonds)
END SUBROUTINE Ends_Charmm
SUBROUTINE Transform_MASS(Residue)
  IMPLICIT NONE 
  TYPE(Resid) :: Residue
  INTEGER :: count,n
  CHARACTER(len=*), PARAMETER  :: line='mass'
  
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     IF(MY_FAM(line,strngs(1))) THEN
        count=count+1
     END IF
  END DO
  ALLOCATE(Resi%mass(2,count))
  count=0
  DO n=1,SIZE(Residue%line)
     CALL MY_Parse(Residue%line(n), strngs)
     IF(MY_FAM(line,strngs(1))) THEN
        count=count+1
        Resi%mass(1,count)=ADJUSTL(strngs(3))
        Resi%mass(2,count)=ADJUSTL(strngs(4))
     END IF
  END DO
END SUBROUTINE Transform_MASS

SUBROUTINE Init_Hash_Charmm
  IMPLICIT NONE 
  Hash_Charmm(1)%Type='bond'
  Hash_Charmm(2)%Type='imph'
  Hash_Charmm(3)%Type='acc'
  Hash_Charmm(4)%Type='acc_'
  Hash_Charmm(5)%Type='don'
  Hash_Charmm(6)%Type='don_'
  Hash_Charmm(7)%Type='dele'
  Hash_Charmm(8)%Type='term'
  WRITE(*,*) 'What'
  ALLOCATE(Hash_Charmm(1)%keys(3))
  WRITE(*,*) 'What'
  Hash_Charmm(1)%keys=(/'bond','doub','trip'/)
  Hash_Charmm(1)%n=2
  ALLOCATE(Hash_Charmm(2)%keys(2))
  Hash_Charmm(2)%keys=(/'impr','imph'/)
  Hash_Charmm(2)%n=4
  ALLOCATE(Hash_Charmm(3)%keys(1))
  Hash_Charmm(3)%n=2
  Hash_Charmm(3)%keys=(/'acce'/)
  ALLOCATE(Hash_Charmm(4)%keys(1))
  Hash_Charmm(4)%n=1
  Hash_Charmm(4)%keys=(/'acce'/)
  ALLOCATE(Hash_Charmm(5)%keys(1))
  Hash_Charmm(5)%n=2
  Hash_Charmm(5)%keys=(/'dono'/)
  ALLOCATE(Hash_Charmm(6)%keys(1))
  Hash_Charmm(6)%n=1
  Hash_Charmm(6)%keys=(/'dono'/)
  ALLOCATE(Hash_Charmm(7)%keys(1))
  Hash_Charmm(7)%Method=-1
  Hash_Charmm(7)%n=1
  Hash_Charmm(7)%keys=(/'dele'/)
  ALLOCATE(Hash_Charmm(8)%keys(1))
  Hash_Charmm(8)%n=1
  Hash_Charmm(8)%keys=(/'term'/)
END SUBROUTINE Init_Hash_Charmm
