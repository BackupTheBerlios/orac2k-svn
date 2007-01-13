!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
    SUBROUTINE Mass_
      INTEGER :: n,nword,c

      c=0
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=MYParse_(Ri % line(n))
         IF(TRIM(strngs(1)) /= 'mass') CYCLE
         IF(nword < 4) CYCLE
         c=c+1
      END DO
      IF(c==0) THEN
         errmsg_f='No Atomic Masses found while reading Charmm Topology'
         CALL Add_Errors(-1, errmsg_f)
         RETURN
      END IF
      ALLOCATE(Res_Char(i_L) % mass (2,c))
      c=0
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=MYParse_(Ri % line(n))
         IF(TRIM(strngs(1)) /= 'mass') CYCLE
         IF(nword < 4) CYCLE
         c=c+1
         Res_Char(i_L) % mass(1,c)=ADJUSTL(strngs(3))
         Res_Char(i_L) % mass(2,c)=ADJUSTL(strngs(4))
      END DO
    END SUBROUTINE Mass_
    SUBROUTINE Charmm_

      ikeys=Hash_Tops__Extract('bond')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % bonds(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % bonds=share
      END IF
      ikeys=Hash_Tops__Extract('imph')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % imph(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % imph=share
      END IF
      ikeys=Hash_Tops__Extract('acc')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % acc(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % acc=share
      END IF
      ikeys=Hash_Tops__Extract('acc_')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % acc_(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % acc_=share
      END IF
      ikeys=Hash_Tops__Extract('don')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % don(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % don=share
      END IF
      ikeys=Hash_Tops__Extract('don_')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % don_(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % don_=share
      END IF
      ikeys=Hash_Tops__Extract('dele')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % dele(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % dele=share
      END IF
      ikeys=Hash_Tops__Extract('term')
      CALL Topol_Charmm_
      IF(ALLOCATED(share)) THEN
         ALLOCATE(Res_Char(i_L) % ends(Ikeys % n, SIZE(share,2)))
         Res_Char(i_L) % ends=share
      END IF
      CALL Group_Charmm_
      IF(ALLOCATED(shareg)) THEN
         ALLOCATE(Res_Char(i_L) % group(SIZE(shareg)))
         DO n=1,SIZE(shareg)
            ALLOCATE(Res_Char(i_L) % group (n) % g (SIZE(shareg (n) % g)))
            Res_Char(i_L) % group (n) % g=shareg (n) % g
         END DO
      END IF
      CALL Ends_Charmm_
      CALL Print_Errors()
    END SUBROUTINE Charmm_
    SUBROUTINE Ends_Charmm_
      INTEGER :: count,n,n_extra
      CHARACTER(len=max_char) :: line1,line2
      CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: new_bonds
      
      IF(.NOT. ALLOCATED(Res_Char(i_L) % bonds)) RETURN
      ALLOCATE(new_bonds(2,SIZE(Res_Char(i_L) % bonds,2)))
      n_extra=-1
      count=0
      DO n=1,SIZE(Res_Char(i_L) % bonds,2)
         line1=Res_Char(i_L) % bonds(1,n)
         line2=Res_Char(i_L) % bonds(2,n)
         IF(MY_FAM('+',line1) .OR. MY_FAM('+',line2)) THEN
            n_extra=n
         ELSE
            count=count+1
            new_bonds(1,count)=line1
            new_bonds(2,count)=line2
         END IF
      END DO
      IF(SIZE(Res_Char(i_L) % bonds,2)-count > 1) STOP
      
      IF(n_extra > 0) THEN
         IF(.NOT. ALLOCATED(Res_Char(i_L) % ends)) THEN
            n=n_extra
            line1=TRIM(Res_Char(i_L) % bonds(1,n))
            line2=TRIM(Res_Char(i_L) % bonds(2,n))
            IF(MY_FAM('+',line1)) THEN 
               CALL CHR_Delete('+',line1)
               line1=TRIM(line1)
            END IF
            IF(MY_FAM('+',line2)) THEN 
               CALL CHR_Delete('+',line2)
               line2=TRIM(line2)
            END IF
            ALLOCATE(Res_Char(i_L) % ends(2,1))
            Res_Char(i_L) % ends(1,1)=line2
            Res_Char(i_L) % ends(2,1)=line1
         END IF
         DEALLOCATE(Res_Char(i_L) % bonds)
         ALLOCATE(Res_Char(i_L) % bonds(2,count))
         Res_Char(i_L) % bonds=new_bonds
      END IF
      DEALLOCATE(new_bonds)
    END SUBROUTINE Ends_Charmm_

    SUBROUTINE Group_Charmm_
      CHARACTER(len=max_char), DIMENSION(:), POINTER :: cvect=>NULL()
      CHARACTER(len=max_char) :: line
      INTEGER :: c,n,m,nword,nshare

      IF(ALLOCATED(shareg)) DEALLOCATE(shareg)
      c=0
      DO n=1,SIZE(Ri % line)
         nword=MYParse_(Ri % line(n))
         IF(MY_FAM('grou',strngs(1))) THEN
            c=c+1
         END IF
      END DO
      IF(c == 0) RETURN
      ALLOCATE(shareg(c))
      IF(.NOT. Node_()) STOP
      c=0
      DO n=1,SIZE(Ri % line)
         nword=MYParse_(Ri % line(n))
         IF(MY_FAM('grou',strngs(1))) THEN
            IF(c /= 0) THEN
               nshare=Node__Size()
               IF(nshare == 0) CYCLE
               ALLOCATE(shareg(c) % g(nshare))
               m=0
               DO WHILE(Node__Pop(cvect))
                  m=m+1
                  shareg (c) % g (m)=cvect(1)
               END DO
               IF(.NOT. Node_()) STOP
            END IF
            c=c+1
         ELSE IF(MY_FAM('atom',strngs(1))) THEN
            line=ST_Concat(2,strngs)
            CALL Node__push(line)
         END IF
      END DO
      nshare=Node__Size()
      IF(nshare /= 0) THEN
         ALLOCATE(shareg(c) % g(nshare))
         n=0
         DO WHILE(Node__Pop(cvect))
            n=n+1
            shareg (c) % g (n)=cvect(1)
         END DO
      END IF
      IF(ASSOCIATED(cvect)) DEALLOCATE(cvect)
      
    END SUBROUTINE Group_Charmm_
    SUBROUTINE Topol_Charmm_
      INTEGER :: n,m
      INTEGER :: nword,nshare
      CHARACTER(len=max_char) :: line
      LOGICAL :: ok
      CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: cvect0
      CHARACTER(len=max_char), DIMENSION(:), POINTER :: cvect=>NULL()
  
      IF(Ikeys % Method == 0) THEN
         ALLOCATE(cvect0(Ikeys % n))
      ELSE
         ALLOCATE(cvect0(1))
      END IF
      IF(ALLOCATED(share)) DEALLOCATE(share)
      IF(.NOT. Node_()) STOP
      DO n=1,SIZE(Ri % line)
         nword=MYParse_(Ri % line(n))
         ok=.FALSE.
         DO m=1,SIZE(Ikeys % keys)
            line=TRIM(Ikeys % keys(m))
            IF(MY_FAM(line,strngs(1))) THEN
               ok=.TRUE.
               EXIT
            END IF
         END DO
         IF(ok) THEN

!!$
!!$--- Fix for h-bond donor and accepto in charmm files
!!$
            IF(MY_FAM('acc',strngs(1)) .OR. MY_FAM('don',strngs(1))) THEN
               IF(MOD(nword-1,Ikeys % n) /= 0) THEN
                  CYCLE
               END IF
               IF(Ikeys % n == 1 .AND. nword-1 /= 1) CYCLE
            ELSE 
               IF(MOD(nword-1,Ikeys % n) /= 0) THEN
                  errmsg_f='Wrong number of argument to keyword '//TRIM(strngs(1))//&
                       &' while reading Charmm Topology'
                  CALL Add_Errors(-1, errmsg_f)
               END IF
            END IF

            IF(Ikeys % Method == 0) THEN
               DO m=2,nword,Ikeys % n
                  cvect0=strngs(m:m+Ikeys % n-1)
                  CALL Node__Push(cvect0)
               END DO
            ELSE
               cvect0=ST_Concat(2,strngs)
               CALL Node__Push(cvect0)
            END IF
         END IF
      END DO
      nshare=Node__Size()
      IF(nshare /= 0) THEN
         ALLOCATE(share(Ikeys % n, nshare))
         n=0
         DO WHILE(Node__Pop(cvect))
            n=n+1
            share(:,n)=cvect
         END DO
      END IF
      DEALLOCATE(cvect0)
      IF(ASSOCIATED(cvect)) DEALLOCATE(cvect)
    END SUBROUTINE Topol_Charmm_
