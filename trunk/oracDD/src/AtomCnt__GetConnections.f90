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
SUBROUTINE AtomCnt__GetConnections
  INTEGER :: Res_No,n,m,i,j,i_N,i_P,i_F,i1,i2,ii1,ii2,count_a,count_extra,m1,o,o1
  INTEGER :: extra_Bond(2),new_bond(2),nind_x,n1,m2,i_2,ip
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Res_Bonds
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
  INTEGER, DIMENSION(:), POINTER :: bonds
  CHARACTER(len=max_char) :: label0,label1,lab0
  LOGICAL :: end_of_list,ok
  
  Res_No=0
  IF(.NOT. Node_()) STOP
  DO n=1,SIZE(Secondary)
     IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
     DO m=1,SIZE(Secondary(n) % line)
        Res_No=Res_No+1
        i_P=-1
        i_N=-1
        IF(m /= 1) i_P=inds(n) % i (m-1)
        IF(m /= SIZE(Secondary(n) % line)) i_N=inds(n) % i (m+1)
        i_F=inds(n) % i (m) 
!!$
!!$--- Add bonds already in the residue i_F list
!!$

        DO i=1,SIZE(App_Char(i_F) % bonds,2)
           DO j=1,SIZE(App_Char(i_F) % bonds,1)
              label0=TRIM(App_Char(i_F) % bonds(j,i))
!!$--- Look for '1' or '2' at the beginning of the atom label
              IF(ICHAR(label0(1:1)) == 49 .OR. ICHAR(label0(1:1)) == 50) EXIT
              new_bond(j)=AtomCnt__Find(Res_No,label0)
           END DO
           label0=TRIM(App_Char(i_F) % bonds(1,i))
           label1=TRIM(App_Char(i_F) % bonds(2,i))
!!$--- Look for '1' or '2' at the beginning of the atom label
           IF((ICHAR(label0(1:1)) == 49 .AND. ICHAR(label1(1:1)) == 50)&
                & .OR. (ICHAR(label0(1:1)) == 50 &
                &.AND. ICHAR(label1(1:1)) == 49)) THEN
              m2=-1
              DO n1=1,SIZE(IndPa)
                 IF(IndPa(n1) % one == m) THEN
                    m2=IndPa(n1) % two
                 END IF
              END DO
              m1=-1
              DO n1=1,SIZE(IndPa)
                 IF(IndPa(n1) % two == m) THEN
                    m1=IndPa(n1) % one
                 END IF
              END DO
              IF(m2 /= -1) THEN
                 DO j=1,SIZE(App_Char(i_F) % bonds,1)
                    label0=TRIM(App_Char(i_F) % bonds(j,i))
                    label1=label0(2:)
                    IF(ICHAR(label0(1:1)) == 49) THEN
                       new_bond(j)=AtomCnt__Find(Res_No,label1)
                    ELSE IF(ICHAR(label0(1:1)) == 50) THEN
                       i_2=inds(n) % i (m2) 
                       new_bond(j)=AtomCnt__Find(m2,label1)
                    END IF
                 END DO
                 CALL Node__Push(new_bond)
              ELSE IF(m1 == -1) THEN
                 WRITE(lab0,'(i4)') m
                 label0=TRIM(App_Char(i_F) % bonds(1,i))
                 label1=TRIM(App_Char(i_F) % bonds(2,i))
                 errmsg_f='Inter-residue connection '//TRIM(label0)&
                      &//' '//TRIM(label1)//' not found. Residue No. '&
                      &//TRIM(lab0)//' not on the list of linked residues.'
                 CALL Add_Errors(-1,errmsg_f)
                 CALL Print_Errors()
              END IF
           ELSE
              CALL Node__Push(new_bond)
           END IF
        END DO

!!$
!!$--- Add extra bonds if connections to other residues exist
!!$
        count_extra=0
        IF(ALLOCATED(App_Char(i_F)%ends)) THEN
           IF((.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(1,1))))&
                & .AND. i_P == -1) THEN 
              errmsg_f='While creating atom topology, residue '//&
                   &TRIM(App_Char(i_F) % Type)//' appers to be a&
                   &t the beginning of & 
                   &the sequence, but a dangling bond was found.'
              CALL Add_errors(-1,errmsg_f)
              CALL Print_Errors()
           ELSE IF((.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(2&
                &,1)))) .AND. i_N == -1) THEN 
              errmsg_f='While creating atom topology, residue '//&
                   &TRIM(App_Char(i_F) % Type)//' appers to be at the end of &
                   &the sequence, but a dangling bond was found.'
              CALL Add_errors(-1,errmsg_f)
              CALL Print_Errors()
           END IF
           IF(.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(2,1)))) THEN
              extra_bond(1)=AtomCnt__Find(Res_No,App_Char(i_F) % ends(2,1))
              extra_bond(2)=AtomCnt__Find(Res_No+1,App_Char(i_N) % ends(1,1))
              count_extra=1
              CALL Node__Push(extra_bond)
           END IF
        END IF
     END DO
  END DO

!!$
!!$--- Get the connections for each atom
!!$              

  ip=Node__Size()
  ALLOCATE(Res_bonds(2,ip))
  ALLOCATE(ind_x(SIZE(AtomCnts)))
  ind_x=0
  ip=0
  DO WHILE(Node__Pop(bonds))
     ip=ip+1
     i1=Bonds(1)
     i2=Bonds(2)
     Res_Bonds(1,ip)=i1
     Res_Bonds(2,ip)=i2
     ind_x(i1)=ind_x(i1)+1
     ind_x(i2)=ind_x(i2)+1
  END DO
  DO i=1,SIZE(AtomCnts)
     ALLOCATE(AtomCnts(i) % cnt(ind_x(i)))
  END DO
  ind_x=0
  DO i=1,SIZE(Res_Bonds,2)
     i1=Res_Bonds(1,i)
     i2=Res_Bonds(2,i)
     ind_x(i1)=ind_x(i1)+1
     ind_x(i2)=ind_x(i2)+1
     AtomCnts(i1) % cnt(ind_x(i1))=i2
     AtomCnts(i2) % cnt(ind_x(i2))=i1
  END DO
  DEALLOCATE(ind_x)
  DEALLOCATE(Res_bonds)
END SUBROUTINE AtomCnt__GetConnections
