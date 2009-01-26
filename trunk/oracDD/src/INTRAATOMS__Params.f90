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
FUNCTION IndIntraBox_n0_() RESULT(out)
  LOGICAL :: out
  INTEGER :: n,m
  INTEGER, ALLOCATABLE :: via(:)
  out=.FALSE.
  IF(.NOT. IndIntraBox_()) RETURN
  
  CALL IntraPot( Prm % Bonds, Tpg % Bonds, Indx_Bonds, Param_Bonds)
  CALL IntraPot( Prm % Constr, Tpg % Bonds, Indx_Constr, Param_Constr)
  CALL IntraPot( Prm % Angles, Tpg % Angles, Indx_Angles, Param_Angles)
  CALL IntraPot( Prm % Imph, Tpg % Imph, Indx_Imph, Param_Imph)
  out=.TRUE.
END FUNCTION IndIntraBox_n0_
FUNCTION IndIntraBox_n1_() RESULT(out)
  LOGICAL :: out
  INTEGER :: n,m
  out=.FALSE.
  
  IF(.NOT. IndIntraBox_()) RETURN
  
  CALL IntraPot(Prm % Dihed, Tpg % Dihed, Indx_Dihed, Param_Dihed)
  CALL IntraPot14(Tpg % Int14, Indx_Int14, Param_Int14)
  
  out=.TRUE.
END FUNCTION IndIntraBox_n1_

SUBROUTINE IntraPot(Prm, Tpg, Indx, Param)
  TYPE(SystemPrm__Chain) :: Prm(:)
  INTEGER :: Tpg(:,:)
  INTEGER, ALLOCATABLE :: Indx(:,:)
  TYPE(IntraParam), ALLOCATABLE :: Param(:)
  
  INTEGER, ALLOCATABLE :: via(:),ind_o(:)
  INTEGER :: n_Tpg,m_Prm,count0,count1,n,nn,ng,n1,n2,nnn


  m_Prm=SIZE(Prm)
  n_Tpg=SIZE(Tpg,1)

  ALLOCATE(via(n_Tpg),ind_o(m_Prm))

  ind_o=0
  count0=0
  DO nn=1,m_Prm
     n=Prm (nn) % pt
     IF(n < 0) CYCLE
     via=Tpg(:,n)
     n1=COUNT(oks(via(:)))
     IF(n1 == 0) CYCLE
     n2=COUNT(okt(via(:)))
     IF(n1+n2 == n_Tpg) THEN
        count0=count0+1
        ind_o(count0)=nn
     END IF
  END DO
  IF(ALLOCATED(Indx)) DEALLOCATE(Indx)
  IF(ALLOCATED(Param)) DEALLOCATE(Param)
  
  ALLOCATE(Indx(n_Tpg,count0),Param(count0))
  
  DO nnn=1,count0
     nn=ind_o(nnn)
     n=Prm (nn) % pt
     via=Tpg (:,n)
     ng=SIZE(Prm (nn) % g)
     ALLOCATE(Param(nnn) % pot(ng+1))
     Indx(:,nnn)=Intra_Conv(via)
     Param(nnn) % pot(1:ng)=Prm(nn) % g
     Param(nnn) % pot(ng+1)=1.0_8
     IF(COUNT(knwn(via(:)) == -1) /= 0) Param(nnn) % pot(ng+1)=0.0_8
  END DO
!!$  CALL MPI_ALLREDUCE(MPI_IN_PLACE,count0,1,MPI_INTEGER4,MPI_SUM,PI_Comm_Cart,ierr)
END SUBROUTINE IntraPot
SUBROUTINE IntraPot14(Tpg, Indx, Param)
  INTEGER :: Tpg(:,:)
  INTEGER, ALLOCATABLE :: Indx(:,:)
  TYPE(IntraParam), ALLOCATABLE :: Param(:)
  
  INTEGER, ALLOCATABLE :: via(:),ind_o(:)
  INTEGER :: n_Tpg,m_Tpg,count0,n,nn,ng,n1,n2
  
  
  n_Tpg=SIZE(Tpg,1)
  m_Tpg=SIZE(Tpg,2)
  
  ALLOCATE(via(n_Tpg),ind_o(m_Tpg))
  count0=0
  DO n=1,m_Tpg
     via=Tpg(:,n)
     n1=COUNT(oks(via(:)))
     IF(n1 == 0) CYCLE
     n2=COUNT(okt(via(:)))
     IF(n1+n2 == n_Tpg) THEN
        count0=count0+1
        ind_o(count0)=n
     END IF
  END DO
  IF(ALLOCATED(Indx)) DEALLOCATE(Indx)
  IF(ALLOCATED(Param)) DEALLOCATE(Param)
  ALLOCATE(Indx(n_Tpg,count0),Param(count0))
  
  DO nn=1,count0
     n=ind_o(nn)
     via=Tpg (:,n)
     Indx(:,nn)=Intra_Conv(via)
     ALLOCATE(Param(nn) % pot(1))
     Param(nn) % pot(1)=1.0_8
     IF(COUNT(knwn(via(:)) == -1) /= 0) Param(nn) % pot(1)=0.0_8
  END DO
!!$  CALL MPI_ALLREDUCE(MPI_IN_PLACE,count0,1,MPI_INTEGER4,MPI_SUM,PI_Comm_Cart,ierr)
END SUBROUTINE IntraPot14

FUNCTION IndIntraBox_() RESULT(out)
  LOGICAL :: out
  INTEGER :: n,m,count_a_p,mm
  
  
  count_a_p=0
  count_a_p=COUNT(Groupa(Atoms(:) % Grp_No) % knwn == 1)
  
  IF(ALLOCATED(oks)) DEALLOCATE(oks)
  IF(ALLOCATED(okt)) DEALLOCATE(okt)
  
  ALLOCATE(oks(SIZE(Atoms)),okt(SIZE(Atoms)))
  
  count_a_p=0 
  oks=.FALSE.
  okt=.FALSE.
  DO n=1,SIZE(Atoms)
     m=Atoms(n) % knwn
     IF(m == 1) THEN
        count_a_p=count_a_p+1
     END IF
     IF(m == 1) oks(n)=.TRUE.
     IF(m == 2) okt(n)=.TRUE.

  END DO
  out=count_a_p /= 0
  IF(.NOT. out) THEN
     errmsg_f='No Primary Atoms for Intramolecular interaction found in the unit box'
     CALL Add_Errors(-1,errmsg_f)
  END IF
END FUNCTION IndIntraBox_
SUBROUTINE Labelling
  INTEGER :: n,nc,l
  LOGICAL, ALLOCATABLE :: ok_tot(:)
  INTEGER, SAVE :: count_label

  IF(.NOT. IndIntraBox_()) RETURN
  ALLOCATE(ok_tot(SIZE(Atoms)))
  ok_tot=.FALSE.
  count_Label=0
  
  CALL MyIntraPot(Prm % Dihed, Tpg % Dihed)
  CALL MyIntraPot( Prm % Angles, Tpg % Angles)
  CALL MyIntraPot( Prm % Bonds, Tpg % Bonds)
  CALL MyIntraPot( Prm % Imph, Tpg % Imph)

  nc=0
  DO n=1,natom
     l=Intra_t(n)
     IF(Atoms(l) % knwn == 1) THEN
        nc =nc +1
        Intra_p(nc)=l
     END IF
  END DO
     
CONTAINS
  SUBROUTINE MyIntraPot(Prm, Tpg)
  TYPE(SystemPrm__Chain) :: Prm(:)
  INTEGER :: Tpg(:,:)
  INTEGER, ALLOCATABLE :: via(:)
  INTEGER :: n_Tpg,m_Prm,count1,n,nn,ng,n1,n2,nnn,m,mm

  m_Prm=SIZE(Prm)
  n_Tpg=SIZE(Tpg,1)

  ALLOCATE(via(n_Tpg))

  DO nn=1,m_Prm
     n=Prm (nn) % pt
     IF(n < 0) CYCLE
     via=Tpg(:,n)
     n1=COUNT(oks(via(:)))
     IF(n1 == 0) CYCLE
     n2=COUNT(okt(via(:)))
     IF(n1+n2 == n_Tpg) THEN
        DO mm=1,n_Tpg
           m=via(mm)
           IF(.NOT. ok_tot(m)) THEN
              ok_tot(m)=.TRUE.
              count_Label=count_Label+1
              Intra_T(count_Label)=m
           END IF
        END DO
     END IF
  END DO
END SUBROUTINE MyIntraPot

END SUBROUTINE Labelling
