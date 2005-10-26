!!$***********************************************************************
!!$   Time-stamp: <2005-10-26 23:59:46 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Feb 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module/class is part of the program ORAC ----*
MODULE Class_Gcharges
  USE xerror_mod
  USE Class_Gauss_Param, ONLY: Gauss_Para, Param_print=>Print
  TYPE Gauss_Charges
     INTEGER :: label=-1
     INTEGER :: nato,natom
     INTEGER, DIMENSION(:), POINTER :: ind
     REAL(8), DIMENSION(:), POINTER :: chg
     REAL(8), DIMENSION(:), POINTER :: phi
     REAL(8), DIMENSION(:), POINTER :: fx,fy,fz
  END TYPE Gauss_Charges
  PRIVATE :: copy, Zero, Create
CONTAINS
  FUNCTION Init(n,beta) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    INTEGER :: n
    CHARACTER(7) :: beta(n)
    TYPE(Gauss_Charges) :: out
    TYPE(Gauss_Para) :: Param

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n_new
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    Param=Param_Print()
    IF(.NOT. Param % initialized) CALL abort_now('Gauss_Param not init&
         &ialized.')

    out % label = 0
    out % nato =n
    n_new=0
    DO i=1,n
       IF(beta(i) == Param % Type) THEN
          n_new=n_new+1
       END IF
    END DO
    out % natom = n_new
    ALLOCATE(out % ind(n_new))

    n_new=0
    DO i=1,n
       IF(beta(i) == Param % Type) THEN
          n_new=n_new+1
          out % ind(n_new)=i
       END IF
    END DO
    ALLOCATE(out % fx(n_new), out % fy(n_new), out % fz(n_new))
    ALLOCATE(out % chg(n_new))
    ALLOCATE(out % phi(n_new))
  END FUNCTION Init
  FUNCTION copy(in) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
    TYPE(Gauss_Charges) :: in,out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: n
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF( in % label == -1) CALL abort_now('Gauss_Charges not initialize&
         &d.')
    n=in % natom
    IF(n >= 0) THEN
       IF(.NOT. ASSOCIATED(out % ind)) THEN
          ALLOCATE(out % ind(n))
          ALLOCATE(out % chg(n))
          ALLOCATE(out % phi(n))
          ALLOCATE(out % fx(n), out % fy(n), out % fz(n))
       END IF
       out=in
    END IF
  END FUNCTION copy
  FUNCTION create(in) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
    TYPE(Gauss_Charges) :: in,out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: n
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF( in % label == -1) CALL abort_now('Gauss_Charges not initialize&
         &d.')
    n=in % natom
    IF(n >= 0) THEN
       IF(.NOT. ASSOCIATED(out % ind)) THEN
          out % label = 0
          out % nato = in % nato
          out % natom = in % natom
          ALLOCATE(out % ind(n))
          out % ind = in % ind
          ALLOCATE(out % chg(n))
          out % chg = in % chg
          ALLOCATE(out % phi(n))
          ALLOCATE(out % fx(n), out % fy(n), out % fz(n))
       END IF
    END IF
  END FUNCTION create
  FUNCTION Make_copy(in) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
    TYPE(Gauss_Charges) :: in,out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: n
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF( in % label == -1) CALL abort_now('Gauss_Charges not initialize&
         &d.')
    n=in % natom
    IF(n >= 0) THEN
       IF(.NOT. ASSOCIATED(out % ind)) THEN
          ALLOCATE(out % ind(n))
          ALLOCATE(out % chg(n))
          ALLOCATE(out % phi(n))
          ALLOCATE(out % fx(n), out % fy(n), out % fz(n))
       END IF
       out=in
    END IF
  END FUNCTION Make_copy
  SUBROUTINE Zero (in)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
    TYPE(Gauss_Charges) :: in
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: n
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF( in % label == -1) CALL abort_now('Gauss_Charges not initialize&
         &d.')

    n=in % natom
    in % fx(1:n) = 0.0D0
    in % fy(1:n) = 0.0D0
    in % fz(1:n) = 0.0D0
    in % phi(1:n) = 0.0D0
    in % chg(1:n) = 0.0D0

  END SUBROUTINE Zero
  FUNCTION Get_Density(in) RESULT (out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    TYPE(Gauss_Charges) :: out,in
!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: n
    LOGICAL, SAVE :: first_call=.TRUE.
    TYPE(Gauss_Para) :: den
    REAL(8) :: fact1,fact2,fact3
 
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF( in % label == -1) CALL abort_now('Gauss_Charges not initialize&
         &d.')

    den=Param_Print()

    IF(den % initialized) THEN
       out = Copy(in)

       fact1=den % n_plus/DBLE(out % natom)
       fact2=den % n_minus/DBLE(out % natom)
       fact3=den % efact/den % kt

       out % chg = fact3 * out % chg * in % phi
       out % chg = fact1 * DEXP(-out % chg) - fact2 * DEXP(out % chg)
    ELSE
       CALL abort_now('Gaussian parameters not initialized.')
    END IF

  END FUNCTION Get_Density

  FUNCTION Get_Phi(in) RESULT (out)
!!$======================== DECLARATIONS ================================*
    USE Module_Fourier
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    TYPE(Gauss_Charges) :: out,in
    CHARACTER(80) :: label

    TYPE(Fourier_Out) :: Pme_Out

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: n,i,ii
 
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    
    n=in % nato

    CALL Modify_Charges(in % chg, in % ind)

    Pme_Out=Do_Fourier()

    out=create(in)

    n=in % natom
    DO ii=1,n
       i=in % ind(ii)
       out % fx(ii) = Pme_Out % fx(i)
       out % fy(ii) = Pme_Out % fy(i)
       out % fz(ii) = Pme_Out % fz(i)
       out % phi(ii)= Pme_Out % phi(i)
    END DO
  END FUNCTION Get_Phi
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Class_Gcharges
