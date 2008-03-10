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

MODULE NUMERIC_KINDS
   !NAMED CONSTANTS USED TO DEFINE VARIABLE TYPES
   INTEGER, PARAMETER :: L1K = 1        !1 byte integer
   INTEGER, PARAMETER :: L2K = 2        !2 byte integer
   INTEGER, PARAMETER :: L4K = 4        !4 byte integer

   INTEGER, PARAMETER :: I1K = SELECTED_INT_KIND(2)                        !1 byte integer
   INTEGER, PARAMETER :: I2K = SELECTED_INT_KIND(4)                        !2 byte integer
   INTEGER, PARAMETER :: I4K = SELECTED_INT_KIND(9)                        !4 byte integer
   INTEGER, PARAMETER :: R4K = KIND(1.0)                                   !single precision
   INTEGER, PARAMETER :: R8K = SELECTED_REAL_KIND(2*PRECISION(1.0_R4K))    !double precision
END MODULE NUMERIC_KINDS
MODULE MATH_USEFUL_DP
   !Some useful mathematical constants in double precision
   USE NUMERIC_KINDS

   REAL(KIND=R8K) zero,one,two,three,four,five,six,seven,eight,nine,ten
   PARAMETER (zero=0.0D+00,one=1.0D+00,two=2.0D+00,three=3.0D+00)
   PARAMETER (four=4.0D+00,five=5.0D+00,six=6.0D+00,seven=7.0D+00)
   PARAMETER (eight=8.0D+00,nine=9.0D+00,ten=10.0D+00)

   REAL(KIND=R8K) half,third
   PARAMETER (HALF=0.5D+00,THIRD=one/three)

   REAL(KIND=R8K) PI,TWOPI,FOURPI,HALFPI
   PARAMETER (PI=3.141592653589793238D+00)
   PARAMETER (HALFPI=HALF*PI,TWOPI=TWO*PI,FOURPI=FOUR*PI)

   !Angle conversion
   REAL(KIND=R8K) RAD2DEG,DEG2RAD
   PARAMETER (RAD2DEG=180.0D+00/PI)  !Convert radians to degress
   PARAMETER (DEG2RAD=PI/180.0D+00)  !Convert degrees to radians

   !Time conversion (based on sidereal year)
   REAL(KIND=R8K) M2H,H2M,H2D,D2H
   REAL(KIND=R8K) S2H,H2S,S2D,D2S,S2Y,Y2S
   REAL(KIND=R8K) M2Y,Y2M,H2Y,Y2H,D2Y,Y2D,D2M,M2D
   PARAMETER (H2S=3600.0D+00,S2H=one/H2S)    !second <-> hour
   PARAMETER (D2S=24.0D+00*H2S,S2D=one/D2S)  !second <-> day
   PARAMETER (Y2S=3.155815D+07,S2Y=one/Y2S)  !second <-> year
   PARAMETER (D2M=1440.0D+00,M2D=one/D2M)    !minute <-> day
   PARAMETER (Y2M=Y2S/60.0D+00,M2Y=one/Y2M)  !minute <-> year
   PARAMETER (H2M=60.0D+00,M2H=one/H2M)      !hour <-> minute
   PARAMETER (D2H=24.0D+00,H2D=one/D2H)      !hour <-> day
   PARAMETER (Y2H=Y2S*S2H,H2Y=one/Y2H)       !hour <-> year

END MODULE MATH_USEFUL_DP
MODULE OTHER_USEFUL
  USE NUMERIC_KINDS
  LOGICAL(L4K) YES,NO
  PARAMETER (YES=.TRUE.,NO=.FALSE.)

  INTEGER(I4K), PARAMETER :: DIAGNOSTIC = 0
  INTEGER(I4K), PARAMETER :: WARNING = 1
  INTEGER(I4K), PARAMETER :: FATAL = -1
  INTEGER(I4K), PARAMETER :: KillProg = -1

END MODULE OTHER_USEFUL
MODULE STRPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  !FLIB Interface

  INTERFACE SP_CATNUM
     MODULE PROCEDURE CATI4N
     MODULE PROCEDURE CATR4N
     MODULE PROCEDURE CATR8N
  END INTERFACE

  INTERFACE SP_GETNUM
     MODULE PROCEDURE STR_GETI4N
     MODULE PROCEDURE STR_GETR4N
     MODULE PROCEDURE STR_GETR8N
  END INTERFACE

  INTERFACE SP_PUTNUM
     MODULE PROCEDURE PUTI4N     
     MODULE PROCEDURE PUTR4N     
     MODULE PROCEDURE PUTR8N     
  END INTERFACE

  INTERFACE SP_WRTNUM
     MODULE PROCEDURE WRTI4N
     MODULE PROCEDURE WRTR4N
     MODULE PROCEDURE WRTR8N
  END INTERFACE
CONTAINS
      SUBROUTINE WRTI4N(i4n,fmat, str, iflg)
!     Attempt to write an integer number to STR using
!     the format specifier FMAT.
!
!     INPUT:
!         I4N: an integer                              [I4]
!        FMAT: integer format specifier such as '(I4)' [A*]
!
!     OUTPUT:
!        STR: character variable                       [A*]
!       IFLG: status flag                              [I4]
!             -1: task failed.
!              0: task successful
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 13-NOV-1997
!     REVISIONS HISTORY:
!        18-NOV-1997 RD STEWART
!           Changed name from PUTSTR to I4NSTR
!        01-JAN-1999 RD STEWART
!           Changed comments to Fortran 90 style and
!           renamed routine from STRI4N to STR_WRTI4N.
!           Also,USE STRPAK to declare i/o variables.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      INTEGER(I4K) I4N,IFLG,IOCHK,ILEN,I
      CHARACTER(LEN=*) STR,fmat
      INTENT(in) :: i4n,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

!     Erase any information in STR.
      ilen=LEN(str)
      DO i=1,ilen
        str(i:i)=' '
      ENDDO

!     Attempt to write I4N to STR.
      WRITE(STR,fmt=FMAT,IOSTAT=iochk) I4N
      IF (iochk.NE.0) THEN
!       Invalid operation
        iflg=-1
      ELSE
!       Success!
        iflg=0
      ENDIF

      RETURN
    END SUBROUTINE WRTI4N

      SUBROUTINE WRTR4N(r4n,fmat, str, iflg)
!     Attempt to write a single precision number to STR using
!     the format specifier FMAT.
!
!     INPUT:
!         R4N: a single precision number                 [R4]
!        FMAT: integer format specifier such as '(F8.3)' [A*]
!
!     OUTPUT:
!        STR: character variable                         [A*]
!       IFLG: status flag                                [I4]
!             -1: task failed.
!              0: task successful
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 13-NOV-1997
!     REVISIONS HISTORY:
!        18-NOV-1997 RD STEWART
!           Changed name from PUTSTR to I4NSTR
!        01-JAN-1999 RD STEWART
!           Changed comments to Fortran 90 style, renamed
!           routine from STRR4N to STR_WRTR4N.  Also,
!           use NUMERIC_KINDS to declare i/o variables.
!
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      REAL(R4K) R4N
      INTEGER(I4K) IFLG,IOCHK,ILEN,I
      CHARACTER(LEN=*) STR,fmat
      INTENT(in) :: r4n,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

!     Erase any information in STR.
      ilen=LEN(str)
      DO i=1,ilen
        str(i:i)=' '
      ENDDO

!     Attempt to write I4N to STR.
      WRITE(STR,fmt=FMAT,IOSTAT=iochk) R4N
      IF (iochk.NE.0) THEN
!       Invalid operation
        iflg=-1
      ELSE
!       Success!
        iflg=0
      ENDIF

      RETURN
    END SUBROUTINE WRTR4N

      SUBROUTINE WRTR8N(r8n,fmat, str, iflg)
!     Attempt to write a double precision number to STR using
!     the format specifier FMAT.
!
!     INPUT:
!         R8N: a single precision number                 [R8]
!        FMAT: integer format specifier such as '(F8.3)' [A*]
!
!     OUTPUT:
!        STR: character variable                         [A*]
!       IFLG: status flag                                [I4]
!             -1: task failed.
!              0: task successful
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 13-NOV-1997
!     REVISIONS HISTORY:
!        18-NOV-1997 RD STEWART
!           Changed name from PUTSTR to I4NSTR
!        01-JAN-1999 RD STEWART
!           Changed comments to Fortran 90 style, renamed
!           routine from STRR8N to STR_WRTR8N.  Also,
!           USE STRPAK to declare i/o variables.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      REAL(R8K) R8N
      INTEGER(I4K) IFLG,IOCHK,ILEN,I
      CHARACTER(LEN=*) STR,fmat
      INTENT(in) :: r8n,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

!     Erase any information in STR.
      ilen=LEN(str)
      DO i=1,ilen
        str(i:i)=' '
      ENDDO

!     Attempt to write I4N to STR.
      WRITE(STR,fmt=FMAT,IOSTAT=iochk) R8N
      IF (iochk.NE.0) THEN
!       Invalid operation
        iflg=-1
      ELSE
!       Success!
        iflg=0
      ENDIF

      RETURN
    END SUBROUTINE WRTR8N
      SUBROUTINE CATI4N(i4n,fmat,div, str)
!     Concatenate the character variable DIV to STR and then
!     write the integer number number I4N to the
!     end of STR // DIV.  Trailing blanks in STR are ignored,
!     but all characters in DIV are significant.  FMAT is the
!     same abbreviated format specifier used in the PUTI4N,
!     PUTR4N, and PUTR8N routines.
!
!     INPUT:
!        I4N: an integer number                      [I4]
!        DIV: a character variable                   [A*]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!
!     OUTPUT:
!        STR: modified string of characters          [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       02-JAN-1999 RD STEWART
!         Change comments to Fortran 90 style.  To simplify
!         future maintenance and upgrades to the CATXXX
!         series of routines, the CATR4N and CATI4N
!         are dummy routines that simply call the CATR8N
!         routine.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) div,str,fmat
      INTEGER(I4K) i4n
      REAL(R8K) r8n
      INTENT(in) :: i4n,fmat,div
      INTENT(inout) :: str

      r8n=i4n
      CALL CATR8N(r8n,fmat,div, str)

      RETURN
    END SUBROUTINE CATI4N

      SUBROUTINE CATR4N(r4n,fmat,div, str)
!     Concatenate the character variable DIV to STR and then
!     write the single precision number R4N to the
!     end of STR // DIV.  Trailing blanks in STR are ignored,
!     but all characters in DIV are significant.  FMAT is the
!     same abbreviated format specifier used in the PUTI4N,
!     PUTR4N, and PUTR8N routines.
!
!     INPUT:
!        R4N: a single precision number              [R4]
!        DIV: a character variable                   [A*]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!
!     OUTPUT:
!        STR: modified string of characters          [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       02-JAN-1999 RD STEWART
!         Change comments to Fortran 90 style.  To simplify
!         future maintenance and upgrades to the CATXXX
!         series of routines, the CATR4N and CATI4N
!         are dummy routines that simply call the CATR8N
!         routine.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) div,str,fmat
      REAL(R4K) r4n
      REAL(R8K) r8n
      INTENT(in) :: r4n,fmat,div
      INTENT(inout) :: str

      r8n=r4n
      CALL CATR8N(r8n,fmat,div, str)

      RETURN
    END SUBROUTINE CATR4N

      SUBROUTINE CATR8N(r8n,fmat,div, str)
!     Concatenate the character variable DIV to STR and then
!     write the double precision number R8N to the
!     end of STR // DIV.  Trailing blanks in STR are ignored,
!     but all characters in DIV are significant.  FMAT is the
!     same abbreviated format specifier used in the PUTI4N,
!     PUTR4N, and PUTR8N routines.
!
!     INPUT:
!        R8N: a double precision number              [R8]
!        DIV: a character variable                   [A*]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!
!     OUTPUT:
!        STR: modified string of characters          [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       02-JAN-1999 RD STEWART
!         Change comments to Fortran 90 style and update
!         to reflect changes in other FLIB routines.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) div,str,fmat
      INTEGER(I4K) IFLG,I1,I2
      REAL(R8K) r8n
      INTENT(in) :: r8n,fmat,div
      INTENT(inout) :: str

      I1=MAX(1,LENTRIM(str))+1               !ignore trailing comments in STR
      I2=LEN(div)+I1                         !all characters in DIV are significant
      CALL STR_FILL(' ',str(i1:))
      CALL STR_COPY(i1,div, str)             !copy DIV to string
      CALL PUTR8N(r8n,i2,-1,fmat, str,iflg)  !write number to end of STR

      RETURN
    END SUBROUTINE CATR8N
      SUBROUTINE STR_GETI4N(str, i4n,iflg)
!     Attempt to convert the sequence of characters in STR
!     into an INTEGER number.
!
!     INPUT:
!         STR: string containing an integer number   [A*]
!
!     OUTPUT:
!         I4N: an integer                            [I4]
!        IFLG: status flag                           [I4]
!              -1: task failed.
!               0: task successful
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-NOV-1997
!     REVISIONS HISTORY:
!       01-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style and
!         renamed from STRI4N to STR_GETI4N
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER(I4K) iflg,i4n
      REAL(R8K) r8n
      INTENT(in) :: str
      INTENT(out) :: i4n,iflg

      CALL STR_GETR8N(str, r8n,iflg)
      i4n=NINT(r8n)  !round to nearest integer

      RETURN
    END SUBROUTINE STR_GETI4N

      SUBROUTINE STR_GETR4N(str, r4n,iflg)
!     Attempt to convert the sequence of characters in STR
!     into an single precision number.
!
!     INPUT:
!         STR: string containing r4n number    [A*]
!
!     OUTPUT:
!         R4N: a single precision number       [R4]
!        IFLG: status flag                     [I4]
!              -1: task failed.
!               0: task successful
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-NOV-1997
!     REVISIONS HISTORY:
!       01-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style and
!         renamed from STRR4N to STR_GETR4N
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER iflg
      REAL(R4K) r4n
      REAL(R8K) r8n
      INTENT(in) :: str
      INTENT(out) :: r4n,iflg

      CALL STR_GETR8N(str, r8n,iflg)
      r4n=r8n

      RETURN
    END SUBROUTINE STR_GETR4N

      SUBROUTINE STR_GETR8N(str, r8n,iflg)
!     Attempt to convert the sequence of characters in STR
!     into a double precision number.
!
!     INPUT:
!         str: string containing R8N number       [A*]
!
!     OUTPUT:
!         R8N: a double precision number          [R8]
!        IFLG: status flag                        [I4]
!               0: task successful
!              -1: task failed.
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-NOV-1997
!     REVISIONS HISTORY:
!       01-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style and
!         renamed from STRR8N to STR_GETR8N
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      CHARACTER(LEN=7) fmat
      INTEGER iflg,ilen
      REAL(R8K) R8N
      INTENT(in) :: str
      INTENT(out) :: r8n,iflg

      ilen=LENTRIM(str)
      IF (ilen.LT.1) THEN
        iflg=-1
      ELSE
        fmat='(F     '
        WRITE(fmat(3:7),'(I2,A)') ILEN,'.0)'
        READ(str,fmt=fmat,iostat=iflg) r8n
        IF (iflg.NE.0) iflg=-1
      ENDIF

      RETURN
    END SUBROUTINE STR_GETR8N

      SUBROUTINE PUTI4N(i4n,loc,nj,fmat, str,iflg)
!     Write an integer number to STR at position LOC
!     using left (NJ=-1), center (NJ=0), or right (NJ=1)
!     justification.  If NJ=2, the number is aligned so
!     that the decimal point (if it existed) would be
!     located at position=LOC.  FMAT is an abbreviated
!     form of the standard format specifiers used in
!     Fortran (see below).
!
!     INPUT:
!         I4N: an integer number                     [I4]
!         LOC: location in STR to put I4N            [I4]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!          NJ: type of justification to use          [I4]
!              -1: left justify at LOC
!               0: center number about LOC
!               1: right justify number at LOC
!               2: place # so that decimal place
!                  is located at LOC.
!
!     OUTPUT:
!         STR: character variable                   [A*]
!        IFLG: status flag                          [I4]
!               0: task successful.
!              -1: undecodable format specifier
!              -2: bad trouble in WRTR8N
!              -3: bad trouble in PUTR8N
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       13-JUL-1998 RD STEWART
!         Changed I=LOC-IL+1 to I=LOC-IL so that justification is
!         consistent with PUTR4N and PUTR8N routines.
!       02-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style.  Converted
!         to a dummy routine for the PUTR8N routine to so
!         that it is easier to maintain/updated the PUTXXN
!         series of routines.
!
!     COMMENTS: This routine is a dummy routine that calls the
!               PUTR8N routine.
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) i4n,loc,nj,iflg
      REAL(R8K) r8n
      INTENT(in) :: i4n,loc,nj,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

      r8n=i4n
      CALL PUTR8N(r8n,loc,nj,fmat, str,iflg)

      RETURN
    END SUBROUTINE PUTI4N

      SUBROUTINE PUTR4N(r4n,loc,nj,fmat, str,iflg)
!     Write an single precision number, R4N, to STR at
!     position LOC using left (NJ=-1), center (NJ=0),
!     or right (NJ=1) justification.  If NJ=2, the
!     number is aligned so that the decimal point is
!     located at position=LOC.  FMAT is an abbreviated
!     form of the standard format specifiers used in
!     Fortran (see below).
!
!     INPUT:
!         R4N: a single precision number             [R4]
!         LOC: location in STR to put I4N            [I4]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!          NJ: type of justification to use          [I4]
!              -1: left justify at LOC
!               0: center number about LOC
!               1: right justify number at LOC
!               2: place # so that decimal place
!                  is located at LOC.
!
!     OUTPUT:
!         STR: character variable                   [A*]
!        IFLG: status flag                          [I4]
!               0: task successful.
!              -1: undecodable format specifier
!              -2: bad trouble in WRTR8N
!              -3: bad trouble in PUTR8N
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       13-JUL-1998 RD STEWART
!         Bug fix: invalid NMAT for some values of FMAT
!       02-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style.  Converted
!         to a dummy routine for the PUTR8N routine to so
!         that it is easier to maintain/updated the PUTXXN
!         series of routines.
!
!     COMMENTS: This routine is a dummy routine that calls the
!               PUTR8N routine.
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) loc,nj,iflg
      REAL(R4K) r4n
      REAL(R8K) r8n
      INTENT(in) :: r4n,loc,nj,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

      r8n=r4n
      CALL PUTR8N(r8n,loc,nj,fmat, str,iflg)

      RETURN
    END SUBROUTINE PUTR4N

      SUBROUTINE PUTR8N(r8n,loc,nj,fmat, str,iflg)
!     Write an double precision number, r8n, to STR at
!     position LOC using left (NJ=-1), center (NJ=0),
!     or right (NJ=1) justification.  If NJ=2, the
!     number is aligned so that the decimal point is
!     located at position=LOC.  FMAT is an abbreviated
!     form of the standard format specifiers used in
!     Fortran (see below).
!
!     INPUT:
!         R8N: a single precision number             [R4]
!         LOC: location in STR to put I4N            [I4]
!        FMAT: format specifier                      [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer format
!              Gn = F or E format, as appropriate.
!          NJ: type of justification to use          [I4]
!              -1: left justify at LOC
!               0: center number about LOC
!               1: right justify number at LOC
!               2: place # so that decimal place
!                  is located at LOC.
!
!     OUTPUT:
!         STR: character variable                   [A*]
!        IFLG: status flag                          [I4]
!               0: task successful.
!              -1: undecodable format specifier
!              -2: bad trouble in WRTR4N
!              -3: bad trouble in PUTR4N
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       13-JUL-1998 RD STEWART
!         Bug fix: invalid NMAT for some values of FMAT
!       02-JAN-1999 RD STEWART
!         Convert comments to Fortran 90 style and updated
!         calls to reflect changes in other FLIB routines.
!
!     COMMENTS: This routine is not at fast as directly printing
!               a number to a device, but it gives a lot more
!               control over the placement of the number (useful
!               for building tables and nice looking output).
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      CHARACTER(LEN=80) NUM
      CHARACTER(LEN=10) NMAT
      REAL(R8K) r8n
      LOGICAL mimick
      INTEGER loc,nj,iflg
      INTEGER I,IP,IL
      INTENT(in) :: r8n,loc,nj,fmat
      INTENT(inout) :: str
      INTENT(out) :: iflg

!     Erase workspace.
      CALL STR_FILL(' ',num)
      CALL FMT_D4RN(fmat, nmat,iflg)
      IF (iflg.LT.0) THEN
        !invalid format specification string
        iflg=-1
        RETURN
      ENDIF
      mimick=(iflg.EQ.1)   !Mimick integer format

!     Transfer r8n to NUM
      CALL WRTR8N(r8n,nmat, num,iflg)
      IF (iflg.NE.0) THEN
        IFLG=-2
        RETURN
      ELSE
        CALL STRIP(num)  !remove blanks
        IP=INDEX(num,'.')
        IF (mimick.AND.(ip.GT.0)) THEN
!         Remove decimal point, and all digits to right, to
!         mimick integer format.
          DO i=ip,LEN(num)
            num(i:i)=' '
          ENDDO
        ENDIF
        IL=LENTRIM(num)
        IF (IL.LT.1) THEN
          IFLG=-3
          RETURN
        ENDIF
      ENDIF

!     Copy NUM into STR with left, right, or center justification.
      IF (nj.EQ.2) THEN
!       Place decimal point at LOC.
        I=LOC-IP+1
        CALL SP_PUTSTR(num(:il),i,-1, str,iflg)
      ELSE
        CALL SP_PUTSTR(num(:il),loc,nj, str,iflg)
      ENDIF

      RETURN
    END SUBROUTINE PUTR8N

      SUBROUTINE CHR_DELETE(ch,str)
!     Delete all occurances of character CH in STR.
!
!     INPUT:
!        ch: a character                  [A1]
!       str: character variable           [A*]
!
!     OUTPUT:
!       str: modified character variable  [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Routine renamed
!          from DELCH to CHR_ERASE.
!
!     COMMENTS: The length of STR on return, minus trailing blanks,
!               will always be less than or equal to the original
!               length minus trailing blanks.
!
      CHARACTER(LEN=1) ch,buf
      CHARACTER(LEN=*) str
      INTEGER ILEN,I,J
      INTENT(in) :: ch
      INTENT(inout) :: str

      ILEN=LEN(str)
      j=0
      DO i=1,ILEN
        IF (str(i:i).NE.ch) THEN
!         Keep this character
          j=j+1
          buf=str(i:i)
          str(j:j) = buf
        ENDIF
      ENDDO
      DO i=j+1,ILEN
        str(i:i)=' '
      ENDDO

      RETURN
    END SUBROUTINE CHR_DELETE

      FUNCTION CHR_ILC(ch)
!     Return CHR_ILC as TRUE if CH is a lowercase letter.
!
!     INPUT:
!        ch: a character [A1]
!
!     OUTPUT:
!       CHR_ILC: is CH a lowercase letter? [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUN-1998
!     REVISIONS HISTORY:
!       29-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISLCL to CHR_ILC
!
!     COMMENTS:
!
      CHARACTER(LEN=1) ch
      INTEGER iasc
      INTENT(in) :: ch
      LOGICAL CHR_ILC

      iasc = ICHAR(ch)
      CHR_ILC = ((iasc.GE.97).AND.(iasc.LE.122))

      RETURN
    END FUNCTION CHR_ILC

FUNCTION CHR_IMAT(ch1,ch2)
!     CHR_IMAT is returned as TRUE if CH1 is an upper-
!     or lowercase match for CH2.
!
!     INPUT:
!        ch1: a character                      [A1]
!        ch2: a character                      [A1]
!
!     OUTPUT:
!       CHR_IMAT: logical variable indicating  [L4]
!                 whether CH1 is a case
!                 insensitive match to CH2
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 07-JUL-1997
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       31-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISMAT to CHR_IMAT.  Removed (trival)
!          option for case sensitive match between CH1 and CH2.
!          Made algorithm more computationally efficient.
!
!     COMMENTS:
!
      CHARACTER(LEN=1) ch1,ch2
      LOGICAL Ok
      INTEGER i1,i2
      INTENT(in) :: ch1,ch2
      LOGICAL CHR_IMAT

      ok=(ch1.EQ.ch2)
      IF (.not.ok) THEN
!       Check for case insensitive match. NOTE: the ASCII
!       code for an uppercase letter can be found by adding
!       32 to the ASCII code for the corresponding
!       lowercase letter.
        i1 = ICHAR(ch1)
        i2 = ICHAR(ch2)
        IF ((i1.GE.65).AND.(i1.LE.90)) THEN
          !ch1 is uppercase
          ok=((i1+32).EQ.i2)
        ELSEIF ((i1.GE.97).AND.(i1.LE.122)) THEN
          !ch1 is lowercase
          ok=((i1-32).EQ.i2)
        ENDIF
      ENDIF
      CHR_IMAT = ok

      RETURN
    END FUNCTION CHR_IMAT

      FUNCTION CHR_IUC(ch)
!     Return CHR_IUC as TRUE if CH is an uppercase letter.
!
!     INPUT:
!        ch: a character                    [A1]
!
!     OUTPUT:
!       CHR_IUC: is CH an uppercase letter? [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUN-1998
!     REVISIONS HISTORY:
!       29-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISLCL to CHR_ILC
!
!     COMMENTS:
      CHARACTER(LEN=1) ch
      INTEGER iasc
      INTENT(in) :: ch
      LOGICAL CHR_IUC

      iasc = ICHAR(ch)
      CHR_IUC = ((iasc.GE.65).AND.(iasc.LE.90))

      RETURN
    END FUNCTION CHR_IUC

      FUNCTION CHR_ILET(ch)
!     Return CHR_ILET as TRUE if CH is an uppercase or
!     lowercase letter.
!
!     INPUT:
!        ch: a character              [A1]
!
!     OUTPUT:
!       CHR_ILET: is CH a letter?     [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-JAN-1995
!     REVISIONS HISTORY:
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       29-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISLET to CHR_ILET
!
!     COMMENTS:
      CHARACTER(LEN=1) ch
      LOGICAL Ok
      INTEGER iasc
      INTENT(in) :: ch
      LOGICAL CHR_ILET

      iasc = ICHAR(ch)
      ok = .false.
      IF ((iasc.GE.65).AND.(iasc.LE.90)) THEN
!       Uppercase Character
        ok = .true.
      ELSEIF ((iasc.GE.97).AND.(iasc.LE.122)) THEN
!       Lowercase Character
        ok = .true.
      ENDIF
      CHR_ILET = Ok

      RETURN
    END FUNCTION CHR_ILET

      FUNCTION CHR_ILON(ch)
!     Return CHR_ILON as TRUE if CH is an upper- or lowercase
!     letter or a number (0 to 9).
!
!     INPUT:
!        ch: a character variable                [A1]
!
!     OUTPUT:
!       CHR_ILON: is CH a letter or number?      [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-JAN-1995
!     REVISIONS HISTORY:
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       29-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISLON to CHR_ILON
!
!     COMMENTS:
!
      CHARACTER(LEN=1) ch
      LOGICAL Ok
      INTEGER iasc
      INTENT(in) :: ch
      LOGICAL CHR_ILON

      iasc = ICHAR(ch)
      ok = .false.
      IF ((iasc.GE.65).AND.(iasc.LE.90)) THEN
!       Uppercase Character
        ok = .true.
      ELSEIF ((iasc.GE.97).AND.(iasc.LE.122)) THEN
!       Lowercase Character
        ok = .true.
      ELSEIF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
!       Number
        ok = .true.
      ENDIF
      CHR_ILON = Ok

      RETURN
    END FUNCTION CHR_ILON

      FUNCTION CHR_INUM(ch)
!     CHR_INUM is returned = TRUE if CH is a number (0 to 9).
!
!     INPUT:
!        ch: a character                       [A1]
!
!     OUTPUT:
!       CHR_INUM: is CH a number?              [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 11-NOV-1993
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       29-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISNUM to CHR_INUM
!
!     COMMENTS:
!
      CHARACTER(LEN=1) ch
      INTEGER iasc
      LOGICAL lnum
      INTENT(in) :: ch
      LOGICAL CHR_INUM

      iasc = ICHAR(ch)
      lnum = .FALSE.
      IF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
        lnum = .TRUE.
      ENDIF
      CHR_INUM = lnum

      RETURN
    END FUNCTION CHR_INUM

      FUNCTION CHR_IPRN(ch)
!     Return CHR_IPRN as TRUE if CH is a printable
!     character, i.e., a character with an ASCII
!     code > 31.
!
!     INPUT:
!             ch: a character                  [A1]
!
!     OUTPUT:
!       CHR_IPRN: is CH a printable character? [L4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       30-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISPRN to CHR_IPRN.  Also, blanks
!          spaces are now counted as a printble character.
!
!     COMMENTS:
      CHARACTER(LEN=1) ch
      INTENT(in) :: ch
      LOGICAL CHR_IPRN

      CHR_IPRN=(ICHAR(ch).GE.32)

      RETURN
    END FUNCTION CHR_IPRN

      SUBROUTINE CHR_RMI(co,cn, str)
!     Search for multiple, contiguous instances of CO in STR
!     and replace them with a single instance of CN.
!
!     INPUT:
!        CO: old character                [A1]
!        CN: new character                [A1]
!       STR: character variable           [A*]
!
!     OUTPUT:
!       STR: modified character variable  [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       28-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from SARMI to CHR_RMI
!
!     COMMENTS: The length of STR on return is always less
!               than or equal to the original.
!
      CHARACTER(LEN=1) cn,co,buf
      CHARACTER(LEN=*) str
      INTEGER Tmax,ilen,i
      LOGICAL srt
      INTENT(in) :: cn,co
      INTENT(inout) :: str

!     Compute length of string
      Tmax=LEN(str)

!     Srt indicates the start of a replacement sequence.
      srt=.true.
      ilen = 0
      DO i=1,Tmax
        IF (str(i:i).EQ.co) THEN
          IF (srt) THEN
!           Perform (macro) substitution.
            srt=.false.
            ilen = ilen+1
            str(ilen:ilen) = cn
          ENDIF
        ELSE
!         Add current character to str at position ilen.
!         srt is set  true so that the next occurance of a
!         character in co will begin the removal procedure.
!         NOTE: buf is used to avoid self-referencing a
!         character variable a no-no in F77.
          ilen = ilen+1
          buf = str(i:i)
          str(ilen:ilen) = buf
          srt = .true.
        ENDIF
      ENDDO

!     Fill remainder of str with blanks.
      DO i=ilen+1,Tmax
        str(i:i) = ' '
      ENDDO

      RETURN
    END SUBROUTINE CHR_RMI

      SUBROUTINE CHR_SAR(co,cn, str)
!     Search And Replace (SAR) all instances of character CO in STR
!     with the character CN.
!
!     INPUT:
!        CO: a character                     [A1]
!        CN: replacement character           [A1]
!       STR: character variable              [A*]
!
!     OUTPUT:
!       STR: modified character variable     [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 13-NOV-1997
!     REVISIONS HISTORY:
!       27-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from SAR to CHR_SAR
!
!     COMMENTS:
!
      CHARACTER(LEN=1) cn,co
      CHARACTER(LEN=*) str
      INTEGER ilen,i
      INTENT(in) :: cn,co
      INTENT(inout) :: str

      ilen=LEN(str)
      DO i=1,ilen
        IF (str(i:i).EQ.co) str(i:i)=cn
      ENDDO

      RETURN
    END SUBROUTINE CHR_SAR

      SUBROUTINE CHR_X2L(i,j,str)
!     Exchange character str(i:i) with str(j:j).
!
!     INPUT:
!         i: position in str                      [I4]
!         j: position in str                      [I4]
!       str: a character variable                 [A*]
!
!     OUTPUT:
!        str: str(i:i) interchanged with str(j:j) [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-FEB-1996
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       27-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from SAC to CHR_X2L
!
!     COMMENTS:
!
      INTEGER i,j,ilen
      CHARACTER(LEN=*) str
      CHARACTER(LEN=1) swap
      INTENT(in) :: i,j
      INTENT(inout) :: str

      ilen=LEN(str)
      IF ((i.GT.0).AND.(j.GT.0)) THEN
        IF ((i.NE.j).AND.(i.LE.ilen).AND.(j.LE.ilen)) THEN
!         Swap characters
          swap = str(j:j)
          str(j:j) = str(i:i)
          str(i:i) = swap
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE CHR_X2L

      FUNCTION LENTRIM(str)
!
!     Compute the length of a character variable minus trailing
!     blanks and any unprintable characters (ASCII code < 33
!     or ASCII > 254).
!
!     INPUT:
!       str: a string of characters [A*]
!
!     OUTPUT:
!        LENTRIM: length of str minus trailing blanks and [I4]
!                 unprintable characters.
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-MAY-1992
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      LOGICAL again
      INTEGER asc,nchr,itmp
      CHARACTER(LEN=*) str
      INTENT(in) :: str
      INTEGER(I4K) LENTRIM

      nchr = LEN(str)
      again = (nchr.GT.0)
      DO WHILE (again)
        asc = ICHAR(str(nchr:nchr))
        IF (asc.LT.33) THEN
!          non-printable characters & blank
           itmp = nchr - 1
           nchr = itmp
        ELSEIF (asc.GT.254) THEN
!          null character
           itmp = nchr - 1
           nchr = itmp
        ELSE
          again = .false.
        ENDIF
        IF (nchr.LT.1) THEN
!         Nothing left to check
          again = .false.
        ENDIF
      ENDDO
      LENTRIM = nchr

      RETURN
    END FUNCTION LENTRIM

      SUBROUTINE SQUEEZE(str)
!     Replace sequences of multiple blank spaces in STR with a single
!     blank space, e.g., STR= 'Rob   Stewart' will be returned as
!     'Rob Stewart').
!
!     INPUT:
!       str: character variable  [A*]
!
!     OUTPUT:
!       str: modified character variable  [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-MAY-1992
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str
      CHARACTER(LEN=1) ch
      INTEGER iasc,icnt,ilen,i
      LOGICAL space
      INTENT(inout) :: str

      ilen = LEN(str)
      IF (ilen.GT.1) THEN
        icnt = 1
        space = .false.
        DO I=1,ilen
          iasc = ICHAR(str(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
            ch = str(i:i)
            str(icnt:icnt) = ch
            icnt = icnt + 1
            SPACE = .true.
          ELSEIF (space) THEN
            ch = str(i:i)
            str(icnt:icnt) = ch
            icnt = icnt + 1
            space = .false.
          ENDIF
        ENDDO
        DO I=icnt,ilen
          str(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
    END SUBROUTINE SQUEEZE

      SUBROUTINE STR_COPY(loc,str1, str2)
!     Copy the character variable STR1 into STR2 beginning
!     at position LOC in STR2 [i.e., str2(loc:) = str1].  Some
!     characters in STR1 may not be copied to STR2 if STR2 is
!     not declared large enough to hold all of (STR1+STR2).
!
!     INPUT:
!         loc: location to begin copying str1 into str2       [I4]
!        str1: a character variable                           [A*]
!
!     OUTPUT:
!        str2: modified character variable . Characters from  [A*]
!              str2(1:loc-1) are not modified.
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 24-DEC-1993
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       28-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Also, renamed
!          routine from CpyStr to STR_COPY.
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str2,str1
      INTEGER loc,ilt1,i1,i,j,ic
      INTENT(in) :: loc,str1
      INTENT(inout) :: str2

      i1 = LEN(str2)
      IF (i1.GT.0) THEN
        ilt1 = LEN(str1)
        IF (ilt1.GT.0) THEN
          ic = MAX0(loc,1)
          i = 1
          j = ic
          DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
            str2(j:j) = str1(i:i)
            i = i + 1
            j = ic + (i-1)
          ENDDO
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE STR_COPY

      SUBROUTINE STR_FILL(ch, str)
!     Overwrite all characters in STR with the character in CH.
!
!     INPUT:
!        ch: a character          [A1]
!       str: character variable   [A*]
!
!     OUTPUT:
!        str: set str = ch        [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Also, renamed
!          routine from FillStr to STR_FILL.
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str
      CHARACTER(LEN=1) ch
      INTEGER ilen,i
      INTENT(in) :: ch
      INTENT(inout) :: str

      ilen = LEN(str)
      DO i=1,ilen
        str(i:i) = ch
      ENDDO

      RETURN
    END SUBROUTINE STR_FILL

      SUBROUTINE STR_SHFT(nspc,str)
!     Shift all of the characters in STR to the left (nspc < 0)
!     or right (nspc > 0) by NSPC spaces.  Blanks are added
!     to STR to pad the newly exposed elements of STR.
!
!     INPUT:
!       nspc: number times string should be shifted to the  [I4]
!             left or right.
!        str: a character variable                          [A*]
!
!     OUTPUT:
!        str: modified character variable                   [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       28-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Also, renamed
!          routine from ShtStr to STR_SHFT.

!     COMMENTS: This routine is similar to the Fortran 90
!               instrinsic routines ADJUSTL and ADJUSTR,
!               except the STR_SHFT routine gives more
!               control over the process.
!
      CHARACTER(LEN=*) str
      INTEGER nspc
      INTEGER num,ilen,i,j
      CHARACTER(LEN=1) ch
      INTENT(in) :: nspc
      INTENT(inout) :: str

      ilen = LEN(str)
      num = ABS(nspc)
      IF (nspc.LT.0) THEN
!       Shift characters to the left
        DO i=1,num
          DO j=2,ilen
!           NOTE: ch is used as to directly avoid self-referencing
!           a character variable.
            ch = str(j:j)
            str(j-1:j-1) = ch
          ENDDO
          str(ilen:ilen) = ' '
        ENDDO
      ELSE
!       Shift characters to the right
        DO i=1,num
          DO j=ilen-1,1,-1
            ch = str(j:j)
            str(j+1:j+1) = ch
          ENDDO
          str(1:1) = ' '
        ENDDO
      ENDIF

      RETURN
    END SUBROUTINE STR_SHFT

      SUBROUTINE STR_TRIM(ch,str)
!     Trim (delete) all characters in STR beyond the first occurance
!     of CH in STR, i.e., look for the 1st occurance of CH in STR and
!     then replace all characters, including the search character, by
!     blank spaces.
!
!     INPUT:
!        ch: character variable  [A1]
!       str: character variable  [A*]
!
!     OUTPUT:
!       str: modified character variable  [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-MAR-1996
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       22-MAR-1998 RD STEWART
!          BUG FIX: With some compilers, an array bounds error is
!          produced when IX=0.
!       27-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from TrimStr to STR_TRIM
!
!     COMMENTS: The search for character CH in STR *is*
!               case sensitive.
!
      CHARACTER(LEN=*) str
      CHARACTER(LEN=1) ch
      INTEGER i,ix,ilen
      INTENT(in) :: ch
      INTENT(inout) :: str

      ix=INDEX(str,ch)
      IF (ix.GT.0) THEN
        ilen=LEN(str)
        DO i=ix,ILEN
          str(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
    END SUBROUTINE STR_TRIM


      SUBROUTINE STRIP(str)
!     Remove ("strip") all blank spaces and unprintable characters
!     (ASCII code < 33 or ASCII > 254) from STR.
!
!     INPUT:
!       str: character variable            [A*]
!
!     OUTPUT:
!       str: modified string of characters [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-MAY-1992
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str
      INTEGER ilen,iasc,icnt,i
      CHARACTER(LEN=1) ch
      INTENT(inout) :: str

      ilen = LEN(str)
      IF (ilen.GT.1) THEN
        icnt = 1
        DO i=1,ilen
          iasc = ICHAR(str(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
!           Keep character
            ch = str(i:i)
            str(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
!       Fill remainder of str with blanks
        DO i=icnt,ilen
          str(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
    END SUBROUTINE STRIP

      SUBROUTINE TRANLC(str)
!     Translate all uppercase characters in STR to their equivalent
!     lowercase counterpart. Non-alphabetic characters are not
!     altered.
!
!     INPUT:
!       str: character variable            [A*]
!
!     OUTPUT:
!       str: modified string of characters [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-MAY-1992
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str
      INTEGER iasc,i,ilen
      INTENT(inout) :: str

      ilen = LEN(str)
      DO I=1,ilen
        iasc = ICHAR(str(i:i))
        IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
          str(i:i) = CHAR(iasc+32)
        ENDIF
      ENDDO

      RETURN
    END SUBROUTINE TRANLC

      SUBROUTINE TRANUC(str)
!     Translate all lowercase characters in STR to their equivalent
!     uppercase counterpart. Non-alphabetic characters are not
!     altered.
!
!     INPUT:
!       str: character variable            [A*]
!
!     OUTPUT:
!       str: modified string of characters [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-MAY-1992
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90
!
!     COMMENTS:
!
      CHARACTER(LEN=*) str
      INTEGER iasc,i,ilen
      INTENT(inout) :: str

      ilen = LEN(str)
      DO i=1,ilen
        iasc = ICHAR(str(i:i))
        IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          str(i:i) = CHAR(iasc-32)
        ENDIF
      ENDDO

      RETURN
    END SUBROUTINE TRANUC


      SUBROUTINE SP_CATSTR(str1,div, str2)
!     Add (or concatenate) STR1 to STR2 with string
!     DIV between them, i.e., STR2=(STR2+DIV+STR1).
!
!     INPUT:
!        STR1: a character variable          [A*]
!         DIV: a character variable          [A*]
!        STR2: a character variable          [A*]
!
!     OUTPUT:
!        STR2: modified character variable   [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-DEC-1995
!     REVISIONS HISTORY:
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,div,str2
      INTEGER(I4K) ilen,JLEN
      INTENT(in) :: str1,div
      INTENT(inout) :: str2

      ILEN=LENTRIM(str2)+1
      CALL STR_FILL(' ',str2(ilen:))
      JLEN=MAX(1,LENTRIM(str1))
      IF (ILEN.GE.1) THEN
        CALL STR_COPY(ilen,div // str1(:jlen), str2)
      ELSE
        CALL STR_COPY(ilen,str1(:jlen), str2)
      ENDIF

      RETURN
    END SUBROUTINE SP_CATSTR

      FUNCTION CHR_CNT(str1,str2)
!     Return as CHR_CNT the number of times a character in STR1
!     matches at least one character in STR2.  The CHR_CNT routine
!     is not a symmetric operation in the sense that CHR_CNT(str1,str2)
!     does not necessarily return the same number as CHR_CNT(str2,str1).
!
!     INPUT:
!        STR1: a character variable           [A*]
!        STR2: a character variable           [A*]
!
!     OUTPUT:
!        CHR_CNT: an integer variable         [I4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 17-JAN-1995
!     REVISIONS HISTORY:
!       06-MAR-1996
!         Modified algorithm so that trailing blanks in
!         STR1 are ignored.
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       02-JAN-1999 RD STEWART
!         Updated comments to Fortran 90, changed to routine to
!         an integer function.
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2
      INTEGER cnt
      INTEGER L1,L2,i,j
      LOGICAL again
      INTENT(in) :: str1,str2

      L1 = MAX(1,LENTRIM(str1))
      L2 = LEN(str2)
      cnt = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            cnt = cnt+1
            again = .false.
          ELSEIF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO
      CHR_CNT=CNT

      RETURN
    END FUNCTION CHR_CNT

      FUNCTION CHR_REP(str1,str2)
!     Return as CHR_REP the number of times that a character STR1
!     also appears in STR2.  The CHR_REP routine is a symmetric
!     operation in the sense that CHR_REP(str1,str2) returns
!     the same number as CHR_REP(str2,str1).
!
!     INPUT:
!        STR1: a character variable           [A*]
!        STR2: a character variable           [A*]
!
!     OUTPUT:
!        CHR_CNT: an integer variable         [I4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 17-JAN-1995
!     REVISIONS HISTORY:
!     CREATION DATE: 17-JAN-1995
!     REVISIONS HISTORY:
!       06-MAR-1996
!         Modified algorithm so that trailing blanks in
!         STR1 are ignored.
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       02-JAN-1999 RD STEWART
!         Updated comments to Fortran 90, changed to routine to
!         an integer function.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2
      INTEGER cnt
      INTEGER L1,L2,i,j
      LOGICAL again
      INTENT(in) :: str1,str2

      L1 = MAX(1,LENTRIM(str1))
      L2 = LEN(str2)
      cnt = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            cnt = cnt+1
          ENDIF
          IF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO
      CHR_REP=CNT

      RETURN
    END FUNCTION CHR_REP

      SUBROUTINE FMT_D4RN(str, fmat,iflg)
!     Decode the short-hand, real number format specification
!     string, STR, into a standard Fortran format specifier
!     and return it as FMAT [e.g., '(F80.3)'].
!
!     INPUT:
!        STR: format specifier                          [A*]
!              Fn = float with n decimal places
!              En = scientific notation with n+1
!                   significant digits.
!               I = integer-like format
!              Gn = F or E format, as appropriate
!
!     OUTPUT:
!        FMAT: character variable                      [A*]
!        IFLG: status flag                             [I4]
!               1: mimick integer format requested.
!               0: task successful.
!              -1: FMAT declared too small.
!              -2: invalid STR
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JAN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: This routine is called by the PUTI4N, PUTR4N and
!               PUTR8N routines, and it is mainly useful as a
!               convenience for the programmer.
!
      CHARACTER(LEN=*) str,fmat
      INTEGER iflg
      INTEGER i,j,k,l,il,i4n
      INTENT(in) :: str
      INTENT(out) :: fmat,iflg

      IF (LEN(fmat).LT.6) THEN
        !insufficient FMAT workspace
        IFLG=-1
        RETURN
      ELSEIF (LENTRIM(str).LT.1) THEN
        !invalid specifier
        IFLG=-2
        RETURN
      ENDIF

!     Erase existing information in FMAT.
      CALL STR_FILL(' ',fmat)

!     Decode short-hand format specifier.
      I=MAX(INDEX(str,'E'),INDEX(str,'e'))
      J=MAX(INDEX(str,'F'),INDEX(str,'f'))
      K=MAX(INDEX(str,'I'),INDEX(str,'i'))
      L=MAX(INDEX(str,'G'),INDEX(str,'g'))
      IL=MAX(I,J,K,L)+1
      CALL STR_GETI4N(str(IL:), i4n,iflg)
      IF (iflg.NE.0) THEN
        !default to 7 decimal places
        i4n=7
      ELSEIF (i4n.GT.70) THEN
        !too many decimal places. default to
        !to maximum allowed.
        i4n=70
      ENDIF

!     Make format specifier.
      IF (i.GT.0) THEN
!       Use exponential format
        FMAT='(1PE80.  )'
        CALL WRTI4N(i4n,'(I2)', fmat(8:9), iflg)
        iflg=0
      ELSEIF (j.GT.0) THEN
!       Use floating point format
        FMAT='(F80.  )'
        CALL WRTI4N(i4n,'(I2)', fmat(6:7), iflg)
        iflg=0
      ELSEIF (k.GT.0) THEN
!       Use floating point format to mimick integer format.
        FMAT='(F80.0)'
        iflg=1
      ELSEIF (l.GT.0) THEN
!       Use G-format
        FMAT='(G80.  )'
        CALL WRTI4N(i4n,'(I2)', fmat(6:7), iflg)
        iflg=0
      ELSE
!       Default to exponential format
        FMAT='(1PE80.7)'
        iflg=0
      ENDIF
      CALL STRIP(fmat)

      RETURN
    END SUBROUTINE FMT_D4RN

      SUBROUTINE PARSE(iopt,lst,str, nxt, tok,iflg)
!     Find the first sequence of characters in STR such that
!     all of the characters in this sequence are in LST.  The
!     search begins at position NXT in STR [i.e., at STR(nxt:)].
!     If the search is successful (iflg=0), the identified sequence
!     of characters is returned as TOK.  On return, NXT is incremented
!     to the location where another sequence of characters might
!     begin.  Trailing blanks in LST are ignored unless (1)
!     LENTRIM(lst)=0 or (2) the blank space is followed by
!     a non-blank character.
!
!     INPUT:
!         IOPT: task options                                  [I4]
!               0: case insensitive search (default)
!               1: case sensitive search
!          LST: a character variable                          [A*]
!          STR: a character variable                          [A*]
!          NXT: position in STR                               [I4]
!
!     OUTPUT:
!          TOK: a subset of the characters in STR (a "token") [A*]
!          NXT: incremented to position of next possible      [I4]
!               group of characters in STR.
!         IFLG: status flag  [I4]
!                0: sequence found
!               -1: sequence not found
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-MAR-1996
!     REVISIONS HISTORY:
!
!       09-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       19-JUN-1998 RD STEWART
!          Bug fix: changed NXT=NXT+1 to NXT=I2+1
!       31-DEC-1998 RD STEWART
!          Updated comments to Fortran 90.  Other misc.
!          changes to make it compatible with updated
!          FLIB90 routines.
!
!     COMMENTS: PARSE and TOKEN are different methods of
!               obtaining a subset of characters from
!               a character variable.
!
      USE NUMERIC_KINDS
      USE OTHER_USEFUL
      CHARACTER(LEN=*) lst,str,tok
      INTEGER iopt,nxt,iflg
      INTEGER i1,i2,ilen,jlen
      LOGICAL again
      INTENT(in) :: iopt,lst,str
      INTENT(inout) :: nxt
      INTENT(out) :: tok,iflg

!     Erase any characters in TOK
      CALL STR_FILL(' ',tok)

!     Compute length of STR and trap bad values of NXT.
      ilen=LEN(str)
      IF (nxt.GT.ilen) THEN
        iflg=-1   !sequence not found
        nxt=1     !Reset to 1st position in string.
        RETURN
      ELSEIF (nxt.LE.0) THEN
!       Trap invalid value of NXT
        nxt=1
      ENDIF
      jlen=MAX(1,LENTRIM(lst))

!     Get rid of any delimiters at the beginning of STR(NXT:)
      DO WHILE ((.NOT.STR_IMAT(iopt,str(nxt:nxt),lst(:jlen)).AND.(nxt.LT.ilen)))
        nxt=nxt+1
      ENDDO

!     Process token.
      IF (STR_IMAT(iopt,str(nxt:nxt),lst)) THEN
        i1=nxt    !1st letter in token.
!       Find last letter in token.
        again=(nxt.LT.ilen)
        DO WHILE (again)
          nxt=nxt+1
          again=(STR_IMAT(iopt,str(nxt:nxt),lst(:jlen)).AND.(nxt.LT.ilen))
        ENDDO
        i2=nxt-1
        CALL STR_COPY(1,str(i1:i2),tok)
        iflg=0
      ELSE
        iflg=-1   !sequence not found
        nxt=1     !Reset nxt
      ENDIF

      RETURN
    END SUBROUTINE PARSE


      SUBROUTINE SP_PUTSTR(str1,loc,nj, str2,iflg)
!     Copy STR1 into STR2 at position LOC using left (NJ=-1),
!     center (NJ=0), or right (NJ=1) justification.  Trailing
!     blanks in STR1 are ignored (not copied) into STR2, except
!     for the special case when STR1=' '.
!
!     INPUT:
!        STR1: character string to copy into STR2  [A*]
!         LOC: location in STR to put I4N          [I4]
!          NJ: type of justification to use        [I4]
!              -1: left justify STR1 at LOC
!               0: center STR1 about LOC
!               1: right justify STR1 at LOC
!
!     OUTPUT:
!        STR2: modified character variable         [A*]
!        IFLG: status flag                         [I4]
!               0: task successful.
!              -1: task failed
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!       02-JAN-1999 RD STEWART
!         Changed comments to Fortran 90 style.
!
!     COMMENTS:
!
      CHARACTER(LEN=*) STR1,STR2
      INTEGER LOC,NJ,IFLG
      INTEGER il,j1,ilen
      INTENT(in) :: str1,loc,nj
      INTENT(inout) :: str2
      INTENT(out) :: iflg

      ilen=MAX(1,LENTRIM(str1))
      IF (ilen.GT.1) THEN
        IF (nj.EQ.(-1)) THEN
!         left justify at LOC
          IL=LOC
        ELSEIF (nj.EQ.0) THEN
!         Center about LOC.
          IL=(LOC-ILEN/2)
        ELSE
!         Right justify at LOC.
          IL=LOC-ILEN+1
        ENDIF
      ELSE
        !value of NJ doesn't matter
        IL=LOC
      ENDIF

      j1=1
      IF (IL.LT.1) THEN
        iflg=-1
        j1=2-il
        il=1
        IF (j1.LE.0) RETURN  !no part of STR1 copied to STR2
      ELSEIF (il+ilen.GT.LEN(str2)) THEN
        iflg=-1
      ELSE
        iflg=0
      ENDIF
      CALL STR_COPY(il,str1(j1:ilen), str2)

      RETURN
    END SUBROUTINE SP_PUTSTR

      SUBROUTINE STR_ADD(str1, str2)
!     Add (or concatenate) STR1 to the end of STR2.
!     Blank spaces at the end of STR2 are ignored.
!     However, blank spaces at the beginning of
!     STR1 are significant (i.e., copied into STR2).
!
!     INPUT:
!        STR1: a character variable          [A*]
!        STR2: a character variable          [A*]
!
!     OUTPUT:
!        STR2: modified character variable   [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 27-SEP-1994
!     REVISIONS HISTORY:
!       22-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Routine
!          renamed to STR_ADD.
!
!     COMMENTS:
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2
      INTEGER ilen
      INTENT(in) :: str1
      INTENT(inout) :: str2

      ilen=LENTRIM(str2)+1
      CALL str_copy(ilen,str1,str2)

      RETURN
    END SUBROUTINE STR_ADD

      SUBROUTINE STR_DAS(str1, str2,iflg)
!     Find the first case insensitive match for STR1 in STR2
!     and then delete it.
!
!     INPUT:
!        STR1: a character variable          [A*]
!        STR2: a character variable          [A*]
!
!     OUTPUT:
!        STR2: modified character variable   [A*]
!        IFLG: status flag                   [I4]
!              0: STR1 located and deleted
!             -1: STR1 not found in STR2
!             -2: bad trouble in STR_FAM
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-DEC-1998
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2
      INTEGER j1,j2,nspc,iflg
      INTENT(in) :: str1
      INTENT(inout) :: str2
      INTENT(out) :: iflg

      CALL STR_FAM(str1,str2, j1,j2,iflg)
      IF (iflg.NE.0) THEN
        !STR1 not found in STR2
        iflg=-1
      ELSE
        nspc=j2-j1 + 1
        IF (nspc < 1) THEN
          iflg=-2
        ELSE
          !Delete STR1 from STR2
          CALL STR_SHFT(-nspc, str2(j1:))
          iflg=0
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE STR_DAS

      SUBROUTINE STR_DXS(str1, str2,iflg)
!     Find the first exact match for STR1 in STR2
!     and then delete it.
!
!     INPUT:
!        STR1: a character variable          [A*]
!        STR2: a character variable          [A*]
!
!     OUTPUT:
!        STR2: modified character variable   [A*]
!        IFLG: status flag                   [I4]
!              0: STR1 located and deleted
!             -1: STR1 not found in STR2
!             -2: bad trouble in STR_FXM
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-DEC-1998
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2
      INTEGER j1,j2,nspc,iflg
      INTENT(in) :: str1
      INTENT(inout) :: str2
      INTENT(out) :: iflg

      CALL STR_FXM(str1,str2, j1,j2,iflg)
      IF (iflg.NE.0) THEN
        !STR1 not found in STR2
        iflg=-1
      ELSE
        nspc=j2-j1 + 1
        IF (nspc < 1) THEN
          iflg=-2
        ELSE
          !Delete STR1 from STR2
          CALL STR_SHFT(-nspc, str2(j1:))
          iflg=0
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE STR_DXS

      SUBROUTINE STR_FAM(str1,str2, j1,j2,iflg)
!     Attempt to locate a case insensitive match for STR1 in STR2.  If
!     the search is successfully (iflg=0), the location of STR1 in
!     STR2 is returned as J1,J2, i.e., STR1 occurs at STR2(J1:J2).
!     Trailing blanks in STR1 are ignored.
!
!     INPUT:
!       STR1: a chracter variable                              [A*]
!       STR2: a chracter variable                              [A*]
!
!     OUTPUT:
!         J1: location of the 1st letter of STR1 in STR2       [I4]
!         J2: location of the last letter of STR1 in STR2      [I4]
!         IFLG: status flag                                    [I4]
!               -1: sub-string STR1 not found in STR2
!                0: found STR1 in STR2
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-NOV-1997
!     REVISIONS HISTORY:
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS      !interface statements
      USE OTHER_USEFUL
      CHARACTER(LEN=*) STR1,STR2
      INTEGER J1,J2,IFLG
      INTEGER i,j,ilen,jlen
      LOGICAL mat,again,more,found
      INTENT(in) :: str1,str2
      INTENT(out) :: j1,j2,iflg

      !Initialize output parameters.
      iflg=0
      j1=1
      j2=1

!     Check for bad input.
      ilen=LENTRIM(str1)  !length minus trailing blanks
      jlen=LEN(str2)
      IF (jlen.EQ.0) THEN
        !STR2 is a null string
        iflg=-1
        RETURN
       ELSEIF (ilen.EQ.0) THEN
         ilen=LEN(str1)
         IF (ilen.EQ.0) THEN
           !STR1 is a null string
           iflg=-1
           RETURN
         ENDIF
       ENDIF

!     Initial condition
      j1=1
      j2=1
      iflg=-1   !no sub-string found

!     Begin search for sub-string
      j=0
      again=(jlen.GE.ilen)   !STR2 must be longer than STR1 to find a match.
      DO WHILE (again)
        j=j+1
        found=no
        mat=CHR_IMAT(str1(1:1),str2(j:j))
        IF (mat) THEN
!         ith letter in STR1 matches the jth letter
!         in STR2. Set position of 1st letter and
!         continue until all of STR1 is found, a character
!         in STR1 doesn't match one in STR2, or the
!         end of STR1 or STR2 is encountered.
          j1=j
          i=1
          found=(ilen.EQ.1)
          more=((.not.found).AND.((j.LT.jlen).AND.(i.LT.ilen)))
          DO WHILE (more)
            i=i+1
            j=j+1
            mat=CHR_IMAT(str1(i:i),str2(j:j))
            found=(mat.AND.(i.EQ.ilen))
            more=((mat).AND.((i.LT.ilen).AND.(j.LT.jlen)))
          ENDDO
          IF (found) THEN
            j2=j
            iflg=0
          ELSE
!           Reset index to STR1
            j=j1
          ENDIF
        ENDIF
        again=((.not.found).AND.(j.LT.jlen))
      ENDDO
      IF (iflg.NE.0) THEN
        j1=1
        j2=1
      ENDIF

      RETURN
    END SUBROUTINE STR_FAM

      SUBROUTINE STR_FXM(str1,str2, j1,j2,iflg)
!     Attempt to locate an exact match for STR1 in STR2.  If the
!     search is successfully (iflg=0), the location of STR1 in
!     STR2 is returned as J1,J2, i.e., STR1 occurs at STR2(J1:J2).
!     Trailing blanks in STR1 are ignored.
!
!     INPUT:
!       STR1: a chracter variable                              [A*]
!       STR2: a chracter variable                              [A*]
!
!     OUTPUT:
!         J1: location of the 1st letter of STR1 in STR2       [I4]
!         J2: location of the last letter of STR1 in STR2      [I4]
!         IFLG: status flag                                    [I4]
!               -1: sub-string STR1 not found in STR2
!                0: found STR1 in STR2
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-NOV-1997
!     REVISIONS HISTORY:
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.
!
!     COMMENTS:
!
      USE NUMERIC_KINDS      !interface statements
      USE OTHER_USEFUL
      CHARACTER(LEN=*) STR1,STR2
      INTEGER J1,J2,IFLG
      INTEGER i,j,ilen,jlen
      LOGICAL mat,again,more,found
      INTENT(in) :: str1,str2
      INTENT(out) :: j1,j2,iflg

      !Initialize output parameters.
      iflg=0
      j1=1
      j2=1

!     Check for bad input.
      ilen=LENTRIM(str1)  !length minus trailing blanks
      jlen=LEN(str2)
      IF (jlen.EQ.0) THEN
        !STR2 is a null string
        iflg=-1
        RETURN
       ELSEIF (ilen.EQ.0) THEN
         ilen=LEN(str1)
         IF (ilen.EQ.0) THEN
           !STR1 is a null string
           iflg=-1
           RETURN
         ENDIF
       ENDIF

!     Initial condition
      iflg=-1   !no sub-string found

!     Begin search for sub-string
      j=0
      again=(jlen.GE.ilen)   !STR2 must be longer than STR1 to find a match.
      DO WHILE (again)
        j=j+1
        found=no
        mat=(str1(1:1).EQ.str2(j:j))
        IF (mat) THEN
!         ith letter in STR1 matches the jth letter
!         in STR2. Set position of 1st letter and
!         continue until all of STR1 is found, a character
!         in STR1 doesn't match one in STR2, or the
!         end of STR1 or STR2 is encountered.
          j1=j
          i=1
          found=(ilen.EQ.1)
          more=((.not.found).AND.((j.LT.jlen).AND.(i.LT.ilen)))
          DO WHILE (more)
            i=i+1
            j=j+1
            mat=(str1(i:i).EQ.str2(j:j))
            found=(mat.AND.(i.EQ.ilen))
            more=((mat).AND.((i.LT.ilen).AND.(j.LT.jlen)))
          ENDDO
          IF (found) THEN
            j2=j
            iflg=0
          ELSE
!           Reset index to STR1
            j=j1
          ENDIF
        ENDIF
        again=((.not.found).AND.(j.LT.jlen))
      ENDDO
      IF (iflg.NE.0) THEN
        j1=1
        j2=1
      ENDIF

      RETURN
    END SUBROUTINE STR_FXM

      FUNCTION STR_IMAT(iopt,ch,str)
!     Return STR_IMAT as TRUE if CH is a case insensitive
!     (iopt=0) or case sensitive (iopt=1) match for one of
!     the characters listed in STR.  Trailing blanks in STR
!     are ignored unless (1) LENTRIM(str)=0 or (2) the blank
!     is followed by one or more non-blank characters.
!
!     INPUT:
!       iopt: task options                        [I4]
!             1: case sensitive comparison
!             0: case insensitive (default
!         ch: a character                         [A1]
!        str: character variable                  [A*]
!
!     OUTPUT:
!       STR_IMAT: logical variable indicating     [L4]
!                 whether CH is a case sensitive
!                 or insensitive match to one of
!                 the characters in STR.
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 31-DEC-1998
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE NUMERIC_KINDS
USE OTHER_USEFUL
USE MATH_USEFUL_DP

  INTEGER iopt,ilen,i
  CHARACTER(LEN=1) CH
  CHARACTER(LEN=*) STR
  LOGICAL found,again,STR_IMAT
  INTENT(in) :: iopt,ch,str

    ilen=MAX(1,LENTRIM(str))
    found=no
    again=(LEN(ch)>0)
    i=0
    DO WHILE (again)
      i=i+1
      found=((ch.EQ.str(i:i)).OR.((iopt.NE.1).AND.CHR_IMAT(ch,str(i:i))))
      again=((i.LT.ilen).AND.(.not.found))
    ENDDO
    STR_IMAT=found

RETURN
END FUNCTION STR_IMAT

      SUBROUTINE STR_INSERT(loc,str1,str2)
!     Insert STR1 into STR2 beginning at position LOC in STR2.
!     Existing characters in STR2 will be shifted to the right.
!     Trailing blanks in STR1 are significant.
!
!     INPUT:
!         LOC: a location in STR2             [I4]
!        STR1: character variable             [A*]
!        STR2: character variable             [A*]
!
!     OUTPUT:
!        STR2: modified character variable    [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-DEC-1993
!     REVISIONS HISTORY:
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Routine
!          renamed from INSSTR to STR_INSERT.
!
!     COMMENTS:  The length of STR2 on return is always >=
!                the length of STR2 on entry.
!
      CHARACTER(LEN=*) str1,str2
      INTEGER loc,LI,ic,LT
      INTENT(in) :: loc,str1
      INTENT(out) :: str2

      lt = LEN(str2)
      IF (lt.GT.0) THEN
        li = LEN(str1)
        ic = MAX0(loc,1)
        CALL STR_SHFT(li,str2(ic:lt))
        CALL STR_COPY(ic,str1,str2)
      ENDIF

      RETURN
    END SUBROUTINE STR_INSERT

      SUBROUTINE STR_RXS(str1,str2, str3,iflg)
!     Find the first exact match of STR1 in STR3
!     and then replace it with str2.
!
!     INPUT:
!        STR1: a character variable          [A*]
!        STR2: a character variable          [A*]
!        STR3: a character variable          [A*]
!
!     OUTPUT:
!        STR3: modified character variable   [A*]
!        IFLG: status flag                   [I4]
!              0: STR1 located and deleted
!             -1: STR1 not found in STR2
!             -2: bad trouble in STR_FXM
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-DEC-1998
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2,str3
      INTEGER j1,j2,nspc,iflg
      INTENT(in) :: str1,str2
      INTENT(inout) :: str3
      INTENT(out) :: iflg

      CALL STR_FXM(str1,str3, j1,j2,iflg)
      IF (iflg.NE.0) THEN
        !STR1 not found in STR3
        iflg=-1
      ELSE
        nspc=j2-j1 + 1
        IF (nspc < 1) THEN
          iflg=-2
        ELSE
          !Delete STR1 from STR3
          CALL STR_SHFT(-nspc, str3(j1:))
          !Insert STR2 into STR3 at position j1
          CALL STR_INSERT(j1,str2, str3)
          iflg=0
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE STR_RXS

      SUBROUTINE STR_RAS(str1,str2, str3,iflg)
!     Find the first case insensitive match of STR1 in STR3
!     and then replace it with str2.
!
!     INPUT:
!        STR1: a character variable          [A*]
!        STR2: a character variable          [A*]
!        STR3: a character variable          [A*]
!
!     OUTPUT:
!        STR3: modified character variable   [A*]
!        IFLG: status flag                   [I4]
!              0: STR1 located and deleted
!             -1: STR1 not found in STR2
!             -2: bad trouble in STR_FAM
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-DEC-1998
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str1,str2,str3
      INTEGER j1,j2,nspc,iflg
      INTENT(in) :: str1,str2
      INTENT(inout) :: str3
      INTENT(out) :: iflg

      CALL STR_FAM(str1,str3, j1,j2,iflg)
      IF (iflg.NE.0) THEN
        !STR1 not found in STR3
        iflg=-1
      ELSE
        nspc=j2-j1 + 1
        IF (nspc < 1) THEN
          iflg=-2
        ELSE
          !Delete STR1 from STR3
          CALL STR_SHFT(-nspc, str3(j1:))
          !Insert STR2 into STR3 at position j1
          CALL STR_INSERT(j1,str2, str3)
          iflg=0
        ENDIF
      ENDIF

      RETURN
    END SUBROUTINE STR_RAS

      FUNCTION STR_SCAN(STR,LST)
!     Return as STR_SCAN the position of the first character
!     in STR that matches one of the characters in LST.  Trailing
!     blanks in LST are ignored except when (1) LENTRIM(lst) = 0
!     and (2) when the blank is followed by one or more non-blank
!     characters (e.g., LST='< >').
!
!     INPUT:
!        LST: a character variable          [A*]
!        STR: a character variable          [A*]
!
!     OUTPUT:
!        STR_SCAN: an integer varialbe      [I4]
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-MAR-1996
!     REVISIONS HISTORY:
!       09-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       02-JAN-1999 RD STEWART
!          Updated to comments to Fortran 90.
!
!     COMMENTS: This routine is more or less equivalent to the
!               Fortran 90 intrinsic function SCAN.
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) lst,str
      INTEGER ilen,jlen,i,j,cnt,jmax
      INTEGER(I4K) STR_SCAN
      INTENT(in) :: str,lst

      ilen=MAX(1,LENTRIM(lst))
      jlen=LEN(str)
      jmax=jlen
      cnt=0
      DO i=1,ilen
!       position of lst(i:i) in str
        j=INDEX(str,lst(i:i))
        IF (j.GT.0) THEN
!         Found a match
          jmax=MIN(jmax,j)
          cnt=cnt+1
        ENDIF
      ENDDO
      IF (cnt.LT.1) THEN
        jmax=0
      ENDIF
      STR_SCAN=jmax

      RETURN
    END FUNCTION STR_SCAN

      SUBROUTINE STR_XLET(str)
!      Interchange characters in str so that, for example,
!      when str='trebor' on entry it becomes str='robert'
!      on return.  Trailing blanks in STR are ignored.
!
!     INPUT:
!        str: a string of characters representing a str [A*]
!
!     OUTPUT:
!        str: modified string of characters [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-FEB-1996
!     REVISIONS HISTORY:
!
!       09-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Routine
!          renamed from XCWORD to STR_XLET
!
!     COMMENTS: A 2nd call to STR_XLET will restore the original
!               sequence of letters.
!
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER ilen,I,J,imid
      INTENT(inout) :: str

      ilen= LENTRIM(str)
      IF (ilen.GT.0) THEN
        imid = ilen/2
        DO i=1,imid
          j = ilen-i + 1
          CALL CHR_X2L(i,j,str)
        ENDDO
      ENDIF

      RETURN
    END SUBROUTINE STR_XLET

      SUBROUTINE STR_XWORD(str)
!      Search for words (groups of characters separated by one or
!      more spaces or punctuation marks) in STR and then interchange
!      the characters in each word so that, for example,
!      'trebor trawets si a dog!' becomes 'robert stewart is a god!'
!      Note: A second call to STR_XWORD will restore the original
!      sequence of letters.
!
!     INPUT:
!        str: character variable          [A*]
!
!     OUTPUT:
!        str: modified character variable [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-FEB-1996
!     REVISIONS HISTORY:
!
!       09-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       10-JUL-1998 RD STEWART
!          Renamed routine from PLSOW to XCSOW
!       31-DEC-1998 RD STEWART
!          Updated comments to Fortran 90 style and renamed
!          routine to STR_XWORD.
!
!     COMMENTS: Not a particularly useful routine, but writing it
!               was *fun* and it can also be used a very
!               simplistic method of encrypting text.
!
      CHARACTER(LEN=*) str
      CHARACTER(LEN=80) word,pg
      INTEGER iopt,col,iflg,i1,i2,ilen
      INTENT(inout) :: str

      iopt=0  !case insensitive
      iflg=0
      col = 1
      DO WHILE (iflg.EQ.0)
        word = ' '
        CALL TOKEN(iopt,'":;,. !?''''',str, col, word,iflg)
        IF (iflg.EQ.0) THEN
          pg = ' '
          pg = word
          CALL STR_XLET(pg)
          ilen=LENTRIM(word)
          i1 = INDEX(str,word(:ilen))
          i2 = i1 + ilen-1
          str(i1:i2) = pg(:ilen)
        ENDIF
      ENDDO

      RETURN
    END SUBROUTINE STR_XWORD

      SUBROUTINE TOKEN(iopt,lst,str, nxt, tok,iflg)
!     Find the first sequence of characters in STR delimited
!     by one of the characters in LST.  The search begins at
!     position NXT in STR [i.e., at STR(nxt:)]. If the search
!     is successful (iflg=0), the identified sequence of
!     characters is returned as TOK.  On return, NXT is incremented
!     to the location where another sequence of characters might
!     begin.  Trailing blanks in LST are ignored unless (1)
!     LENTRIM(lst)=0 or (2) the blank space is followed by
!     a non-blank character.
!
!     INPUT:
!         IOPT: task options                                  [I4]
!               0: case insensitive search (default)
!               1: case sensitive search
!          LST: a character variable                          [A*]
!          STR: a character variable                          [A*]
!          NXT: position in STR                               [I4]
!
!     OUTPUT:
!          TOK: a subset of the characters in STR (a "token") [A*]
!          NXT: incremented to position of next possible      [I4]
!               group of characters in STR.
!         IFLG: status flag                                   [I4]
!                0: token found
!               -1: token not found
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 06-MAR-1996
!     REVISIONS HISTORY:
!
!       09-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       19-JUN-1998 RD STEWART
!          Bug fix: changed NXT=NXT+1 to NXT=I2+1
!       31-DEC-1998 RD STEWART
!          Updated comments to Fortran 90 style + misc.
!          changes to make it compatible with FLIB90.
!
!     COMMENTS: PARSE and TOKEN are different methods of
!               obtaining a subset of characters from a
!               another string of characters.
!
      USE NUMERIC_KINDS
      USE OTHER_USEFUL
      CHARACTER(LEN=*) lst,str,tok
      INTEGER iopt,nxt,iflg
      INTEGER i1,i2,ilen,jlen
      LOGICAL again
      INTENT(in) :: iopt,lst,str
      INTENT(inout) :: nxt
      INTENT(out) :: tok,iflg

!     Erase any characters in TOK
      CALL STR_FILL(' ',tok)

!     Compute length of STR and trap bad values of NXT.
      ilen=LEN(str)
      IF (nxt >= ilen) THEN
        iflg=-1   !no token found.
        nxt=1     !Reset to 1st position in string.
        RETURN
      ELSEIF (nxt.LE.0) THEN
!       Trap invalid value of NXT
        nxt=1
      ENDIF
      jlen=MAX(1,LENTRIM(lst))

!     Get rid of any delimiters at the beginning of STR(NXT:)
      DO WHILE (STR_IMAT(iopt,str(nxt:nxt),lst(:jlen)).AND.(nxt.LT.ilen))
        nxt=nxt+1
      ENDDO

!     Process token.
      IF (.not.STR_IMAT(iopt,str(nxt:nxt),lst(:jlen))) THEN
        i1=nxt    !1st letter in token.
!       Find last letter in token.
        again=(nxt<ilen)
        DO WHILE (again)
          nxt=nxt+1
          again=(.not.STR_IMAT(iopt,str(nxt:nxt),lst(:jlen)).AND.(nxt < ilen))
        ENDDO
        i2=nxt-1
        IF(nxt == ilen) i2=nxt
        CALL STR_COPY(1,str(i1:i2),tok)
        iflg=0
      ELSE
        iflg=-1   !Didn't find any tokens in STR
        nxt=1     !Reset nxt
      ENDIF

      RETURN
    END SUBROUTINE TOKEN
    SUBROUTINE CHANNEL(io)
!     Find an unused file handle (or 'channel') that can be
!     used for file i/o.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     INPUT:
!
!     OUTPUT:
!           io: an used device number that can be used    [I4]
!                for file i/o
!
!     AUTHOR: Peter Gans
!             School of Chemistry,
!             University of Leeds,
!             Leeds LS2 9JT
!             England
!
!     CREATION DATE: 11-MAY-1995
!     REVISIONS HISTORY:
!       06-DEC-1995 RD STEWART
!          Added to FLIB package and modified so that search
!          for a free file handle begins at unit=10.
!
!       05-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!
!       29-NOV-1998 RD STEWART
!          Updated to meet F90 language standard.
!
!     COMMENTS:
!
      INTEGER io
      LOGICAL opened
      INTENT(out) :: io

      opened=.TRUE.
      io=9
      DO WHILE (opened)
         io=io+1
         inquire(io,opened=opened)
      ENDDO

      RETURN
    END SUBROUTINE CHANNEL

END MODULE STRPAK
MODULE CONSTANTS
  INTEGER, PARAMETER :: max_pars=200,max_char_tree = 80,&
       & max_char_long = 12000, max_data=12000,max_char=120,max_atm=10
  CHARACTER(len=1), DIMENSION(2), PARAMETER :: Comms=(/'!','#'/)
  CHARACTER(len=max_char), DIMENSION(9), PARAMETER :: Used=(/'ATOM  ','H&
       &ETATM','CONECT','SSBOND','CRYST1','SEQRES','HELIX ','SHEET ','CI&
       &SPEP'/)
  REAL(8), PARAMETER :: Huge=1.0D10,Tiny=1.0D-10
END MODULE CONSTANTS
MODULE Print_Defs
  IMPLICIT none
  INTEGER, SAVE :: kprint=6
END MODULE Print_Defs
MODULE ERROR_List
  USE CONSTANTS, ONLY: max_char,max_data
  USE Print_Defs
  IMPLICIT NONE 
  INTEGER, SAVE :: counter=0, count_f=0, count_w=0
  INTEGER, PARAMETER :: max_err=4*max_char, max_err_long=12000
  INTEGER, PARAMETER :: Max_errors=50
  TYPE Orac_Errors
     INTEGER :: tag
     INTEGER :: count
     CHARACTER(len=max_err_long) :: err_text
     TYPE(Orac_Errors), POINTER :: next
  END TYPE Orac_Errors
  TYPE(Orac_Errors),  POINTER, SAVE :: root, current
CONTAINS
  SUBROUTINE Start()
    IMPLICIT none
    IF(ASSOCIATED(root)) THEN
       CALL Cleanup()
    END IF
    ALLOCATE(current)
    NULLIFY(current%next)
    root=>current
  END SUBROUTINE Start
  SUBROUTINE Add(No,error)
    IMPLICIT none
    CHARACTER(len=*), INTENT(IN) :: error
    INTEGER, INTENT(IN) :: No
    IF(counter == 0) THEN
       CALL Start()
    END IF
    IF(counter > Max_errors) RETURN
    counter=counter+1
    IF(No < 0) THEN
       count_f=count_f+1
       current % count = count_f
    ELSE
       count_w=count_w+1
       current % count = count_w
    END IF
    current % tag = No
    current % err_text = error
    ALLOCATE(current % next)
    current=>current % next
    NULLIFY(current % next)
  END SUBROUTINE Add
  SUBROUTINE Cleanup()
    IMPLICIT none
    TYPE(Orac_Errors), POINTER :: dummy
    
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       dummy=>current % next
       DEALLOCATE(current)
       current=>dummy
    END DO
    NULLIFY(dummy, root)
    count_w=0
  END SUBROUTINE Cleanup
END MODULE ERROR_List
MODULE TYPES
  USE CONSTANTS
  TYPE List
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: g
  END TYPE List
  TYPE KHASH
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: keys
     INTEGER :: n,Method=0
  END TYPE KHASH
  TYPE Param
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: line
  END TYPE Param
  TYPE Docs
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
  END TYPE Docs
  TYPE Commands_D
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
     TYPE(Docs), DIMENSION(:), POINTER :: keys
  END TYPE Commands_D
  TYPE Environment_D
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
     TYPE(commands_D), DIMENSION(:), POINTER :: Comm
  END TYPE Environment_D

  TYPE link_tble
     CHARACTER(len=max_char) :: key
     TYPE(link_tble), DIMENSION(:), POINTER :: n_key
  END TYPE link_tble

  TYPE keys
     CHARACTER(len=max_char) :: label
     CHARACTER(len=max_char), DIMENSION(:), POINTER :: line
     TYPE(keys), DIMENSION(:), POINTER :: next
  END TYPE keys
  TYPE PATCH
     CHARACTER(8) :: Type
     CHARACTER(8) :: New_Res,pres,res
     CHARACTER(8) :: Res_l(2)
     INTEGER :: one,two
  END TYPE PATCH
  TYPE Chain
     INTEGER, DIMENSION(:), ALLOCATABLE :: g
  END TYPE Chain
  TYPE ChainR8
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE ChainR8
  TYPE ChainR8C
     CHARACTER(len=max_char) :: lab
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE ChainR8C
END MODULE TYPES
MODULE Node
  USE CONSTANTS
  PRIVATE
  PUBLIC :: Node_,Node__Delete,Node__Push,Node__Pop&
       &,Node__Size
  TYPE :: NODES
     CHARACTER(len=1), DIMENSION(:), POINTER  :: vector=>NULL()
     TYPE(NODES), POINTER :: next=>NULL()
  END TYPE NODES

  TYPE(NODES), POINTER :: root=>NULL()
  TYPE(NODES), POINTER :: current=>NULL()
  INTEGER, SAVE :: counter=0
  INTERFACE Node__Push
     MODULE PROCEDURE Node_I__Push
     MODULE PROCEDURE Node_I_1__Push
     MODULE PROCEDURE Node_R8__Push
     MODULE PROCEDURE Node_R8_1__Push
     MODULE PROCEDURE Node_CL__Push
     MODULE PROCEDURE Node_CL_1__Push
  END INTERFACE
  INTERFACE Node__Pop
     MODULE PROCEDURE Node_I_1__Pop
     MODULE PROCEDURE Node_I__Pop
     MODULE PROCEDURE Node_R8__Pop
     MODULE PROCEDURE Node_R8_1__Pop
     MODULE PROCEDURE Node_CL__Pop
     MODULE PROCEDURE Node_CL_1__Pop
  END INTERFACE
CONTAINS
  FUNCTION Node_() RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    IF(counter /= 0) THEN
       counter=0
       CALL Node__Delete()
    END IF
    ALLOCATE(root)
    NULLIFY(root%next)
    current=>root
    out=.TRUE.
  END FUNCTION Node_

  FUNCTION Node__Size() RESULT(out)
    INTEGER :: out
    out=counter
  END FUNCTION Node__Size

  SUBROUTINE Node__delete()
    CALL Remove(root)
    NULLIFY(root)
    NULLIFY(current)
  END SUBROUTINE Node__delete
  RECURSIVE SUBROUTINE Remove(old_node)
    TYPE(NODES), POINTER :: old_node
    IF(.NOT. ASSOCIATED(old_node) ) RETURN
    CALL Remove(old_node % next)
    IF(ASSOCIATED(old_node % vector)) DEALLOCATE(old_node % vector)
    DEALLOCATE(old_node)
  END SUBROUTINE Remove

!!$
!!$--- Overloading 
!!$

  SUBROUTINE Node_C__Push(vect)
    IMPLICIT none
    CHARACTER(len=1), DIMENSION(:), INTENT(IN) :: vect

    TYPE(NODES), POINTER :: new_node
    counter=counter+1
    ALLOCATE(new_node)
    ALLOCATE(new_node % vector(SIZE(vect)))
    new_node % vector = vect
    NULLIFY(new_node % next)

    current % next => new_node
    current => current % next

  END SUBROUTINE Node_C__Push

  FUNCTION Node_C__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=1), DIMENSION(:), POINTER :: vect
    INTEGER, SAVE :: count0=0

    IF(count0 == 0) THEN
       current=>root
    END IF
    
    count0=count0+1
    out=ASSOCIATED(current % next )
    IF(.NOT. out) THEN
       count0=0
       RETURN
    END IF
 
    current=>current % next
    IF(ASSOCIATED(vect)) DEALLOCATE(Vect)
    ALLOCATE(Vect(SIZE(current % vector)))
    vect=current % vector

  END FUNCTION Node_C__Pop

!!$
!!$--- Hack to get around bug on the TRANSFER() function of the pgi compiler!
!!$--- version 6.2-5
!!$

  FUNCTION MyTransfer_Str2Char(vect) RESULT(out)
    CHARACTER(len=Max_Char), DIMENSION(:) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: out
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: cvect0
    INTEGER :: n,m,p

    n=Max_char*SIZE(vect)
    IF(ALLOCATED(cvect0)) DEALLOCATE(cvect0)
    ALLOCATE(cvect0(n))
    p=0
    DO n=1,SIZE(vect)
       DO m=1,LEN(vect(n))
          p=p+1
          cvect0(p)=vect(n)(m:m)
       END DO
    END DO
    out=>cvect0
  END FUNCTION MyTransfer_Str2Char

  FUNCTION MyTransfer_Char2Str(cvect) RESULT(out)
    CHARACTER(len=1), DIMENSION(:) :: cvect

    CHARACTER(len=Max_Char), DIMENSION(:), POINTER :: out
    CHARACTER(len=Max_Char), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: vect0
    INTEGER :: n,m,p

    n=SIZE(cvect)/Max_Char
    IF(ALLOCATED(vect0)) DEALLOCATE(vect0)
    ALLOCATE(vect0(n))

    p=0
    DO n=1,SIZE(vect0)
       DO m=1,LEN(vect0(n))
          p=p+1
          vect0(n)(m:m)=cvect(p)
       END DO
    END DO
    out=>vect0
  END FUNCTION MyTransfer_Char2Str
!!$
!!$--- Hack to get around bug in the pgi compiler end!
!!$

!!$
!!$--- Overloaded routines
!!$

  SUBROUTINE Node_CL__Push(vect)
    CHARACTER(len=max_Char), DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect

    cvect=>MYTransfer_Str2Char(vect)

    CALL Node_C__Push(cvect)

  END SUBROUTINE Node_CL__Push

  SUBROUTINE Node_CL_1__Push(vect)
    CHARACTER(len=max_Char), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    CHARACTER(len=max_Char) :: vect0(1)

    vect0(1)=vect
    cvect=>MYTransfer_Str2Char(vect0)

    CALL Node_C__Push(cvect)

  END SUBROUTINE Node_CL_1__Push

  FUNCTION Node_CL__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: vect
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)

    IF(.NOT. out) RETURN

    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    aux=>MyTransfer_Char2Str(cvect)
    ALLOCATE(vect(SIZE(aux)))
    vect=aux
    DEALLOCATE(cvect)
  END FUNCTION Node_CL__Pop
  FUNCTION Node_CL_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=max_char), INTENT(OUT) :: vect
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN

    aux=>MyTransfer_Char2Str(cvect)
    vect=aux(1)
    DEALLOCATE(cvect)
  END FUNCTION Node_CL_1__Pop

  SUBROUTINE Node_R8__Push(vect)
    REAL(8), DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect(length))
    cvect=TRANSFER(vect,cvect)
    CALL Node_C__Push(cvect)
    DEALLOCATE(cvect)
  END SUBROUTINE Node_R8__Push

  SUBROUTINE Node_R8_1__Push(vect)
    REAL(8), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect(length))
    cvect=TRANSFER(vect,cvect)
    CALL Node_C__Push(cvect)
    DEALLOCATE(cvect)
  END SUBROUTINE Node_R8_1__Push
  FUNCTION Node_R8__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    REAL(8), DIMENSION(:), POINTER :: vect
    REAL(8), DIMENSION(:), POINTER  :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    ALLOCATE(vect(SIZE(TRANSFER(cvect,aux))))
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_R8__Pop
  FUNCTION Node_R8_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    REAL(8), INTENT(OUT) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_R8_1__Pop
  SUBROUTINE Node_I__Push(vect)
    INTEGER, DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect(length))
    cvect=TRANSFER(vect,cvect)
    CALL Node_C__Push(cvect)
    DEALLOCATE(cvect)
  END SUBROUTINE Node_I__Push
  SUBROUTINE Node_I_1__Push(vect)
    INTEGER, INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect(length))
    cvect=TRANSFER(vect,cvect)
    CALL Node_C__Push(cvect)
    DEALLOCATE(cvect)
  END SUBROUTINE Node_I_1__Push
  FUNCTION Node_I__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    INTEGER, DIMENSION(:), POINTER :: vect
    INTEGER, DIMENSION(:), POINTER  :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    IF(ASSOCIATED(vect)) DEALLOCATE(vect)

    ALLOCATE(vect(SIZE(TRANSFER(cvect,aux))))
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_I__Pop
  FUNCTION Node_I_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    INTEGER, INTENT(OUT) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_I_1__Pop
END MODULE Node
MODULE Errors
  USE Node
  USE ERROR_List
  USE TYPES, ONLY: list
  IMPLICIT NONE 
  PRIVATE
  PUBLIC abort_now, abort_later, warning, add, print_errors,&
       & error_args,error_unr ,error_file,error_other&
       &,Setup_Errors,errmsg_f,errmsg_w, Print_warnings

  CHARACTER(LEN=45), SAVE :: err_arg_no='.. Arguments to command must be at least '
  TYPE fatal_err
     CHARACTER(74), DIMENSION (3) :: top= &
          &(/ '**************************************************************************'&
          &,'*                                                                        *'&
          &,'*-- The following FATAL errors were found:         ----------------------*'/)
     CHARACTER(74), DIMENSION (:), ALLOCATABLE :: body
     CHARACTER(74), DIMENSION(1) :: intrabodies = &
          (/ '*                                                                        *'/)
     CHARACTER(74), DIMENSION(2) :: bottom= &
          (/ '*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE fatal_err

  TYPE warning_err
     CHARACTER(74), DIMENSION (3) :: top= &
          (/ '**************************************************************************' &
          &,'*                                                                        *'&
          &,'*-- The following WARNING errors were found:         --------------------*'/)
     CHARACTER(74), DIMENSION(1) :: intrabodies = &
          (/ '*                                                                        *'/)
     CHARACTER(74), DIMENSION (:), ALLOCATABLE :: body
     CHARACTER(74), DIMENSION(2) :: bottom= &
          (/ '*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE warning_err

  TYPE(List), SAVE :: error_args,error_unr,error_file,error_other
  CHARACTER(len=max_data) :: errmsg_f,errmsg_w
CONTAINS

  SUBROUTINE abort_now(msg)
    !*******************************************************************
    !                                                                  *
    !  Abort NOW.                                                      *
    !  msg1: Main Error message                                        *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    CHARACTER(*)  :: msg

    TYPE(fatal_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()

    n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string

!!$-- Make the error box

    nn=0
    nlast=0
    IF(.NOT. Node_()) STOP
    DO WHILE(nlast /= n)
       nn_old=nn          
       nlast=MIN(nn+nl,n)
       nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
       i_nn=nn-nn_old-1
       vector0=' '
       vector0(1:1)='*'
       vector0(74:74)='*'
       vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
       CALL Node__Push(vector0)
    END DO
        
    m=Node__Size()
    
    ALLOCATE(msg2 % body(m))
    i=0
    DO WHILE(Node__Pop(vector))
       i=i+1
       msg2%body(i)=vector(1)
    END DO
    
    WRITE(kprint,'(2x,a)') msg2%top
    WRITE(kprint,'(2x,a)') msg2%body(1:m)
    WRITE(kprint,'(2x,a)') msg2%bottom
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
    
    STOP

  END SUBROUTINE abort_now
  SUBROUTINE abort_later(msg)
!!$!*******************************************************************
!!$!                                                                  *
!!$!  Abort.                                                          *
!!$!  msg1: Main Error message                                        *
!!$!                                                                  *
!!$!*******************************************************************
    IMPLICIT NONE
    CHARACTER(*), OPTIONAL  :: msg

    TYPE(fatal_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()
    
    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string

!!$-- Make the error box

       nn=0
       nlast=0
       IF(.NOT. Node_()) STOP
       DO WHILE(nlast /= n)
          nn_old=nn          
          nlast=MIN(nn+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-nn_old-1
          vector0=' '
          vector0(1:1)='*'
          vector0(74:74)='*'
          vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
          CALL Node__Push(vector0)
       END DO
       
       m=Node__Size()
       ALLOCATE(msg2 % body(m))
       i=0
       DO WHILE(Node__Pop(vector))
          i=i+1
          msg2%body(i)=vector(1)
       END DO

       IF(counter == 0) WRITE(kprint,'(2x,a)') msg2%top
       WRITE(kprint,'(2x,a)') msg2%intrabodies
       WRITE(kprint,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(kprint,'(2x,a)') msg2%bottom
    END IF
    counter=counter+1
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
  END SUBROUTINE abort_later
  SUBROUTINE Warning(msg, reset)
!!$!*******************************************************************
!!$!                                                                  *
!!$!  Warn.                                                           *
!!$!  msg1: Main Error message                                        *
!!$!                                                                  *
!!$!*******************************************************************
    IMPLICIT NONE
    CHARACTER(*), OPTIONAL  :: msg
    INTEGER, OPTIONAL  :: reset

    TYPE(warning_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()


    IF(PRESENT(reset)) THEN
       counter=0    
       RETURN
    END IF
    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string
       
!!$-- Make the error box


       nn=0
       nlast=0
       IF(.NOT. Node_()) STOP
       DO WHILE(nlast /= n)
          nn_old=nn          
          nlast=MIN(nn+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-nn_old-1
          vector0=' '
          vector0(1:1)='*'
          vector0(74:74)='*'
          vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
          CALL Node__Push(vector0)
       END DO
       
       m=Node__Size()
       ALLOCATE(msg2 % body(m))
       i=0
       DO WHILE(Node__Pop(vector))
          i=i+1
          msg2%body(i)=vector(1)
       END DO

       IF(counter == 0) WRITE(kprint,'(2x,a)') msg2%top
       WRITE(kprint,'(2x,a)') msg2%intrabodies
       WRITE(kprint,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(kprint,'(2x,a)') msg2%bottom
    END IF
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
    counter=counter+1
  END SUBROUTINE Warning
  SUBROUTINE Print_Warnings()
    LOGICAL :: warnings_found
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,No_Warn,count

    error=' '
    No_Warn=0
    warnings_found=.FALSE.
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       No=current % tag
       count=current % count
       WRITE(dummy,'(i3,'') --'')') count
       Error=TRIM(dummy)//' '//current % err_text
       IF(No > 0) THEN
          current % tag = 0
          No_Warn=No_Warn+1
          IF(No_Warn == 1) CALL Warning(' ',0)
          CALL Warning(Error)
          warnings_found=.TRUE.
       END IF
       current=>current % next
    END DO
    IF(warnings_found) THEN
       CALL Warning(); count_w=0
    END IF
  END SUBROUTINE Print_Warnings
  SUBROUTINE Print_Errors()
    LOGICAL :: stop_run
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,count

    
    error=' '
    IF(.NOT. ASSOCIATED(root)) RETURN
    stop_run=.FALSE.
    current=>root    
    DO WHILE(ASSOCIATED(current % next))
       No=current % tag
       count=current % count
       WRITE(dummy,'(i3,'') --'')') count
       Error=TRIM(dummy)//' '//current % err_text
       IF(No < 0) THEN
          stop_run=.TRUE.
          CALL Abort_Later(Error)
       END IF
       current=>current % next
    END DO
    IF(stop_run) THEN
       CALL Abort_Later
       STOP
    END IF
  END SUBROUTINE Print_Errors
  SUBROUTINE Setup_Errors
    ALLOCATE(error_args % g (4),error_unr % g (4), error_file % g (1)&
         &, error_other % g (1)) 

    error_args % g (1)= 'OPEN keyword not found '
    error_args % g (2)= 'Number of arguments must be at least '
    error_args % g (3)= 'Number of arguments must not exceed '
    error_args % g (4)= 'Number of arguments must be only '
    error_unr % g (1)='Unrecognized command  ---> '
    error_unr % g (2)='Unrecognized subcommand -> '
    error_unr % g (3)='Unrecognized keyword ----> '
    error_unr % g (4)='UNSUPPORTED  COMMAND ----> '
    error_file % g (1)='...or missing &END'
    error_other % g (1)=' file not found '

  END SUBROUTINE Setup_Errors
END MODULE Errors
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
MODULE Strings

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 18:29:49 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 13 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Errors, ONLY: Add_Errors=> Add, errmsg_f, Print_Errors
  USE CONSTANTS, ONLY: max_pars, max_char, Comms
  USE TYPES
  USE STRPAK
  INTERFACE Myread
     MODULE PROCEDURE MyRead_I4
     MODULE PROCEDURE MyRead_R8
  END INTERFACE
  INTERFACE Myputnum
     MODULE PROCEDURE MyPutNum_i4
     MODULE PROCEDURE MyPutNum_r8
  END INTERFACE
CONTAINS
  FUNCTION MyPutnum_i4(n) RESULT(out)
    INTEGER  :: n
    CHARACTER(len=max_char) :: out
    CHARACTER(len=max_char) :: str1
    INTEGER :: iflag

    CALL STR_Fill(' ',str1)
    CALL PUTI4N(n,1,-1,'I',str1,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert integer to character '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    out=TRIM(str1)
  END FUNCTION MyPutnum_i4
  FUNCTION MyPutnum_R8(r) RESULT(out)
    REAL(8) :: r
    CHARACTER(len=max_char) :: out
    CHARACTER(len=max_char) :: str1
    INTEGER :: iflag

    CALL PUTR8N(r,1,-1,'R5',str1,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert real to character '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    out=TRIM(str1)
  END FUNCTION MyPutnum_R8
  SUBROUTINE Myread_I4(str1,out)
    CHARACTER(len=*) :: str1
    INTEGER :: out
    
    INTEGER :: iflag
    
    CALL STR_Geti4n(str1,out,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert to integer: '''//TRIM(str1)//''' '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  END SUBROUTINE Myread_I4
  SUBROUTINE Myread_R8(str1, out)
    CHARACTER(len=*) :: str1
    REAL(8) :: out
    
    INTEGER :: iflag
    
    CALL STR_Getr8n(str1,out,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert to double: '''//TRIM(str1)//''' '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  END SUBROUTINE Myread_R8

  FUNCTION My_FAM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fam
    INTEGER :: j1,j2,iflag

    My_FAM=.FALSE.
    CALL STR_FAM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fam=.TRUE.
  END FUNCTION My_FAM
  FUNCTION My_FXM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fxm
    INTEGER :: j1,j2,iflag

    My_FXM=.FALSE.
    CALL STR_FXM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fxm=.TRUE.
  END FUNCTION My_FXM
  FUNCTION ST_Concat(i,strngs)
    IMPLICIT NONE 
    CHARACTER(len=max_char) :: ST_Concat
    INTEGER :: i
    CHARACTER(len=max_pars), DIMENSION(:) :: strngs

    INTEGER :: n,m,pad
    CHARACTER(len=max_char) :: line,aux 

    CALL STR_FILL(' ',line)
    
    pad=1
    m=0
    DO n=i,SIZE(strngs)
       m=m+LEN_TRIM(strngs(n))+pad
    END DO
    IF(m > LEN(line)) THEN
       WRITE(aux,'(i3)') LEN(line)
       errmsg_f='Cannot contatenate strings. Final string length larger&
            & than hard bound '//TRIM(aux)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO n=i,SIZE(strngs)
       line=TRIM(line)//' '//TRIM(strngs(n))
    END DO
    ST_Concat=line
  END FUNCTION ST_Concat
  FUNCTION ST_Clean_Line(line)
    IMPLICIT NONE 
    LOGICAL :: ST_Clean_Line
    CHARACTER(len=max_char) :: line
    INTEGER :: i

    ST_Clean_Line=.FALSE.
    DO i=1,SIZE(Comms)
       CALL STR_TRIM(Comms(i),line)
    END DO
    IF(LEN_TRIM(line) == 0) ST_Clean_Line=.TRUE.
  END FUNCTION ST_Clean_Line
END MODULE Strings
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
MODULE ReadStore
  USE Node
  USE Constants
  USE STRPAK
  USE Strings, ONLY: My_Fam
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors,errmsg_f
  IMPLICIT NONE 
  PRIVATE
  PUBLIC ReadStore_, RS__String, ReadStore__Delete
  CHARACTER(len=max_char), ALLOCATABLE, SAVE :: RS__string(:)
CONTAINS
  FUNCTION ReadStore_(line) RESULT(out)
    LOGICAL :: out
    CHARACTER(len=*) :: line
    LOGICAL :: ok
    CHARACTER(len=max_char) :: linea
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: line1=>NULL()
    INTEGER :: n,io, count_out,count_a,iopt

    IF(ALLOCATED(RS__string)) DEALLOCATE(RS__String)
    out=.TRUE.
    CALL Channel(io)
    INQUIRE(file=line,EXIST=ok)
    IF(ok) THEN
       OPEN(unit=io,file=line,form='FORMATTED',status='OLD')
    ELSE
       out=.FALSE.
       errmsg_f='Trying to read file '''//TRIM(line)//''',  which does not exist'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    IF(.NOT. Node_()) STOP
    DO 
       READ(unit=io,fmt='(a80)',IOSTAT=iopt) linea
       IF(iopt /= 0) EXIT
       ok=.FALSE.
       DO n=1,SIZE(Used)
          IF(My_Fam(TRIM(Used(n)),linea)) ok=.TRUE.
       END DO
       IF(.NOT. ok) CYCLE
       CALL Node__Push(linea)
    END DO
    count_out=Node__Size()
    ALLOCATE(RS__string(count_out))
    
    count_a=0
    DO WHILE(Node__Pop(line1))
       count_A=count_A+1
       RS__string(count_a)=ADJUSTL(TRIM(line1(1)))
    END DO
    CALL Node__Delete()
    CLOSE(io)
  END FUNCTION ReadStore_
  SUBROUTINE ReadStore__Delete
    IF(ALLOCATED(RS__string)) DEALLOCATE(RS__String)
  END SUBROUTINE ReadStore__Delete
END MODULE ReadStore
