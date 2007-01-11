
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

MODULE MATH_USEFUL_SP
   !Some useful mathematical constants in single precision
   USE NUMERIC_KINDS

   REAL(KIND=R4K) zero,one,two,three,four,five,six,seven,eight,nine,ten
   PARAMETER (zero=0.0E+00,one=1.0E+00,two=2.0E+00,three=3.0E+00)
   PARAMETER (four=4.0E+00,five=5.0E+00,six=6.0E+00,seven=7.0E+00)
   PARAMETER (eight=8.0E+00,nine=9.0E+00,ten=10.0E+00)

   REAL(KIND=R4K) half,third
   PARAMETER (HALF=0.5E+00,THIRD=one/three)

   REAL(KIND=R4K) PI,TWOPI,FOURPI,HALFPI
   PARAMETER (PI=3.14159265358979323E+00)
   PARAMETER (HALFPI=HALF*PI,TWOPI=TWO*PI,FOURPI=FOUR*PI)

   !Angle conversion
   REAL(KIND=R4K) RAD2DEG,DEG2RAD
   PARAMETER (RAD2DEG=180.0E+00/PI)  !Convert radians to degress
   PARAMETER (DEG2RAD=PI/180.0E+00)  !Convert degrees to radians

   !Time conversion (based on sidereal year)
   REAL(KIND=R4K) M2H,H2M,D2H,H2D
   REAL(KIND=R4K) S2H,H2S,S2D,D2S,S2Y,Y2S
   REAL(KIND=R4K) M2Y,Y2M,H2Y,Y2H,D2Y,Y2D,D2M,M2D
   PARAMETER (H2S=3600.0E+00,S2H=one/H2S)    !second <-> hour
   PARAMETER (D2S=24.0E+00*H2S,S2D=one/D2S)  !second <-> day
   PARAMETER (Y2S=3.155815E+07,S2Y=one/Y2S)  !second <-> year
   PARAMETER (D2M=1440.0E+00,M2D=one/D2M)    !minute <-> day
   PARAMETER (Y2M=Y2S/60.0E+00,M2Y=one/Y2M)  !minute <-> year
   PARAMETER (H2M=60.0E+00,M2H=one/H2M)      !hour <-> minute
   PARAMETER (D2H=24.0E+00,H2D=one/D2H)      !hour <-> day
   PARAMETER (Y2H=Y2S*S2H,H2Y=one/Y2H)       !hour <-> year

END MODULE MATH_USEFUL_SP

MODULE OTHER_USEFUL
  USE NUMERIC_KINDS
  LOGICAL(L4K) YES,NO
  PARAMETER (YES=.TRUE.,NO=.FALSE.)

  INTEGER(I4K), PARAMETER :: DIAGNOSTIC = 0
  INTEGER(I4K), PARAMETER :: WARNING = 1
  INTEGER(I4K), PARAMETER :: FATAL = -1
  INTEGER(I4K), PARAMETER :: KillProg = -1

END MODULE OTHER_USEFUL

MODULE FLIB_DEVICE
  !COMPUTER, DEVICE, or OPERATING SYSTEM DEPENDENT INFORMATION
  USE NUMERIC_KINDS

  INTEGER(I4K), PARAMETER :: U4S=6 !device identificaion number for screen

END MODULE FLIB_DEVICE


MODULE STRPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  !FLIB Interface
  INTEGER(I4K)  LENTRIM,CHR_CNT,CHR_REP,STR_SCAN
  EXTERNAL LENTRIM,CHR_CNT,CHR_REP,STR_SCAN
  LOGICAL  CHR_ILC,CHR_IUC,CHR_ILET,CHR_ILON,CHR_INUM,CHR_IPRN
  LOGICAL  CHR_IMAT,STR_IMAT
  EXTERNAL CHR_ILC,CHR_IUC,CHR_ILET,CHR_ILON,CHR_INUM,CHR_IPRN
  EXTERNAL CHR_IMAT,STR_IMAT

  INTERFACE SP_CATNUM
    SUBROUTINE CATI4N(i4n,fmat,div, str)
      USE NUMERIC_KINDS
      INTEGER(I4K) i4n
      CHARACTER(LEN=*) fmat,div,str
    END SUBROUTINE CATI4N
    SUBROUTINE CATR4N(r4n,fmat,div, str)
      USE NUMERIC_KINDS
      REAL(R4K) r4n
      CHARACTER(LEN=*) fmat,div,str
    END SUBROUTINE CATR4N
    SUBROUTINE CATR8N(r8n,fmat,div, str)
      USE NUMERIC_KINDS
      REAL(R8K) r8n
      CHARACTER(LEN=*) fmat,div,str
    END SUBROUTINE CATR8N
  END INTERFACE

  INTERFACE SP_GETNUM
    SUBROUTINE STR_GETI4N(str, i4n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) i4n,iflg
      CHARACTER(LEN=*) str
    END SUBROUTINE STR_GETI4N
    SUBROUTINE STR_GETR4N(str, r4n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) iflg
      REAL(R4K) r4n
      CHARACTER(LEN=*) str
    END SUBROUTINE STR_GETR4N
    SUBROUTINE STR_GETR8N(str, r8n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) iflg
      REAL(R8K) r8n
      CHARACTER(LEN=*) str
    END SUBROUTINE STR_GETR8N
  END INTERFACE

  INTERFACE SP_PUTNUM
    SUBROUTINE PUTI4N(i4n,loc,nj,fmat, str,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) i4n,loc,nj,iflg
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE PUTI4N
    SUBROUTINE PUTR4N(r4n,loc,nj,fmat, str,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) loc,nj,iflg
      REAL(R4K) r4n
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE PUTR4N
    SUBROUTINE PUTR8N(r8n,loc,nj,fmat, str,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) loc,nj,iflg
      REAL(R8K) r8n
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE PUTR8N
  END INTERFACE

  INTERFACE SP_WRTNUM
    SUBROUTINE WRTI4N(i4n,fmat, str, iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) i4n,iflg
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE WRTI4N
    SUBROUTINE WRTR4N(r4n,fmat, str, iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) iflg
      REAL(R4K) r4n
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE WRTR4N
    SUBROUTINE WRTR8N(r8n,fmat, str, iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) iflg
      REAL(R8K) r8n
      CHARACTER(LEN=*) str,fmat
    END SUBROUTINE WRTR8N
  END INTERFACE

END MODULE STRPAK


MODULE RANPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  INTEGER(I4K) RP_INT
  EXTERNAL     RP_INT

  !Generic real number in interval (a,b)
  INTERFACE RP_UNI
    FUNCTION RP_U8N(ir)
      USE NUMERIC_KINDS
      REAL(R8K) RP_U8N
      INTEGER(I4K) IR
    END FUNCTION RP_U8N
    FUNCTION RP_R8N(ir,a,b)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) A,B,RP_R8N
    END FUNCTION RP_R8N
    FUNCTION RP_R4N(ir,a,b)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) A,B,RP_R4N
    END FUNCTION RP_R4N
  END INTERFACE

  !Points from a plane
  INTERFACE RP_PLANE
    SUBROUTINE RP_R4PLANE(ir,w,h, x,y)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) w,h,x,y
    END SUBROUTINE RP_R4PLANE
    SUBROUTINE RP_R8PLANE(ir,w,h, x,y)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) w,h,x,y
    END SUBROUTINE RP_R8PLANE
  END INTERFACE

  !Isotropic direction vector
  INTERFACE RP_ISO
    SUBROUTINE RP_R4ISO(ir, u,v,w)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) u,v,w
    END SUBROUTINE RP_R4ISO
    SUBROUTINE RP_R8ISO(ir, u,v,w)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) u,v,w
    END SUBROUTINE RP_R8ISO
  END INTERFACE

  !Points from a disk
  INTERFACE RP_DISK
    SUBROUTINE RP_R4DISK(ir,r, px,py)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) r,px,py
    END SUBROUTINE RP_R4DISK
    SUBROUTINE RP_R8DISK(ir,r, px,py)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) r,px,py
    END SUBROUTINE RP_R8DISK
  END INTERFACE

  !Points from a ring
  INTERFACE RP_RING
    SUBROUTINE RP_R4RING(ir,r1,r2, px,py)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) r1,r2,px,py
    END SUBROUTINE RP_R4RING
    SUBROUTINE RP_R8RING(ir,r1,r2, px,py)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) r1,r2,px,py
    END SUBROUTINE RP_R8RING
  END INTERFACE

  !Points from a sphere
  INTERFACE RP_SPHERE
    SUBROUTINE RP_R4SPH(ir,r, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) r,px,py,pz
    END SUBROUTINE RP_R4SPH
    SUBROUTINE RP_R8SPH(ir,r, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) r,px,py,pz
    END SUBROUTINE RP_R8SPH
  END INTERFACE

  !Points from a spherical shell
  INTERFACE RP_SHELL
    SUBROUTINE RP_R4SHL(ir,R1,R2, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) r1,r2,px,py,pz
    END SUBROUTINE RP_R4SHL
    SUBROUTINE RP_R8SHL(ir,R1,R2, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) r1,r2,px,py,pz
    END SUBROUTINE RP_R8SHL
  END INTERFACE

  !Points from a cylinder
  INTERFACE RP_CYLINDER
    SUBROUTINE RP_R4CYL(ir,R,H, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) r,h,px,py,pz
    END SUBROUTINE RP_R4CYL
    SUBROUTINE RP_R8CYL(ir,R,H, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) r,h,px,py,pz
    END SUBROUTINE RP_R8CYL
  END INTERFACE

  !Points from a box
  INTERFACE RP_BOX
    SUBROUTINE RP_R4BOX(ir,wx,wy,wz, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) wx,wy,wz,px,py,pz
    END SUBROUTINE RP_R4BOX
    SUBROUTINE RP_R8BOX(ir,wx,wy,wz, px,py,pz)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) wx,wy,wz,px,py,pz
    END SUBROUTINE RP_R8BOX
  END INTERFACE

  !Points from a box
  INTERFACE RP_ICC
    SUBROUTINE RP_R4ICC(ir,mu, u,v,w)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R4K) mu,u,v,w
    END SUBROUTINE RP_R4ICC
    SUBROUTINE RP_R8ICC(ir,mu, u,v,w)
      USE NUMERIC_KINDS
      INTEGER(I4K) IR
      REAL(R8K) mu,u,v,w
    END SUBROUTINE RP_R8ICC
  END INTERFACE

END MODULE RANPAK

MODULE RANPAK_WORKSPACE
   USE NUMERIC_KINDS

   INTEGER(KIND=I4K) NRG,IDEF
   PARAMETER (NRG=100)         !Max. # of generators
   PARAMETER (IDEF=666)        !default random number seed
   INTEGER(KIND=I4K) ISD(NRG)  !Seed for nth generator

   REAL(R8K) a,m,Minv
   PARAMETER (a=16807.0D+00,m=2147483647.0D+00)
   PARAMETER (Minv=1.0D+00/m)

END MODULE RANPAK_WORKSPACE


MODULE FIOPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  !Points from a plane
  INTERFACE VI_DATA
    SUBROUTINE VI_I4NDAT(n, i4n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) n,i4n,iflg
    END SUBROUTINE VI_I4NDAT
    SUBROUTINE VI_R4NDAT(n, r4n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) n,iflg
      REAL(R4K) R4N
    END SUBROUTINE VI_R4NDAT
    SUBROUTINE VI_R8NDAT(n, r8n,iflg)
      USE NUMERIC_KINDS
      INTEGER(I4K) n,iflg
      REAL(R8K) R8N
    END SUBROUTINE VI_R8NDAT
  END INTERFACE

END MODULE FIOPAK

MODULE FIOPAK_WORKSPACE
   USE NUMERIC_KINDS

   INTEGER(KIND=I4K) MXFN      !max. length of filename
   PARAMETER (MXFN=350)
   INTEGER(KIND=I4K) MXLL      !max. line length
   PARAMETER (MXLL=350)

   CHARACTER(LEN=MXFN)    FNI,FNO   !filename in/out
   CHARACTER(LEN=MXFN)    BUF       !temporary workspace
   CHARACTER(LEN=MXFN)    FIOBUF    !temporary workspace
   CHARACTER(LEN=MXFN-90) CWD       !current working directory
   CHARACTER(LEN=1)       CWD_dl    !drive letter
   CHARACTER(LEN=90)      FN_nop    !filename no path
   CHARACTER(LEN=75)      FN_base   !filename base
   CHARACTER(LEN=15)      FN_ext    !filename extension
   LOGICAL                UNIX      !unix-style directory (Y/N)

   !VARIABLES USED TO RETRIEVE DATA FROM ASCII FILES
   !USING VI ROUTINES
   !  DON: data block on/off
   !  CON: comment block on/off
   ! NINC: number of data elements added when the data-input
   !       stack size is increased.
   ! NDES: number of data elements on stack
   ! NDMX: dimensioned size of VI_DAT array
   LOGICAL DON,CON
   INTEGER(KIND=I4K) NINC
   INTEGER(KIND=I4K) :: NDMX=0
   INTEGER(KIND=I4K) :: NDES=0
   PARAMETER (Ninc=5000)
   REAL(R8K), DIMENSION(:), ALLOCATABLE :: VI_DAT  !data storage
   REAL(R8K), DIMENSION(:), ALLOCATABLE :: VI_TMP  !temp. data storage

END MODULE FIOPAK_WORKSPACE


MODULE TIMPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  LOGICAL TP_IsLeapYear       !Define interface to leap year function
  EXTERNAL TP_IsLeapYear

  !INTERFACE TO GET TIME-OF-DAY ROUTINES
  INTERFACE TP_GETTOD
    SUBROUTINE TP_TOD1(tod)
      USE NUMERIC_KINDS
      REAL(R8K) TOD
    END SUBROUTINE
    SUBROUTINE TP_TOD2(tod)
      USE NUMERIC_KINDS
      REAL(R4K) TOD
    END SUBROUTINE
    SUBROUTINE TP_TOD3(hour,minute,sec,msec)
      USE NUMERIC_KINDS
      INTEGER(I4K) hour,minute,sec,msec
    END SUBROUTINE
    SUBROUTINE TP_TOD4(hour,minute,sec)
      USE NUMERIC_KINDS
      INTEGER(I4K) hour,minute
      REAL(R8K) sec
    END SUBROUTINE
    SUBROUTINE TP_TOD5(hour,minute,sec)
      USE NUMERIC_KINDS
      INTEGER(I4K) hour,minute
      REAL(R4K) sec
    END SUBROUTINE
    SUBROUTINE TP_TOD6(hour,minute,sec)
      USE NUMERIC_KINDS
      INTEGER(I4K) hour,minute,sec
    END SUBROUTINE
  END INTERFACE

END MODULE TIMPAK

MODULE TIMPAK_WorkSpace
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  INTEGER(I4K) :: TP_DATE(3) = -1  !year, month, day
  INTEGER(I4K) :: TP_TIME(4) = -1  !hour, minute, second, millisecond
  CHARACTER(LEN=8)  :: TP_DATSTR =' '     !Date String = YYYYMMDD
  CHARACTER(LEN=10) :: TP_TIMSTR =' '     !Time String = HHMMSS.SSS

END MODULE TIMPAK_WorkSpace


MODULE PRNPAK
  !Global PrnPak Variables and Parameters
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP
  USE FLIB_DEVICE

  INTERFACE PP_PRNSTR
    SUBROUTINE PP_PRNSTR0(str,adv)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      LOGICAL(L4K) adv
    END SUBROUTINE PP_PRNSTR0
    SUBROUTINE PP_PRNSTR1(str)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
    END SUBROUTINE PP_PRNSTR1
    SUBROUTINE PP_PRNSTR2(io,str,adv)
      USE NUMERIC_KINDS
      INTEGER(I4K) io
      CHARACTER(LEN=*) str
      LOGICAL(L4K) adv
    END SUBROUTINE PP_PRNSTR2
    SUBROUTINE PP_PRNSTR3(io,str)
      USE NUMERIC_KINDS
      INTEGER(I4K) io
      CHARACTER(LEN=*) str
    END SUBROUTINE PP_PRNSTR3
  END INTERFACE

END MODULE PRNPAK

MODULE PRNPAK_GLOBAL
  !Global PrnPak Variables and Parameters
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP
  USE FLIB_DEVICE

  !FILE I/O
  LOGICAL :: ECHO=NO               !Run-time info OFF by default
  INTEGER(I4K) :: U4E=0            !No run-time echo by default
  INTEGER(I4K) :: U4O=U4S          !device identificaion number PrnPak output device

  !FORMATTING OPTIONS
  LOGICAL :: TextWrap = YES                !Word wrap on by default
  INTEGER(I4K) :: Margin_Left = 0          !Begin printing in column 1
  INTEGER(I4K) :: Margin_Right = 80        !Stop printing in column 79
  CHARACTER(LEN=25) :: Page_Fmat='(T1,A)'

END MODULE PRNPAK_GLOBAL
MODULE ERRPAK
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  INTERFACE EP_NewErr
    SUBROUTINE EP_NewErr0(NOR,EID,ECT)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) nor
      INTEGER(I4K) eid,ect
    END SUBROUTINE EP_NewErr0
    SUBROUTINE EP_NewErr1(NOR,EID)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) nor
      INTEGER(I4K) eid
    END SUBROUTINE EP_NewErr1
  END INTERFACE

  INTERFACE EP_AddNum
    !COMPLETE ARGUMENT SPECIFICATION
    SUBROUTINE EP_AddR8N0(str,r8n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      REAL(R8K) r8n
    END SUBROUTINE EP_AddR8N0
    SUBROUTINE EP_AddR4N0(str,r4n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      REAL(R4K) r4n
    END SUBROUTINE EP_AddR4N0
    SUBROUTINE EP_AddI4N0(str,i4n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) i4n
    END SUBROUTINE EP_AddI4N0

    !DEFAULT NUMBER FORMAT
    SUBROUTINE EP_AddR8N1(str,r8n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      REAL(R8K) r8n
    END SUBROUTINE EP_AddR8N1
    SUBROUTINE EP_AddR4N1(str,r4n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      REAL(R4K) r4n
    END SUBROUTINE EP_AddR4N1
    SUBROUTINE EP_AddI4N1(str,i4n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER(I4K) i4n
    END SUBROUTINE EP_AddI4N1

  END INTERFACE


  INTERFACE EP_AddArr
    !COMPLETE ARGUMENT SPECIFICATION
    SUBROUTINE EP_AddR8A0(str,xa,n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) n
      REAL(R8K) xa(*)
    END SUBROUTINE EP_AddR8A0
    SUBROUTINE EP_AddR4A0(str,xa,n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) n
      REAL(R4K) xa(*)
    END SUBROUTINE EP_AddR4A0
    SUBROUTINE EP_AddI4A0(str,xa,n,fmat)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str,fmat
      INTEGER(I4K) n
      INTEGER(I4K) xa(*)
    END SUBROUTINE EP_AddI4A0

   !DEFAULT NUMBER FORMAT
    SUBROUTINE EP_AddR8A1(str,xa,n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER(I4K) n
      REAL(R8K) xa(*)
    END SUBROUTINE EP_AddR8A1
    SUBROUTINE EP_AddR4A1(str,xa,n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER(I4K) n
      REAL(R4K) xa(*)
    END SUBROUTINE EP_AddR4A1
    SUBROUTINE EP_AddI4A1(str,xa,n)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) str
      INTEGER(I4K) n
      INTEGER(I4K) xa(*)
    END SUBROUTINE EP_AddI4A1
  END INTERFACE

  INTERFACE EP_EndErr
    SUBROUTINE EP_EndErr0(nor,iopt)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) nor
      INTEGER(I4K) iopt
    END SUBROUTINE EP_EndErr0
    SUBROUTINE EP_EndErr1(nor)
      USE NUMERIC_KINDS
      CHARACTER(LEN=*) nor
    END SUBROUTINE EP_EndErr1
    SUBROUTINE EP_EndErr2(iopt)
      USE NUMERIC_KINDS
      INTEGER(I4K) iopt
    END SUBROUTINE EP_EndErr2
    SUBROUTINE EP_EndErr3
    END SUBROUTINE EP_EndErr3
  END INTERFACE

END MODULE ERRPAK

MODULE ERRPAK_GLOBAL
  !Global ERRPAK Variables and Parameters
  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP

  INTEGER(I4K), PARAMETER :: MXML = 320  !Max. message length
  INTEGER(I4K) :: EP_ect=0               !Current error type
  LOGICAL :: ActiveError = No
  LOGICAL :: LineFeed = No

END MODULE ERRPAK_GLOBAL


MODULE FLIB

  USE NUMERIC_KINDS
  USE OTHER_USEFUL
  USE MATH_USEFUL_DP
  USE RANPAK
  USE STRPAK
  USE FIOPAK
  USE TIMPAK
  USE PRNPAK
  USE ERRPAK

  CHARACTER(LEN=*), PARAMETER :: FLIB_VERSION='2.02R'
  CHARACTER(LEN=*), PARAMETER :: FLIB_DATE='04-AUG-1999'

END MODULE FLIB

SUBROUTINE FLIB_PRNVER(io, iflg)
!  Print FLIB version control information to unit=io.
!
!  INPUT:
!      io: a valid file or device handle                     [I4]
!
!  OUTPUT:
!     IFLG: task completion flag                             [I4]
!           0: task successful
!          -1: task failed: invalid file/device handle
!          -2: task failed while print to unit=io.
!
!      LIBRARY: [FLIB]
!       STATUS: [BETA]
!
!     AUTHOR: Rob Stewart, Ph.D.
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             (509) 375-6851  Voice
!             (509) 375-6936  Fax
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JUN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE FLIB
INTEGER IO,IFLG,JFLG
LOGICAL ok
INTENT(IN)  :: IO
INTENT(OUT) :: IFLG

  IF (io.NE.6) THEN
    INQUIRE(unit=io, OPENED = ok)
    IF (.not.ok) THEN
      !Specified file/device handle is invalid.
      IFLG=-1
      RETURN
    ENDIF
  ENDIF
  WRITE(unit=io,iostat=jflg,fmt='(A)') &
     'PNNL Fortran Library (FLIB) ' // FLIB_VERSION // '  ' // FLIB_DATE
  IF (jflg.NE.0) THEN
    IFLG=-1
    RETURN
  ENDIF

RETURN
END


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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      END

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
      EXTERNAL CATR8N
      INTENT(in) :: i4n,fmat,div
      INTENT(inout) :: str

      r8n=i4n
      CALL CATR8N(r8n,fmat,div, str)

      RETURN
      END

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
      EXTERNAL CATR8N

      r8n=r4n
      CALL CATR8N(r8n,fmat,div, str)

      RETURN
      END

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
      INTEGER LENTRIM
      EXTERNAL STR_FILL,STR_COPY,PUTR8N,LENTRIM
      INTENT(in) :: r8n,fmat,div
      INTENT(inout) :: str

      I1=MAX(1,LENTRIM(str))+1               !ignore trailing comments in STR
      I2=LEN(div)+I1                         !all characters in DIV are significant
      CALL STR_FILL(' ',str(i1:))
      CALL STR_COPY(i1,div, str)             !copy DIV to string
      CALL PUTR8N(r8n,i2,-1,fmat, str,iflg)  !write number to end of STR

      RETURN
      END

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
      INTEGER LENTRIM
      EXTERNAL LENTRIM,STR_COPY
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
      END

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
      INTEGER(I4K) CHR_CNT,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER(I4K) CHR_REP,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER i,j,k,l,il,i4n,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER(I4K) LENTRIM
      LOGICAL again
      LOGICAL STR_IMAT
      EXTERNAL STR_IMAT,LENTRIM
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
      END

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
      END

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
      END

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
      INTEGER I,IP,IL,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER il,j1,ilen,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER LENTRIM
      EXTERNAL LENTRIM
      INTENT(in) :: str1
      INTENT(inout) :: str2

      ilen=LENTRIM(str2)+1
      CALL str_copy(ilen,str1,str2)

      RETURN
      END

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
      END

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
      END

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
      LOGICAL mat,again,more,found,CHR_IMAT
      INTEGER LENTRIM
      EXTERNAL CHR_IMAT,LENTRIM
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
      END

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
      INTEGER LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      END

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
      END

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
      INTEGER LENTRIM
      EXTERNAL LENTRIM
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
      END

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

  INTEGER iopt,ilen,i,LENTRIM
  CHARACTER(LEN=1) CH
  CHARACTER(LEN=*) STR
  LOGICAL found,again,CHR_IMAT,STR_IMAT
  EXTERNAL LENTRIM,CHR_IMAT
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
END

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
      END

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
      END

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
      END

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
      INTEGER(I4K) STR_SCAN,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER iopt,col,iflg,i1,i2,ilen,LENTRIM
      EXTERNAL LENTRIM
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
      END

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
      INTEGER i1,i2,ilen,jlen,LENTRIM
      LOGICAL again,STR_IMAT
      EXTERNAL STR_IMAT,LENTRIM
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
      END

SUBROUTINE RP_GETRNS(ir, iseed)
!  Return the current random number seed for the IRth generator.
!
!   INPUT:
!         IR: index to a random number generator       [I4]
!
!   OUTPUT:
!      ISEED: integer seed for the IRth generator      [I4]
!             generator
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: DEC. 1, 1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ir,iseed
INTENT(in) :: ir
INTENT(out) :: iseed

   IF ((ir.GT.0).AND.(ir.LE.nrg)) THEN
     iseed = isd(ir)
   ENDIF

RETURN
END

SUBROUTINE RP_SETRNS(ir, iseed)
!  Initialize one or more RANPAK random number generators.
!
!  INPUT:
!         IR: generator index                             [I4]
!             ir <= 0, initialize all generators.
!            NRG >= IR  > 0, set seed for IRth generator
!    ISEED: integer seed for generators                   [I4]
!           iseed < 1: use default seed
!           iseed > 1: use seed ISEED
!
!     OUTPUT:
!           N/A
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: DEC. 1, 1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90 standard
!
!     COMMENTS:
!       After initialization, the value of the jth
!       random number seed can be accessed with
!       the RANGET routine.
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ir,iseed,iset,i
INTENT(in) :: ir,iseed

  ISET=MAX(0,ISEED)
  IF (ISET.LT.1) ISET=IDEF  !Use default seed

  IF (ir.LT.1) THEN
!   Initialize all random number generator seeds to ISET
    DO i=1,nrg
      isd(i) = iset
    ENDDO
  ELSEIF (ir.LE.nrg) THEN
    isd(ir) = iset
  ENDIF

RETURN
END

FUNCTION RP_INT(ir,ilo,ihi)
!
!     Return an INTEGER number randomly selected from
!     the interval [ILO,IHI] using the irth geneator.
!
!     INPUT:
!          ILO: lower limit of sampling interval  [I4]
!          IHI: upper limit of sampling interval  [I4]
!           IR: index to a random number          [I4]
!               generator  NOTE: IR should
!               be between 0 and NRG.
!
!     OUTPUT:
!         I4X: random number in [ILO,IHI]        [I4]
!               NOTE: the interval is inclusive.
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!
!       26-DEC-1997 RD STEWART
!          Add trap for bad generator indices.
!
!     COMMENTS:
!       This routine should work correctly on any computer
!       with a maximum integer greater than or equal
!       to 2**31 - 1.  NOTE: The proper function of
!       the routine can verified by checking to see
!       that the random number seed after the generation
!       of 10000 random numbers is 1,043,618,065
!
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ilo,ihi,ir,jr
REAL(R8K) tmp,seed
INTEGER(I4K) RP_INT
INTENT(in) :: ir,ilo,ihi

jr=MAX(1,(MIN(nrg,ir)))  !trap bad generator indices
seed = isd(jr)
tmp = A*seed
seed = tmp - M*DINT(tmp*Minv)
RP_INT = INT((ihi-ilo+1)*seed*Minv) + ilo
isd(jr) = seed

RETURN
END

FUNCTION RP_R8N(ir, xlo,xhi)
!  Return a "double precision" real number randomly selected
!  from the interval (xlo,xhi).
!
!     INPUT:
!         IR: random number generator index             [I4]
!        XLO: lower limit on the sampling interval.     [R8]
!        XHI: upper limit on the sampling interval.     [R8]
!
!     OUTPUT:
!         RP_R8N: random number in (0,1)        [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       26-DEC-1997 RD STEWART
!          Add trap for bad generator indices.
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS:
!       This routine should work correctly on any computer
!       with a maximum integer greater than or equal
!       to 2**31 - 1.  NOTE: The proper function of
!       the routine can verified by checking to see
!       that the random number seed after the generation
!       of 10000 random numbers is 1,043,618,065
!
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ir,jr
REAL(R8K) xlo,xhi
REAL(R8K) tmp,seed
REAL(R8K) RP_R8N
INTENT(in) :: ir,xlo,xhi

jr=MAX(1,(MIN(nrg,ir)))  !trap bad generator indices
seed = isd(jr)
tmp = A*seed
seed = tmp - M*DINT(tmp*Minv)
RP_R8N = (xhi-xlo)*seed*Minv + xlo
isd(jr) = seed                !update seed

RETURN
END

FUNCTION RP_R4N(ir, xlo,xhi)
!  Return a "double precision" real number randomly selected
!  from the interval (xlo,xhi).
!
!     INPUT:
!         IR: random number generator index             [I4]
!        XLO: lower limit on the sampling interval.     [R4]
!        XHI: upper limit on the sampling interval.     [R4]
!
!     OUTPUT:
!         RP_R4N: random number in (0,1)                [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       26-DEC-1997 RD STEWART
!          Add trap for bad generator indices.
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS:
!       This routine should work correctly on any computer
!       with a maximum integer greater than or equal
!       to 2**31 - 1.  NOTE: The proper function of
!       the routine can verified by checking to see
!       that the random number seed after the generation
!       of 10000 random numbers is 1,043,618,065
!
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ir,jr
REAL(R4K) xlo,xhi
REAL(R8K) tmp,seed
REAL(R4K) RP_R4N
INTENT(in) :: ir,xlo,xhi

jr=MAX(1,(MIN(nrg,ir)))  !trap bad generator indices
seed = isd(jr)
tmp = A*seed
seed = tmp - M*DINT(tmp*Minv)
RP_R4N = (xhi-xlo)*seed*Minv + xlo
isd(jr) = seed                !update seed

RETURN
END

FUNCTION RP_U8N(ir)
!  Return a "double precision" real number randomly selected
!  from the interval (0,1).
!
!     INPUT:
!         IR: random number generator index             [I4]
!
!     OUTPUT:
!         RP_U8N: random number in (0,1)                [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       26-DEC-1997 RD STEWART
!          Add trap for bad generator indices.
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS:
!       This routine should work correctly on any computer
!       with a maximum integer greater than or equal
!       to 2**31 - 1.  NOTE: The proper function of
!       the routine can verified by checking to see
!       that the random number seed after the generation
!       of 10000 random numbers is 1,043,618,065
!
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE RANPAK_WORKSPACE
INTEGER(I4K) ir,jr
REAL(R8K) tmp,seed
REAL(R8K) RP_U8N
INTENT(in) :: ir

jr=MAX(1,(MIN(nrg,ir)))  !trap bad generator indices
seed = isd(jr)
tmp = A*seed
seed = tmp - M*DINT(tmp*Minv)
RP_U8N = seed*Minv
isd(jr) = seed                !update seed

RETURN
END

SUBROUTINE RP_CHECK(iflg)
!  Verify that all of the core Lehmer random number generators
!  used in RANPAK working correctly on the current hardware/software
!  configuration.
!
!  INPUT:
!
!  OUTPUT:
!     IFLG: task completion flag                        [I4]
!           0: all generators are in working order
!          <0: one or more generators failed
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, revised RANCHK report format.
!       31-JAN-1999 RD STWEART
!          Updated to Fortran 90 version of RANPAK
!
!     COMMENTS:
!       The proper function of core generators is verified
!       by checking to see that the random number seed
!       after 10,000 calls is 1,043,618,065.
!
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IFLG
INTEGER(I4K) IREF,NTRY
PARAMETER(NTRY=10000,IREF=1043618065)
REAL(R4K)    RX
REAL(R8K)    DX
INTEGER(I4K) IX,I,ISEED
REAL(R8K) RP_U8N,RP_R8N
REAL(R4K) RP_R4N
INTEGER(I4K) RP_INT
EXTERNAL RP_U8N,RP_R8N,RP_R4N,RP_INT
INTENT(out) :: iflg


! SET ALL RANDOM NUMBER SEEDS INITIALLY TO UNITY.
CALL RP_SETRNS(0,1)

! Test the core RANPAK generators.
IFLG=0
DO i=1,ntry
   dx = RP_U8N(1)
   dx = RP_R8N(2,ZERO,ONE)
   rx = RP_R4N(3,0.0,1.0)
   ix = RP_INT(4,-10,10)
ENDDO
CALL RP_GETRNS(1, iseed)
IF (ISEED.NE.IREF) IFLG=-1
CALL RP_GETRNS(2, iseed)
IF (ISEED.NE.IREF) IFLG=IFLG-10
CALL RP_GETRNS(3, iseed)
IF (ISEED.NE.IREF) IFLG=IFLG-100
CALL RP_GETRNS(3, iseed)
IF (ISEED.NE.IREF) IFLG=IFLG-1000

RETURN
END

SUBROUTINE RP_STAT(io,niter, iflg)
!  Conduct some performance tests as an additional means of
!  verifying that RANPAK generators are working correctly on
!  the a computer system.
!
!  INPUT:
!       IO: file handle for a valid output device           [I4]
!    NITER: number of iterations used in performance tests  [I4]
!
!  OUTPUT:
!     IFLG: task completion flag                            [I4]
!           0: all generators are in working order
!          <0: one or more generators failed
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-DEC-1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, revised RANCHK report format.
!       31-JAN-1999 RD STWEART
!          Updated to Fortran 90 version of RANPAK
!
!     COMMENTS:
!       For more information on random number generators,
!       refer to the classic paper
!
!       Park, S.K. and Miller, K.W., "Random number generators:
!       good ones are hard to find."  Communications of the
!       ACM, Vol 31, No 10 (Oct. 1988).
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IO,IFLG
INTEGER(I4K) NITER,ITER
REAL(R8K)    DX,BIN(-4:5),XAVG,X2
INTEGER(I4K) IX,I
REAL(R8K) RP_U8N
INTEGER(I4K) RP_INT
EXTERNAL RP_INT,RP_U8N
INTENT(in) :: io,niter
INTENT(out) :: iflg

! SET ALL RANDOM NUMBER SEEDS INITIALLY TO UNITY.
CALL RP_SETRNS(0,987654321)

ITER=MAX(1000,NITER)
BIN=zero
xavg=zero
x2=zero
DO i=1,iter
   dx = RP_U8N(1)
   XAVG=XAVG+DX  !1st moment (should be one half)
   X2=X2+DX*DX   !2nd moment (should be one third)
   ix = RP_INT(4,-4,5)
   BIN(ix)=BIN(ix)+one
ENDDO
XAVG=XAVG/ITER
X2=X2/ITER
WRITE(unit=io,fmt='(A,1PG15.6)') '1st moment: ',XAVG
WRITE(unit=io,fmt='(A,1PG15.6)') '2nd moment: ',X2

BIN=BIN/ITER
WRITE(unit=io,fmt='(A)',iostat=iflg) 'RP_INT BIN TEST: '
DO i=-4,5
  WRITE(unit=io,fmt='(I3,1PG15.6)',iostat=iflg) i,bin(i)
ENDDO
WRITE(unit=io,fmt='(A,1PG15.6)',iostat=iflg)    &
    'EA. BIN SHOULD BE APPROXIMATELY: ',ONE/ten
iflg=0

RETURN
END

SUBROUTINE RP_R4PLANE(ir,w,h, x,y)
!  Select x-y cartesian coordinates at random from a
!  finite-sized plane of width W and height H.
!
!    INPUT:
!       IR: random number generator index.            [I4]
!        W: width of the plane                        [R4]
!        H: height of the plane                       [R4]
!
!    OUTPUT:
!        X: x-coordinate in [-w/2,w/2]                [R4]
!        Y: y-coordinate in [-h/2,h/2]                [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, re-wrote algorithm to use the external generator
!          R4X instead of in-line generator.
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: Recommend use of generic RP_PLANE routine.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) ir
REAL(R4K) w,h,x,y
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,w,h
INTENT(out) :: x,y

      x = w*(RP_U8N(ir) - half)
      y = h*(RP_U8N(ir) - half)

RETURN
END

SUBROUTINE RP_R8PLANE(ir,w,h, x,y)
!  Select x-y cartesian coordinates at random from a
!  finite-sized plane of width W and height H.
!
!    INPUT:
!       IR: random number generator index.            [I4]
!        W: width of the plane                        [R8]
!        H: height of the plane                       [R8]
!
!    OUTPUT:
!        X: x-coordinate in [-w/2,w/2]                [R8]
!        Y: y-coordinate in [-h/2,h/2]                [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, re-wrote algorithm to use the external generator
!          R8X instead of in-line generator.
!       29-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: Recommend use of generic RP_PLANE routine.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) ir
REAL(R8K) w,h,x,y
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,w,h
INTENT(out) :: x,y

      x = w*(RP_U8N(ir) - half)
      y = h*(RP_U8N(ir) - half)

RETURN
END

SUBROUTINE RP_R8ISO(ir, u,v,w)
! Select a set of random direction cosines from all 4PI
! steradians.  U,V, and W are returned as the x,y, and z
! components of the direction vector.
!
! INPUT:
!   IR: random number generator index.                   [I4]
!
! OUTPUT:
!       U: direction cosine for x-axis [-1,1]            [R8]
!       V: direction cosine for y-axis [-1,1]            [R8]
!       W: direction cosine for z-axis [-1,1]            [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       22-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.

!     COMMENTS:
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(i4K) ir
REAL(R8K) MU,PHI,U,V,W
REAL(R8K) sint
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir
INTENT(out) :: u,v,w

      mu = TWO*RP_U8N(ir) - ONE
      sint = SQRT(ONE-mu*mu)
      phi = TWOPI*RP_U8N(ir)
      u = COS(phi)*sint
      v = SIN(phi)*sint
      w = mu

RETURN
END

SUBROUTINE RP_R4ISO(ir, u,v,w)
! Select a set of random direction cosines from all 4PI
! steradians.  U,V, and W are returned as the x,y, and z
! components of the direction vector.
!
! INPUT:
!   IR: random number generator index.                   [I4]
!
! OUTPUT:
!       U: direction cosine for x-axis [-1,1]            [R8]
!       V: direction cosine for y-axis [-1,1]            [R8]
!       W: direction cosine for z-axis [-1,1]            [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       22-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.
!
!     COMMENTS:
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(i4K) ir
REAL(R4K) MU,PHI,U,V,W
REAL(R4K) sint
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir
INTENT(out) :: u,v,w

      mu = TWO*RP_U8N(ir) - ONE
      sint = SQRT(ONE-mu*mu)
      phi = TWOPI*RP_U8N(ir)
      u = COS(phi)*sint
      v = SIN(phi)*sint
      w = mu

RETURN
END

SUBROUTINE RP_R8DISK(ir,r, px,py)
!  Select x-y cartesian coordinates at random from a
!  disk of radius R centered about the origen.
!
!  INPUT:
!      IR: random number generator index     [I4]
!       R: radius of disk                    [R8]
!
!  OUTPUT:
!      PX: x-coordinate [0,R]                [R8]
!      PY: y-coordinate [0,R]                [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       22-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: The radius of the disk R *must* be greater than
!               zero or an untrapped error will occur.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) R,PX,PY
REAL(R8K) phi,rad
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,r
INTENT(out) :: px,py

  phi = TWOPI*RP_U8N(ir)
  rad = R*SQRT(RP_U8N(ir))
  px = rad*COS(phi)
  py = rad*SIN(phi)

RETURN
END

SUBROUTINE RP_R4DISK(ir,r, px,py)
!  Select x-y cartesian coordinates at random from a
!  disk of radius R centered about the origen.
!
!  INPUT:
!      IR: random number generator index     [I4]
!       R: radius of disk                    [R4]
!
!  OUTPUT:
!      PX: x-coordinate [0,R]                [R4]
!      PY: y-coordinate [0,R]                [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       22-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: The radius of the disk R *must* be greater than
!               zero or an untrapped error will occur.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) R,PX,PY
REAL(R4K) phi,rad
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,r
INTENT(out) :: px,py

  phi = TWOPI*RP_U8N(ir)
  rad = R*SQRT(RP_U8N(ir))
  px = rad*COS(phi)
  py = rad*SIN(phi)

RETURN
END

SUBROUTINE RP_R8RING(ir,r1,r2, px,py)
!  Select x-y cartesian coordinates at random from a
!  ring with an inner radius R1 and outer radius R2.
!  The ring is centered about the origen.
!
!  INPUT:
!    IR: index to a random number generator.       [I4]
!    R1: inner radius of ring                      [R8]
!    R2: outer radius of ring                      [R8]
!
!  OUTPUT:
!    PX: x-coordinate [R1,R2]                      [R8]
!    PY: y-coordinate [R1,R2]                      [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.
!
!     COMMENTS: R1 and R2 *must* be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) r1,r2,PX,PY
REAL(R8K) phi,rad
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,r1,r2
INTENT(out) :: px,py

  phi = TWOPI*RP_U8N(ir)
  px = R1*R1              !used as dummy variable
  rad = SQRT(RP_U8N(ir)*(R2*R2-px) + px)
  px = rad*COS(phi)
  py = rad*SIN(phi)

RETURN
END

SUBROUTINE RP_R4RING(ir,R1,R2, px,py)
!  Select x-y cartesian coordinates at random from a
!  ring with an inner radius R1 and outer radius R2.
!  The ring is centered about the origen.
!
!  INPUT:
!    IR: index to a random number generator.       [I4]
!    R1: inner radius of ring                      [R4]
!    R2: outer radius of ring                      [R4]
!
!  OUTPUT:
!    PX: x-coordinate [R1,R2]                      [R4]
!    PY: y-coordinate [R1,R2]                      [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.
!
!     COMMENTS: R1 and R2 *must* be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) R1,R2,PX,PY
REAL(R4K) phi,rad
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,r1,r2
INTENT(out) :: px,py

  phi = TWOPI*RP_U8N(ir)
  px = R1*R1              !used as dummy variable
  rad = SQRT(RP_U8N(ir)*(R2*R2-px) + px)
  px = rad*COS(phi)
  py = rad*SIN(phi)

RETURN
END

SUBROUTINE RP_R8SPH(ir,r, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  sphere with radius R centered about the origen.
!
!  INPUT:
!    IR: random number generator index             [I4]
!    R: radius of sphere                           [R8]
!
!  OUTPUT:
!    PX: x-coordinate [0,R]                        [R8]
!    PY: y-coordinate [0,R]                        [R8]
!    PZ: z-coordinate [0,R]                        [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.
!
!     COMMENTS: The radius *must* be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) R,PX,PY,PZ
REAL(R8K) rad,u,v,w
REAL(R8K) RP_U8N
EXTERNAL RP_U8N,RP_R8ISO
INTENT(in) :: ir,r
INTENT(out) :: px,py,pz

!  Select an unit direction vector at random
!  from all 4PI steradians.
   CALL RP_R8ISO(ir, u,v,w)

!  Select Radius
   rad = R*(RP_U8N(ir)**THIRD)
   px = rad*u
   py = rad*v
   pz = rad*w

RETURN
END

SUBROUTINE RP_R4SPH(ir,r, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  sphere with radius R centered about the origen.
!
!  INPUT:
!    IR: random number generator index             [I4]
!    R: radius of sphere                           [R8]
!
!  OUTPUT:
!    PX: x-coordinate [0,R]                        [R8]
!    PY: y-coordinate [0,R]                        [R8]
!    PZ: z-coordinate [0,R]                        [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90.
!
!     COMMENTS: The radius *must* be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) R,PX,PY,PZ
REAL(R4K) rad,u,v,w
REAL(R8K) RP_U8N
EXTERNAL RP_U8N,RP_R4ISO
INTENT(in) :: ir,r
INTENT(out) :: px,py,pz

!  Select an unit direction vector at random
!  from all 4PI steradians.
   CALL RP_R4ISO(ir, u,v,w)

!  Select Radius
   rad = R*(RP_U8N(ir)**THIRD)
   px = rad*u
   py = rad*v
   pz = rad*w

RETURN
END

SUBROUTINE RP_R8SHL(ir,R1,R2, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  spherical shell with inner radius R1 and outer radius
!  R2.  The spherical shell is centered about the origen.
!
!  INPUT:
!    IR: random number generator index.            [I4]
!    R1: inner radius of sphere                    [R8]
!    R2: outer radius of sphere                    [R8]
!
!  OUTPUT:
!    PX: x-coordinate [R1,R2]                      [R8]
!    PY: y-coordinate [R1,R2]                      [R8]
!    PZ: z-coordinate [R1,R2]                      [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: R1 must be less than or equal to R2
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) R1,R2,PX,PY,PZ
REAL(R8K) rad,u,v,w
REAL(R8K) RP_U8N
EXTERNAL RP_U8N,RP_R8ISO
INTENT(in) :: ir,r1,r2
INTENT(out) :: px,py,pz

!  Select an unit direction vector at random
!  from all 4PI steradians.
   CALL RP_R8ISO(ir, u,v,w)

!  Select Radius
   px = R1*R1*R1    !used here as a dummy variable
   rad = (RP_U8N(ir)*(R2*R2*R2 - px) + px)**THIRD
   px = rad*u
   py = rad*v
   pz = rad*w

RETURN
END

SUBROUTINE RP_R4SHL(ir,R1,R2, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  spherical shell with inner radius R1 and outer radius
!  R2.  The spherical shell is centered about the origen.
!
!  INPUT:
!    IR: random number generator index.            [I4]
!    R1: inner radius of sphere                    [R4]
!    R2: outer radius of sphere                    [R4]
!
!  OUTPUT:
!    PX: x-coordinate [R1,R2]                      [R4]
!    PY: y-coordinate [R1,R2]                      [R4]
!    PZ: z-coordinate [R1,R2]                      [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: R1 must be less than or equal to R2
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) R1,R2,PX,PY,PZ
REAL(R4K) rad,u,v,w
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,r1,r2
INTENT(out) :: px,py,pz

!  Select an unit direction vector at random
!  from all 4PI steradians.
   CALL RP_R4ISO(ir, u,v,w)

!  Select Radius
   px = R1*R1*R1    !used here as a dummy variable
   rad = (RP_U8N(ir)*(R2*R2*R2 - px) + px)**THIRD
   px = rad*u
   py = rad*v
   pz = rad*w

RETURN
END

SUBROUTINE RP_R8CYL(ir,R,H, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  cylinder with radius R and height H.  The center of
!  the cylinder is at (0,0,0), and the axis of the cylinder
!  is parallel to the z-axis.
!
!  INPUT:
!    IR: random number generator index.             [I4]
!     R: radius of cylinder                         [R8]
!     H: height of cylinder                         [R8]
!
!  OUTPUT:
!    PX: x-coordinate [0,R]                         [R8]
!    PY: y-coordinate [0,R]                         [R8]
!    PZ: z-coordinate [-H/2,+h/2]                   [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: R and H must be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) R,H,PX,PY,PZ
REAL(R8K) zlo,zhi
REAL(R8K) RP_R8N
EXTERNAL RP_R8N
INTENT(in) :: ir,r,h
INTENT(out) :: px,py,pz

  zhi = half*H
  zlo = -zhi
  pz = RP_R8N(ir,zlo,zhi)
  CALL RP_R8DISK(ir,r, px,py)

RETURN
END

SUBROUTINE RP_R4CYL(ir,R,H, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  cylinder with radius R and height H.  The center of
!  the cylinder is at (0,0,0), and the axis of the cylinder
!  is parallel to the z-axis.
!
!  INPUT:
!    IR: random number generator index.             [I4]
!     R: radius of cylinder                         [R4]
!     H: height of cylinder                         [R4]
!
!  OUTPUT:
!    PX: x-coordinate [0,R]                         [R4]
!    PY: y-coordinate [0,R]                         [R4]
!    PZ: z-coordinate [-H/2,+h/2]                   [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: R and H must be greater than or equal to zero.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) R,H,PX,PY,PZ
REAL(R4K) zlo,zhi
REAL(R4K) RP_R4N
EXTERNAL RP_R4N
INTENT(in) :: ir,r,h
INTENT(out) :: px,py,pz

  zhi = half*H
  zlo = -zhi
  pz = RP_R4N(ir,zlo,zhi)
  CALL RP_R4DISK(ir,r, px,py)

RETURN
END

SUBROUTINE RP_R8BOX(ir,WX,WY,WZ, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  box with inside dimensions WX, WY, WZ.  The box is
!  centered about the origen.
!
!  INPUT:
!    IR: random number generator index.            [I4]
!    WX: x-dimension of box                        [R8]
!    WY: y-dimension of box                        [R8]
!    WZ: z-dimension of box                        [R8]
!
!  OUTPUT:
!    PX: x-coordinate [-WX/2,+WX/2]                [R8]
!    PY: y-coordinate [-WY/2,+WY/2]                [R8]
!    PZ: z-coordinate [-WZ/2,+WZ/2]                [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS:
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) WX,WY,WZ,PX,PY,PZ
REAL(R8K) zlo,zhi
REAL(R8K) RP_R8N
EXTERNAL RP_R8N
INTENT(in) :: ir,wx,wy,wz
INTENT(out) :: px,py,pz

  zhi = half*WX
  zlo = -zhi
  px = RP_R8N(ir,zlo,zhi)

  zhi = half*WY
  zlo = -zhi
  py = RP_R8N(ir,zlo,zhi)

  zhi = half*WZ
  zlo = -zhi
  pz = RP_R8N(ir,zlo,zhi)

RETURN
END

SUBROUTINE RP_R4BOX(ir,WX,WY,WZ, px,py,pz)
!  Select x-y-z cartesian coordinates at random from a
!  box with inside dimensions WX, WY, WZ.  The box is
!  centered about the origen.
!
!  INPUT:
!    IR: random number generator index.            [I4]
!    WX: x-dimension of box                        [R4]
!    WY: y-dimension of box                        [R4]
!    WZ: z-dimension of box                        [R4]
!
!  OUTPUT:
!    PX: x-coordinate [-WX/2,+WX/2]                [R4]
!    PY: y-coordinate [-WY/2,+WY/2]                [R4]
!    PZ: z-coordinate [-WZ/2,+WZ/2]                [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-DEC-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS:
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) WX,WY,WZ,PX,PY,PZ
REAL(R4K) zlo,zhi
REAL(R4K) RP_R4N
EXTERNAL RP_R4N
INTENT(in) :: ir,wx,wy,wz
INTENT(out) :: px,py,pz

  zhi = half*WX
  zlo = -zhi
  px = RP_R4N(ir,zlo,zhi)

  zhi = half*WY
  zlo = -zhi
  py = RP_R4N(ir,zlo,zhi)

  zhi = half*WZ
  zlo = -zhi
  pz = RP_R4N(ir,zlo,zhi)

RETURN
END

SUBROUTINE RP_R8ICC(ir,mu, u,v,w)
! Select a set of random (isotropic) direction cosines from the
! "ice cream cone" sub-tending the angle MU=COS(theta) centered
! about the z-axis.  U,V, and W are returned as the x,y, and z
! components of the unit direction vector.
!
!  INPUT:
!    IR: random number generator index.              [I4]
!    MU: dot product of z-hat axis and a             [R8]
!        vector with polar angle theta.
!        NOTE: MU=COS(theta).
!
!  OUTPUT:
!    U: direction cosine for x-axis [-1,1]           [R8]
!    V: direction cosine for y-axis [-1,1]           [R8]
!    W: direction cosine for z-axis [-1,1]           [R8]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, re-wrote algorithm to use the external generator
!          R8X instead of in-line generator.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: To select a direction vector from all 4PI
!               steradians, set MU = -1.  However, the
!               RP_ISO routine will be slightly faster than
!               the RP_ICC routine.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R8K) MU,U,V,W
REAL(R8K) MUP,SINT,PHI
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,mu
INTENT(out) :: u,v,w

  mup = one - RP_U8N(ir)*(one-mu)
  sint = SQRT(one-mup*mup)
  phi = TWOPI*RP_U8N(ir)
  u = COS(phi)*sint
  v = SIN(phi)*sint
  w = mup

RETURN
END

SUBROUTINE RP_R4ICC(ir,mu, u,v,w)
! Select a set of random (isotropic) direction cosines from the
! "ice cream cone" sub-tending the angle MU=COS(theta) centered
! about the z-axis.  U,V, and W are returned as the x,y, and z
! components of the unit direction vector.
!
!  INPUT:
!    IR: random number generator index.              [I4]
!    MU: dot product of z-hat axis and a             [R4]
!        vector with polar angle theta.
!        NOTE: MU=COS(theta).
!
!  OUTPUT:
!    U: direction cosine for x-axis [-1,1]           [R4]
!    V: direction cosine for y-axis [-1,1]           [R4]
!    W: direction cosine for z-axis [-1,1]           [R4]
!
!      LIBRARY: [RANPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-DEC-1994
!     REVISIONS HISTORY:
!       23-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions.
!          Also, re-wrote algorithm to use the external generator
!          R8X instead of in-line generator.
!       31-JAN-1999 RD STEWART
!          Updated to Fortran 90
!
!     COMMENTS: To select a direction vector from all 4PI
!               steradians, set MU = -1.  However, the
!               RP_ISO routine will be slightly faster than
!               the RP_ICC routine.
!
USE NUMERIC_KINDS
USE MATH_USEFUL_DP
INTEGER(I4K) IR
REAL(R4K) MU,U,V,W
REAL(R4K) MUP,SINT,PHI
REAL(R8K) RP_U8N
EXTERNAL RP_U8N
INTENT(in) :: ir,mu
INTENT(out) :: u,v,w

  mup = one - RP_U8N(ir)*(one-mu)
  sint = SQRT(one-mup*mup)
  phi = TWOPI*RP_U8N(ir)
  u = COS(phi)*sint
  v = SIN(phi)*sint
  w = mup

RETURN
END

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
END

SUBROUTINE FN_ERASE
!  Erase all filename information from active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 1]
!
!     INPUT:
!
!     OUTPUT:
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 01-FEB-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE STRPAK
USE FIOPAK_WORKSPACE

  CALL STR_FILL(' ',fni)
  CALL STR_FILL(' ',fno)
  CALL STR_FILL(' ',buf)
  CALL STR_FILL(' ',cwd)
  CALL STR_FILL(' ',cwd_dl)
  CALL STR_FILL(' ',fn_nop)
  CALL STR_FILL(' ',fn_base)
  CALL STR_FILL(' ',fn_ext)

RETURN
END


SUBROUTINE FILE_DELETE(fn, iflg)
!     Delete file=fn if it exists.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!         FN: name of a file.     [A*]
!
!     OUTPUT:
!       IFLG: status flag  [I4]
!             -2: task failed for an unknown reason.
!             -1: task failed because file does not exist
!              0: task successful.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 22-NOV-1997
!     REVISIONS HISTORY:
!
!       29-NOV-1998 RD STEWART
!          Updated to meet F90 language standard and renamed
!          to better describe purpose (formely called DELFN).
!
!     COMMENTS:
!
  CHARACTER(LEN=*) fn
  INTEGER io,iflg,jflg
  LOGICAL exists
  INTENT(IN) :: FN
  INTENT(OUT) :: IFLG

  INQUIRE(file=fn,EXIST=exists,iostat=jflg)
  IF (exists) THEN
!   Delete file
    CALL CHANNEL(io)
    OPEN(unit=io,file=fn,iostat=jflg,status='old')
    CLOSE(unit=io,status='DELETE',iostat=jflg)
    IF (jflg.NE.0) THEN
      iflg=-2
    ELSE
      iflg=0
    ENDIF
  ELSE
!   File does not exist.
    iflg=-1
  ENDIF

RETURN
END

SUBROUTINE FN_GETFNB(fnb)
!  Return as FNB the base (root) part of the filename in active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!
!     OUTPUT:
!        FNB: name of a file extension        [A*]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 04-MAY-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) fnb
INTENT(OUT) :: FNB

CALL STR_FILL(' ',fnb)
CALL STR_ADD(FN_base,fnb)

RETURN
END


SUBROUTINE FN_GETFNX(ext)
!  Return as EXT the extension of the filename in active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!
!     OUTPUT:
!        ext: name of a file extension        [A*]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 07-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) ext
INTENT(OUT) :: EXT

CALL STR_FILL(' ',ext)
CALL STR_ADD(FN_ext,ext)

RETURN
END

SUBROUTINE FN_SAVE(fn, iflg)
!  Save filename=FN to active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        FN: name of an input file.         [A*]
!
!     OUTPUT:
!       IFLG: status flag                   [I4]
!              0: task successful
!             -1: insufficient workspace
!                 allocated for filename.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!         05-MAR-1998 RD STEWART
!           Modified algorithm to handle filenames without
!           an extension.
!         01-FEB-1999 RD STEWART
!           Updated to Fortran 90 version
!
!     COMMENTS: Routine may produce unpredictable results if
!               directory names are not terminated with a
!               '\' or '/'.
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) FN
CHARACTER(LEN=1) DD
INTEGER(I4K) IFLG,NXT,ILEN,JLEN
INTENT(in) :: fn
INTENT(out) :: iflg


  IF (LEN(fn).GT.MXFN) THEN
    IFLG=-1
    RETURN
  ELSE
    IFLG=0
  ENDIF
  CALL FN_ERASE  !Erase current info
  CALL STR_ADD(fn, fni) !transfer to active memory
  CALL STR_ADD(fn, buf) !copy to workspace buffer

  !Analyze file and path.
  UNIX=((INDEX(fni,':').EQ.0).AND.(INDEX(fni,'/').GT.0))
  IF (unix) THEN
     CALL CHR_SAR(CHAR(92),'/',buf)
     DD='/'
  ELSE
    CALL CHR_SAR('/',CHAR(92), buf)
    DD=CHAR(92)
  ENDIF
  CALL STR_XLET(buf)  !interchange letters to simplify analysis

  nxt=1
  IF (.NOT.((BUF(1:1).EQ.'/').OR.(BUF(1:1).EQ.CHAR(92)))) THEN
    !attempt to locate filename extension
    ilen=INDEX(buf,'.')
    jlen=INDEX(buf,dd)
    IF (jlen.EQ.0) THEN
      jlen=LENTRIM(buf)
    ENDIF
    IF (ilen.GT.0) THEN
      IF (ilen.LT.jlen) THEN
        CALL TOKEN(0,'.',buf, nxt, fn_ext,iflg)
        CALL STR_XLET(fn_ext)
      ENDIF
    ENDIF
    !Attempt to locate filename base
    CALL TOKEN(0,'\/:',buf,nxt, fn_base,iflg)
    CALL STR_XLET(fn_base)
    CALL CHR_DELETE('.',fn_base)

    !Complete filename without path
    CALL STR_ADD(FN_base,FN_nop)
    IF (LENTRIM(fn_ext).GT.0) THEN
      !Add the filename extension.
      CALL STR_ADD('.',fn_nop)
      CALL STR_ADD(fn_ext,fn_nop)
    ENDIF
  ENDIF

 !current working directory
 cwd=buf(nxt:)
 CALL STR_XLET(cwd)
 ilen=INDEX(cwd,':')
 IF (ilen.GT.0) THEN
   IF (unix) THEN
     CALL STR_SHFT(-ilen,cwd) !Erase drive letter info from CWD
   ELSE
     !look for drive letter
     CALL STR_FILL(' ',buf)
     CALL STR_ADD(cwd,buf)
     CALL STR_SHFT(-ilen,cwd) !Erase drive letter info from CWD
     CALL STRIP(buf(:ilen))
     ilen=INDEX(buf,':')
     IF (ilen.GT.1) THEN
       ilen=ilen-1
       CWD_DL=buf(ilen:ilen)
     ENDIF

   ENDIF
 ENDIF

RETURN
END

SUBROUTINE FN_GETCWD(dn)
!    Return as DN the name of the current working
!    directory from active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!
!     OUTPUT:
!        DN: directory name                    [A*]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-NOV-1997
!     REVISIONS HISTORY:
!       28-DEC-1997 RD STEWART
!         check for existing '/' or CHAR(92) in directory name
!         before appending one to the end of DN.
!       20-FEB-1998 RD STEWART
!         change algorithm so that a '/' or CHAR(92) is not added to the
!         directory name when the directory name length = 0.
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) DN
CHARACTER(LEN=1) DD
INTENT(out) :: DN

  CALL STR_FILL(' ',DN)
  IF (unix) THEN
    DD='/'
  ELSE
    DD=CHAR(92)
  ENDIF

  IF (.not.unix) THEN
    IF (LENTRIM(CWD_DL).GT.0) THEN
      CALL STR_ADD(CWD_DL//':',DN)
    ENDIF
  ENDIF
  CALL STR_ADD(cwd,dn)

RETURN
END

SUBROUTINE FN_GETCFN(iopt, fn)
!  Use the filename information stored in active memory
!  to assemble a filename with/without the directory
!  path.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!         IOPT: task options                   [I4]
!              -1: filename w/o path
!       otherwise: path + filename
!
!     OUTPUT:
!        FN: name of file                      [A*]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) IOPT
CHARACTER(LEN=*) FN
INTENT(in) :: iopt
INTENT(out) :: FN

  CALL STR_FILL(' ',FN)
  IF (iopt.NE.(-1)) THEN
    CALL FN_GETCWD(fn)  !get directory name
  ENDIF
  CALL STR_ADD(fn_nop,fn)

RETURN
END

SUBROUTINE FN_SETCWD(dn)
!  Set the current working directory to FN
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        DN: folder/directory name             [A*]
!
!     OUTPUT:
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) DN
CHARACTER(LEN=1) DD
INTEGER(I4K) ILEN
INTENT(in) :: dn

CALL STR_FILL(' ',CWD)
CALL STR_FILL(' ',CWD_DL)
CALL STR_ADD(DN,CWD)

!Analyze file and path.
UNIX=((INDEX(cwd,':').EQ.0).AND.(INDEX(cwd,'/').GT.0))
IF (unix) THEN
  CALL CHR_SAR(CHAR(92),'/', cwd)
  DD='/'
ELSE
  CALL CHR_SAR('/',CHAR(92), cwd)
  DD=CHAR(92)
ENDIF
ilen=LENTRIM(cwd)
IF (ilen.EQ.0) THEN
  return
ELSEIF (cwd(ilen:ilen).NE.DD) THEN
  CALL STR_ADD(DD,cwd)
ENDIF

ilen=INDEX(cwd,':')
IF (ilen.GT.0) THEN
  IF (unix) THEN
    CALL STR_SHFT(-ilen,cwd) !Erase drive letter info from CWD
  ELSE
    !look for drive letter
    CALL STR_FILL(' ',buf)
    CALL STR_ADD(cwd,buf)
    CALL STR_SHFT(-ilen,cwd) !Erase drive letter info from CWD
    CALL STRIP(buf(:ilen))
    ilen=INDEX(buf,':')
    IF (ilen.GT.1) THEN
      ilen=ilen-1
      CWD_DL=buf(ilen:ilen)
    ENDIF
  ENDIF
ENDIF

RETURN
END

SUBROUTINE FN_SETDL(dl)
!  Set the current drive letter to dl.  If DL is any
!  upper- or lowercase letter, then the all instances
!  of '/' in the current working directory are converted
!  to CHAR(92).
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        DL: drive letter             [A1]
!
!     OUTPUT:
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=1) DL
INTEGER(I4K) ilen
INTENT(in) :: dl

IF (.not.CHR_ILET(dl)) THEN
  !invalid drive letter
  RETURN
ENDIF
CWD_DL=DL
CALL CHR_SAR('/',CHAR(92), cwd)
unix=no

ilen=INDEX(cwd,CHAR(92))
DO WHILE(ilen.GT.0)
  IF (cwd(ilen:ilen).EQ.'.') THEN
    CALL STR_SHFT(-1,cwd(ilen:))
  ENDIF
  ilen=ilen-1
ENDDO

RETURN
END

SUBROUTINE FN_SETFNB(base)
!  Set the filename base to base.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        base: name of a file             [A*]
!
!     OUTPUT:
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) base
INTENT(in) :: base

CALL STR_FILL(' ',FN_base)
CALL STR_FILL(' ',FN_nop)
CALL STR_ADD(base,FN_base)

!Complete filename without path
CALL STR_ADD(FN_base,FN_nop)
IF (LENTRIM(fn_ext).GT.0) THEN
  !Add the filename extension.
  CALL STR_ADD('.',fn_nop)
  CALL STR_ADD(fn_ext,fn_nop)
ENDIF

RETURN
END

SUBROUTINE FN_SETFNX(ext)
!  Set the filename base to FN.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        ext: name of a file extension        [A*]
!
!     OUTPUT:
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 10-NOV-1997
!     REVISIONS HISTORY:
!       01-FEB-1999 RD STEWART
!         Updated to Fortran 90 version
!       24-MAR=1999 RD STEWART
!         Delete '.' from the input filename extension.
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
CHARACTER(LEN=*) ext
INTENT(in) :: ext

CALL STR_FILL(' ',FN_ext)
CALL STR_FILL(' ',FN_nop)
CALL STR_ADD(ext,FN_ext)
CALL CHR_DELETE('.',FN_ext) !remove all periods from FN_ext

!Complete filename without path
CALL STR_ADD(FN_base,FN_nop)
IF (LENTRIM(fn_ext).GT.0) THEN
  !Add the filename extension.
  CALL STR_ADD('.',fn_nop)
  CALL STR_ADD(fn_ext,fn_nop)
ENDIF

RETURN
END

SUBROUTINE FIO_LOCSTR(io,opt,sstr, nlp,ostr,iflg)
!  Search through the ASCII file open on unit=io for the
!  first line containing a match for the characters in
!  SSTR.  Trailing blanks in SSTR are ignored.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        IO: device number for an ASCII file                     [I4]
!       OPT: three digit number indicating task option           [A3]
!             0xx: begin search from current location (default)
!             1xx: rewind file before performing search
!             x0x: ignore extra spaces between words (default)
!             x1x: all spaces are ignored.
!             x2x: all spaces are significant
!             xx0: case insensitive (default)
!             xx1: case sensitive
!             EXAMPLE: OPT = '021' = rewind file then perform
!                      a case and space sensitve search.
!      SSTR:  character variable                                [A*]
!
!     OUTPUT:
!        NLP: number of lines successfully processed            [I4]
!       OSTR: last line of text read from unit=io.              [A*]
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: search string not found
!             -2: FIOPAK workspace insufficient to store SSTR
!             -3: unit=io not ready for file i/o
!             -4: task failed while attempting to rewind unit=io.
!             -5: task failed while performing file i/o.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 24-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: The most likely cause of a IFLG=-5 error is
!               that the end-of-file was reached without finding
!               the search string SSTR.
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) IO,NLP,N,J1,J2,IFLG
CHARACTER(LEN=3) OPT
CHARACTER(LEN=*) SSTR,OSTR
LOGICAL ok
INTENT(in) :: io,opt,sstr
INTENT(out) :: nlp,ostr,iflg

IF (LENTRIM(SSTR).GT.LEN(buf)) THEN
  IFLG=-2
  RETURN
ELSE
  INQUIRE(unit=io,OPENED=ok)
  IF (.not.OK) THEN
    IFLG=-3
    RETURN
  ENDIF
ENDIF

IF (OPT(1:1).EQ.'1') THEN
  REWIND(unit=io,iostat=iflg)
  IF (iflg.NE.0) THEN
    IFLG=-4
    RETURN
  ENDIF
ENDIF

!Generate search string
CALL STR_FILL(' ',BUF)
CALL STR_ADD(SSTR,BUF)
IF (OPT(2:2).EQ.'1') THEN
  !ignore all blank spaces
   CALL CHR_DELETE(' ',BUF)
ELSEIF (OPT(2:2).NE.'2') THEN
  !ignore extra blanks between words
  CALL CHR_RMI(' ',' ', BUF)
ENDIF
N=MAX(LENTRIM(buf),1)

NLP=0
DO WHILE (ok)
  READ(unit=io,fmt='(A)',iostat=iflg) OSTR
  IF (iflg.EQ.0) THEN
    NLP=NLP+1

    CALL STR_FILL(' ',FIOBUF)
    CALL STR_ADD(OSTR,FIOBUF)
    IF (OPT(2:2).EQ.'1') THEN
      !ignore all blank spaces
      CALL CHR_DELETE(' ',FIOBUF)
    ELSEIF (OPT(2:2).NE.'2') THEN
      !ignore extra blanks between words
      CALL CHR_RMI(' ',' ', FIOBUF)
    ENDIF

    IF (OPT(3:3).NE.'1') THEN
      !case insensitive match
      CALL STR_FAM(buf(:n),fiobuf, j1,j2,iflg)
    ELSE
      !case sensitive match
      CALL STR_FXM(buf(:n),fiobuf, j1,j2,iflg)
    ENDIF
    OK=(IFLG.NE.0)

  ELSE
    !terminate search
    CALL STR_FILL(' ',OSTR)
    OK=NO
    IFLG=-5
  ENDIF
ENDDO

RETURN
END

SUBROUTINE FIO_ECHO(ii,io,opt, iflg)
!  Echo the contents of the ASCII file open unit=ii to
!  the ASCII file open on unit=io with optional
!  modifications to format/content of text.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        II: device handle for source file                       [I4]
!        IO: device handle for output file                       [I4]
!       OPT: task options                                        [A4]
!            xxx0: append data to output file (default)
!            xxx1: rewind output file
!            xx0x: all letters left unchanged (default)
!            xx1x: convert all letters to uppercase
!            xx2x: convert all letters to lowercase
!            x0xx: include comments and comment blocks (default)
!            x1xx: delete comments, comment blocks, and blank lines
!            0xxx: leave word spacing unchanged (default)
!            1xxx: delete extra spaces
!
!     OUTPUT:
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: unit=ii not ready for file i/o
!             -2: unit=i0 not ready for file i/o
!             -3: unable to rewind unit=ii
!             -4: unable to rewind unit=io
!             -5: read/write error
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 25-MAR-1999
!     REVISIONS HISTORY:
!       20-APR-1999 RD STEWART
!         Changed algorithm so that blank lines are deleted
!         when OPT = 'x1xx'
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) ii,io,iflg,n
CHARACTER(LEN=4) opt
CHARACTER(LEN=MXLL) STR
LOGICAL again,ok
INTENT(in) :: ii,io,opt
INTENT(out) :: iflg

INQUIRE(unit=ii,OPENED=ok)
IF (.not.OK) THEN
  IFLG=-1
  RETURN
ELSE
  !rewind file
  REWIND(unit=ii,iostat=iflg)
  IF (iflg.NE.0) THEN
    IFLG=-3
    RETURN
  ENDIF
ENDIF
IF (io.NE.6) THEN
  INQUIRE(unit=io,OPENED=ok)
  IF (.not.OK) THEN
    IFLG=-2
    RETURN
  ELSEIF (OPT(4:4).EQ.'1') THEN
    !rewind file
    REWIND(unit=io,iostat=iflg)
    IF (iflg.NE.0) THEN
      IFLG=-4
      RETURN
    ENDIF
  ENDIF
ENDIF

again=YES
DO WHILE (again)
  READ(unit=ii,fmt='(A)',iostat=iflg) STR
  IF (iflg.EQ.0) THEN
    N=MAX(1,LENTRIM(str))
    IF (opt(1:1).EQ.'1') THEN
      !Delete extra spaces
      CALL CHR_RMI(' ',' ', str(:n))
    ENDIF
    IF (opt(3:3).EQ.'1') THEN
      !convert to uppercase
      CALL TRANUC(str(:n))
    ELSEIF (opt(3:3).EQ.'2') THEN
      !convert to lowercase
      CALL TRANLC(str(:n))
    ENDIF

    IF (opt(2:2).EQ.'1') THEN
      !Delete comments and comment blocks
      CALL VI_CHKSTR(str, n,iflg)
      IF ((iflg.LE.0).AND.(n.GT.0)) THEN
        WRITE(unit=io,fmt='(A)',iostat=iflg) str(:n)
      ELSE
        !Skip comment
        iflg=0
      ENDIF
    ELSE
      WRITE(unit=io,fmt='(A)',iostat=iflg) str(:n)
    ENDIF

    IF (iflg.NE.0) THEN
      again=no
      iflg=-5
    ENDIF
  ELSE
    !Normal termination
    again=no
    iflg=0
  ENDIF
ENDDO

RETURN
END

SUBROUTINE FIO_LOCVAR(io,opt,var, nlp,str,iflg)
!  Locate variable VAR in the ASCII file open on unit=io.
!  If the variable search is successful, the data block
!  file i/o option is turned ON.  Also, STR is returned
!  as the line of text containing the desired variable.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        IO: device number for an ASCII file                     [I4]
!       OPT: task options                                        [A1]
!            0: begin search at current location (default)
!            1: rewind unit=io before beginning search
!       VAR: variable name                                       [A*]
!
!     OUTPUT:
!        NLP: number of lines processed                         [I4]
!        STR: last line of text processed                       [A*]
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: unit=io not ready for file i/o
!             -2: file rewind failed
!             -3: variable name not found before end of file.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 25-MAR-1999
!     REVISIONS HISTORY:
!       04-MAY-1999 RD STEWART
!         BUG FIX: Verify that variable name is followed by
!                  a blank space or equal sign to prevent
!                  mis-identification of two variables with
!                  similar names (e.g., VB and VBX).
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) IO,NLP,IFLG,J1,J2,N
CHARACTER(LEN=*) VAR,STR
CHARACTER(LEN=1) OPT
LOGICAL OK,AGAIN
INTENT(in) :: io,opt,var
INTENT(out) :: nlp,str,iflg

DON=NO  !Turn data block off
CON=NO  !Turn comment block off
INQUIRE(unit=io,OPENED=ok)
IF (.not.OK) THEN
  IFLG=-1
  RETURN
ELSEIF (OPT(1:1).EQ.'1') THEN
  !rewind file
  REWIND(unit=io,iostat=iflg)
  IF (iflg.NE.0) THEN
    IFLG=-2
    RETURN
  ENDIF
ENDIF

!CASE INSENSITIVE SEARCH FOR VARIABLE.
again=(iflg.EQ.0)
NLP=0
DO WHILE (again)
  READ(unit=io,FMT='(A)',iostat=iflg) STR
  IF (iflg.NE.0) THEN
    !Variable name not found
    IFLG=-3
    RETURN
  ELSE
    NLP=NLP+1
    CALL VI_CHKSTR(str, n,iflg)
    IF ((n.GT.0).AND.(iflg.EQ.0)) THEN
      CALL STR_FAM(var,str(:n), j1,j2,iflg)
      IF ((j2.LT.n).AND.(j2.GT.0)) THEN
        !Variable name must be followed by either a blank space or equal sign
        AGAIN=((IFLG.EQ.0).AND.STR_IMAT(0,str(j2+1:j2+1),' ='))
        AGAIN=(.NOT.AGAIN)
      ELSE
        AGAIN=(IFLG.LT.0)
      ENDIF
    ENDIF
  ENDIF
ENDDO

IF (iflg.EQ.0) THEN
  CALL STR_FILL(' ',STR(n+1:))  !Erase appended comment
  CALL CHR_RMI(' ',' ', str)    !delete extra spaces
ENDIF
DON=YES  !TURN ON DATA I/O BLOCK

RETURN
END

SUBROUTINE VI_ALLOC(N, IFLG)
!  Increase the amount of memory allocated to to the FIOPAK
!  variable input module so that a minmum of N data elements
!  can be stored.
!
!  INPUT:
!    N: minimum number of data storage units to allocate.    [I4]
!
!  OUTPUT:
!     IFLG: task completion flag                             [I4]
!           0: task successful
!          -1: task failed while allocating memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
USE FIOPAK_WORKSPACE
USE MATH_USEFUL_DP
INTEGER(I4K) I,N,M,IFLG
INTENT(in) :: n
INTENT(out) :: iflg

IFLG=0
IF (N < NDMX) RETURN  !Don't need to allocate any more space yet.

IF (NDMX.LT.1) THEN
  !as yet, no memory has been allocated.
  NDMX=MAX(N,NINC)
  DEALLOCATE(VI_DAT, STAT = IFLG)
  ALLOCATE(VI_DAT(NDMX),STAT = IFLG)
  IF (IFLG.NE.0) THEN
    iflg=-1
    RETURN
  ENDIF
  VI_DAT=zero
  NDES=0
ELSE
  !Need to increase space allocated to dose-interval stack.
  ALLOCATE(VI_TMP(NDMX),STAT = IFLG)
  IF (IFLG.NE.0) THEN
    iflg=-1
    RETURN
  ENDIF

  M = MAX(N,NDMX+NINC)
  DO i=1,NDES
    VI_TMP(i)=VI_DAT(i)
  ENDDO

  DEALLOCATE(VI_DAT, STAT = IFLG)
  ALLOCATE(VI_DAT(m),STAT = IFLG)
  IF (IFLG.NE.0) THEN
    IFLG=-1
    RETURN
  ENDIF
  DO i=1,ndes
    VI_DAT(i)=VI_TMP(i)
  ENDDO
  NDMX=M
  DEALLOCATE(VI_TMP, STAT = IFLG)

ENDIF

RETURN
END

SUBROUTINE VI_CHKSTR(str, n,iflg)
!  Compute properties associated with the line of
!  text in STR.  The comment (CON) block properties
!  associated with file i/o are automatically turned
!  on/off as appropriate.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!       STR: character variable                                 [A*]
!
!     OUTPUT:
!          N: length of STR minus appended comments             [I4]
!       IFLG: status flag                                       [I4]
!              3: STR terminates comment block
!              2: STR initiates comment block
!              1: STR is a comment
!              0: STR may contain alphanumeric data
!             -1: STR terminates data input block
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!       18-MAY-1999 RD STEWART
!         Removed Fortran 77 style comment indicator by changing
!         IF ((buf(1:1).EQ.'!').OR.(buf(1:1).EQ.'C'))
!                      to
!         IF (buf(1:1).EQ.'!') THEN
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) N,IFLG,ILEN,JLEN,KLEN
CHARACTER(LEN=*) STR
INTENT(in) :: str
INTENT(out) :: n,iflg

!PROCESS COMMENT BLOCK
N=0
CALL STR_FILL(' ',BUF)
CALL STR_ADD(STR,BUF)
CALL STRIP(buf)
CALL TRANUC(buf)
ILEN=LENTRIM(buf)
KLEN=ILEN-1
IF (CON) THEN
  !Comment block is currently ON.
  IF (BUF(1:2).EQ.'*/') THEN
    IFLG=3  !this line is a comment
    CON=NO  !end of comment block
    RETURN
  ELSEIF (BUF(1:2).EQ.'!>') THEN
    IFLG=3  !this line is a comment
    CON=NO  !end of comment block
    RETURN
  ELSEIF (klen.GT.0) THEN
    IF (BUF(klen:ilen).EQ.'*/') THEN
      IFLG=3  !this line is a comment
      CON=NO  !end of comment block
      RETURN
    ELSEIF (BUF(klen:ilen).EQ.'!>') THEN
      IFLG=3  !this line is a comment
      CON=NO  !end of comment block
      RETURN
    ELSE
      !This line is part of a comment block
      IFLG=1
      RETURN
    ENDIF
  ELSE
    !This line is part of a comment block
    IFLG=1
    RETURN
  ENDIF
ELSEIF (.NOT.CON) THEN
  !COMMENT BLOCK IS TURNED OFF
  IF (BUF(1:2).EQ.'/*') THEN
    !begin comment block
    CON=YES
    IFLG=2
    RETURN
  ELSEIF (BUF(1:2).EQ.'<!') THEN
    !begin comment block
    CON=YES
    IFLG=2
    RETURN
  ENDIF
ENDIF
IF (buf(1:1).EQ.'!') THEN
  !1-line comment
  IFLG=1
  RETURN
ENDIF

!CHECK FOR END OF DATA INPUT BLOCK CONDITION
IF (ILEN.EQ.0) THEN
  !Blank line terminates data input block
  IFLG=-1
  RETURN
ELSEIF (BUF(1:4).EQ.'&END') THEN
  !End of data input block
  IFLG=-1
  RETURN
ELSEIF (BUF(1:1).EQ.'/') THEN
  !End of data input block
  IFLG=-1
  RETURN
ENDIF

!STR MAY CONTAIN DATA
IFLG=0
JLEN=INDEX(STR,'!')
IF (JLEN.GT.0) THEN
  N=JLEN-1
ELSE
  N=LENTRIM(str)
ENDIF

RETURN
END

SUBROUTINE VI_DEALLOC
!  Deallocate all memory assigned to the variable input module.
!
!  INPUT:
!
!  OUTPUT:
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
USE FIOPAK_WORKSPACE
INTEGER(I4K) IFLG

IF (NDMX.GT.0) THEN
  NDES=0
  NDMX=0
  DEALLOCATE(VI_TMP, STAT = IFLG)
  DEALLOCATE(VI_DAT, STAT = IFLG)
ENDIF

RETURN
END

SUBROUTINE VI_ERASE
!  Erase the contents of the variable input module without
!  deallocating memory.
!
!  INPUT:
!
!  OUTPUT:
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
USE FIOPAK_WORKSPACE
USE MATH_USEFUL_DP

IF (NDMX.GT.0) THEN
  NDES=0
  VI_DAT=zero
ENDIF

RETURN
END

SUBROUTINE VI_GETFN(io,opt,var, fn,iflg)
!  Search through the ASCII file open on unit=io for a character
!  string associated with the input variable VAR.  Filenames
!  must not contain any blank spaces and cannot span multiple
!  lines in the file (i.e., the entire filename must appear
!  on a single line of the ASCII file along with VAR.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        IO: device number for an ASCII file                     [I4]
!       OPT: task options                                        [A1]
!            0: begin search at current location (default)
!            1: rewind unit=io before beginning search
!       VAR: variable name                                       [A*]
!
!     OUTPUT:
!         FN: filename                                          [A*]
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: unit=io not ready for file i/o
!             -2: file rewind failed
!             -3: variable name not found before end of file.
!             -4: filename not found
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER IO,NLP,IFLG,J1,J2,NXT
CHARACTER(LEN=*) VAR,FN
CHARACTER(LEN=1) opt
CHARACTER(LEN=MXLL) STR
CHARACTER(LEN=*), PARAMETER :: LST=' `()={}[]|:;''",<>?'
INTENT(in) :: io,opt,var
INTENT(out) :: fn,iflg

!Locate variable name in unit=io.
CALL FIO_LOCVAR(io,opt,var, nlp,str,iflg)
IF (iflg.NE.0) THEN
  CALL STR_FILL(' ',fn)
  RETURN
ENDIF

!LOCATE POSITION OF KEYWORD IN STR
CALL STR_FAM(var,str, j1,j2,iflg) !case insensitive
IF (iflg.NE.0) THEN
  IFLG=-3
  RETURN
ENDIF

!GET FILENAME ASSOCIATED WITH VAR
nxt=j2+1
CALL TOKEN(0,lst,str,nxt, fn,iflg)  !case insensitive
IF (iflg.NE.0) THEN
  IFLG=-4
ENDIF

RETURN
END

SUBROUTINE VI_GETSTR(io,opt,var, str,iflg)
!  Search through the ASCII file open on unit=io for a character
!  string associated with the input variable VAR.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        IO: device number for an ASCII file                     [I4]
!       OPT: task options                                        [A1]
!            0: begin search at current location (default)
!            1: rewind unit=io before beginning search
!       VAR: variable name                                       [A*]
!
!     OUTPUT:
!         FN: filename                                          [A*]
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: unit=io not ready for file i/o
!             -2: file rewind failed
!             -3: variable name not found before end of file.
!             -4: filename not found
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 30-MAR-1999
!     REVISIONS HISTORY:
!       This routine is very similar to the VI_GETFN routine.  The
!       only difference is that the character string associated with
!       VAR can include blank spaces and *special* characters such
!       as single (') or double quotes (").
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER IO,NLP,IFLG,J1,J2,NXT
CHARACTER(LEN=*) VAR,STR
CHARACTER(LEN=1) opt
CHARACTER(LEN=MXLL) TXT
LOGICAL again
INTENT(in) :: io,opt,var
INTENT(out) :: str,iflg

!Locate variable name in unit=io.
CALL STR_FILL(' ',str)
CALL FIO_LOCVAR(io,opt,var, nlp,txt,iflg)
IF (iflg.NE.0) THEN
  IFLG=-1
  RETURN
ENDIF

!LOCATE POSITION OF KEYWORD IN STR
CALL STR_FAM(var,txt, j1,j2,iflg) !case insensitive
IF (iflg.NE.0) THEN
  IFLG=-3
  RETURN
ENDIF

!GET TEXT ASSOCIATED WITH VAR
nxt=j2+1
again=(STR_SCAN(txt(nxt:nxt),' =').GT.0)
DO WHILE (again)
  nxt=nxt+1
  again=((STR_SCAN(txt(nxt:nxt),' =').GT.0).AND.(nxt.LT.LEN(txt)))
ENDDO
CALL STR_ADD(txt(nxt:),STR)

RETURN
END

SUBROUTINE VI_GETVAR(io,opt,var, ndat,iflg)
!  Retrieve the data associated with variable VAR from the
!  ASCII file open on unit=io.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!        IO: device number for an ASCII file                     [I4]
!       OPT: task options                                        [A1]
!            0: begin search at current location (default)
!            1: rewind unit=io before beginning search
!       VAR: variable name                                       [A*]
!
!     OUTPUT:
!       NDAT: number of data points found for VAR               [I4]
!       IFLG: status flag                                       [I4]
!              1: WARNING: abnormal data block termination
!              0: task successful
!             -1: unit=io not ready for file i/o
!             -2: task failed during rewind operation
!             -3: variable name not found before end of file.
!             -4: no data found before end of file.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) IO,NLP,NDAT,IFLG,JFLG,J1,J2,NXT,N
CHARACTER(LEN=*) VAR
CHARACTER(LEN=1) OPT
CHARACTER(LEN=50) NUM
CHARACTER(LEN=MXLL) STR
CHARACTER(LEN=*), PARAMETER :: LST=' `~@()={}[]|:;''",<>?'
LOGICAL OK,AGAIN
REAL(R8K) R8N
INTENT(in) :: io,opt,var
INTENT(out) :: ndat,iflg

INQUIRE(unit=io,OPENED=ok)
IF (.not.OK) THEN
  IFLG=-1
  RETURN
ELSEIF (OPT(1:1).EQ.'1') THEN
  !Rewind unit=io.
  REWIND(unit=io,iostat=iflg)
  IF (iflg.NE.0) THEN
    IFLG=-2
    RETURN
  ENDIF
ENDIF

NDAT=0
!CASE INSENSITIVE SEARCH THROUGH UNIT=IO FOR VARIABLE.
CALL FIO_LOCVAR(io,'0',var, nlp,str,iflg)
IF (iflg.NE.0) THEN
  IFLG=-3
  RETURN
ENDIF

!LOCATE POSITION OF KEYWORD IN STR
CALL STR_FAM(var,str, j1,j2,iflg) !case insensitive
IF (iflg.NE.0) THEN
  IFLG=-3
  RETURN
ENDIF

!GET DATA ASSOCIATED WITH VARIABLE
DON=yes  !Turn on data input block
nxt=j2+1
again=(nxt.GT.1)
DO WHILE (again)
  CALL TOKEN(0,lst,str,nxt, num,iflg)  !case insensitive
  IF (IFLG.EQ.0) THEN
    !PROCESS ALPHANUMERIC SEQUENCE
    CALL VI_CHKSTR(num, n,jflg)
    IF ((jflg.EQ.0).AND.DON) THEN
      !attempt to convert character sequence to a number
      CALL STR_GETR8N(num, r8n,jflg)
      IF (jflg.EQ.0) THEN
        !SAVE DATA TO WORKSPACE
        NDAT=NDAT+1
        CALL VI_SAVDAT(0,r8n, iflg)  !save data
      ELSE
        !Unexpected termination of data entry.
        IFLG=1
        AGAIN=NO
      ENDIF
    ELSEIF (JFLG.GT.0) THEN
      !temporarily turn of data entry because encountered a comment
      DON=NO
    ELSEIF (jflg.LT.0) THEN
      !Normal termination of data entry
      iflg=0
      AGAIN=NO
     ENDIF
  ELSE
    !READ next line of text
    READ(unit=io,FMT='(A)',iostat=iflg) STR
    CALL VI_CHKSTR(str, n,iflg)
    IF (iflg.NE.0) THEN
      DON=NO       !turn off data input
      IF (IFLG.LT.0) THEN
        AGAIN=NO   !data input block terminated
        IFLG=1     !unexpected termination of data entry
      ENDIF
    ELSE
      DON=YES  !Data input on
    ENDIF
    NXT=1
  ENDIF

ENDDO  !END OF DATA INPUT FOR THIS VARIABLE

IF ((ndat.LT.1).AND.(iflg.GE.0)) THEN
  !no data found for this variable
  IFLG=-4
ENDIF

RETURN
END

SUBROUTINE VI_R8NDAT(n, r8n,iflg)
!  Retrieve the nth data element from active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!         N: VI_DAT array index (positive integer)              [I4]
!       R8N: double precision number                            [R8]
!
!     OUTPUT:
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: ERROR: n less than 1
!             -2: ERROR: n too big
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:

USE FIOPAK_WORKSPACE
USE MATH_USEFUL_DP

INTEGER(I4K) N,IFLG
REAL(R8K) R8N
INTENT(in) :: n
INTENT(out) :: r8n,iflg

R8N=zero
IF (N.LT.1) THEN
  IFLG=-1
  RETURN
ELSEIF (N.GT.NDES) THEN
  IFLG=-2
  RETURN
ELSE
  R8N = VI_DAT(n)
  IFLG=0
ENDIF

RETURN
END

SUBROUTINE VI_R4NDAT(n, r4n,iflg)
!  Retrieve the nth data element from active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!         N: VI_DAT array index (positive integer)              [I4]
!       R4N: single precision number                            [R4]
!
!     OUTPUT:
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: ERROR: n less than 1
!             -2: ERROR: n too big
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:

USE FIOPAK_WORKSPACE
USE MATH_USEFUL_DP
INTEGER(I4K) N,IFLG
REAL(R4K) R4N
INTENT(in) :: n
INTENT(out) :: r4n,iflg

R4N=zero
IF (N.LT.1) THEN
  IFLG=-1
  RETURN
ELSEIF (N.GT.NDES) THEN
  IFLG=-2
  RETURN
ELSE
  R4N = VI_DAT(n)
  IFLG=0
ENDIF

RETURN
END

SUBROUTINE VI_I4NDAT(n, i4n,iflg)
!  Retrieve the nth data element from active memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!         N: VI_DAT array index (positive integer)              [I4]
!       I4N: single precision number                            [I4]
!
!     OUTPUT:
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: ERROR: n less than 1
!             -2: ERROR: n too big
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:

USE FIOPAK_WORKSPACE
USE MATH_USEFUL_DP
INTEGER(I4K) N,IFLG
INTEGER(I4K) I4N
INTENT(in) :: n
INTENT(out) :: i4n,iflg

I4N=zero
IF (N.LT.1) THEN
  IFLG=-1
  RETURN
ELSEIF (N.GT.NDES) THEN
  IFLG=-2
  RETURN
ELSE
  I4N = VI_DAT(n)
  IFLG=0
ENDIF

RETURN
END

SUBROUTINE VI_SAVDAT(iopt,r8n, iflg)
!  Save data associated with the active variable to FIOPAK
!  module-level memory.
!
!      LIBRARY: [FIOPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     INPUT:
!      IOPT: task options                                        [I4]
!           -1: erase current contents of storage array
!            0: add r8n to existing data
!       R8N: double precision number                             [R8]
!
!     OUTPUT:
!       IFLG: status flag                                       [I4]
!              0: task successful
!             -1: data storage task failed
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 26-MAR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE STRPAK
USE FIOPAK_WORKSPACE
INTEGER(I4K) IOPT,IFLG
REAL(R8K) R8N
INTENT(in) :: iopt,r8n
INTENT(out) :: iflg

IF (NDES+1.GT.NDMX) THEN
  CALL VI_ALLOC(NDES, IFLG)
  IF (iflg.NE.0) THEN
    IFLG=-1
  ENDIF
ENDIF

IF (iopt.LT.0) THEN
  NDES=0
  VI_DAT=ZERO
ENDIF

NDES=NDES+1
VI_DAT(ndes)=r8n

RETURN
END


SUBROUTINE TP_TOD1(tod)
!  Return the time of day (TOD) according to the TimPak clock
!  settings.  TOD is returned as a negative number if the
!  TimPak clock has not been set.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!         TOD: time of day (h)                           [R8]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
REAL(R8K) TOD
INTENT(out) :: TOD

 IF (TP_time(1).LT.0) THEN
   !Initialize to current time of day
   CALL TP_SETSYS
 ENDIF

 TOD =   TP_time(1)                                    &    !hours
        + (TP_time(2)                                  &    !minutes
        +  (TP_time(3)                                 &    !seconds
        +   0.001D+00*TP_time(4))/60.0D+00)/60.0D+00        !milliseconds

RETURN
END

SUBROUTINE TP_TOD2(tod)
!  Return the time of day (TOD) according to the TimPak clock
!  settings.  TOD is returned as a negative number if the
!  TimPak clock has not been set.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!         TOD: time of day (h)                           [R4]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
REAL(R4K) TOD
REAL(R8K) R8N
INTENT(out) :: TOD

  CALL TP_TOD1(r8n)
  TOD = r8n

RETURN
END

SUBROUTINE TP_TOD3(hour,minute,sec,msec)
!  Return the hour, minutes, seconds, and milliseconds
!  of the day according to the TimPak clock settings.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!       hour: hour of the day (0-23)                     [I4]
!     minute: minutes of the hour (0-59)                 [I4]
!        sec: seconds of the minute (0-59)               [I4]
!       msec: milliseconds of the second (0-999)         [I4]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
INTEGER(I4K) hour,minute,sec,msec
INTENT(out) :: hour,minute,sec,msec

  IF (TP_time(1).LT.0) THEN
    !Initialize to current time of day
    CALL TP_SETSYS
  ENDIF

  HOUR=TP_TIME(1)
  MINUTE=TP_TIME(2)
  SEC=TP_TIME(3)
  MSEC=TP_TIME(4)

RETURN
END

SUBROUTINE TP_TOD4(hour,minute,sec)
!  Return the hour, minutes, seconds of the day according
!  to the TimPak clock settings.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!       hour: hour of the day (0-23)                     [I4]
!     minute: minutes of the hour (0-59)                 [I4]
!        sec: seconds of the minute                      [R8]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
INTEGER(I4K) hour,minute
REAL(R8K) sec
INTENT(out) :: hour,minute,sec

  IF (TP_time(1).LT.0) THEN
    !Initialize to current time of day
    CALL TP_SETSYS
  ENDIF

  HOUR=TP_TIME(1)
  MINUTE=TP_TIME(2)
  SEC=TP_TIME(3) + TP_TIME(4)*0.001D+00

RETURN
END

SUBROUTINE TP_TOD5(hour,minute,sec)
!  Return the hour, minutes, seconds of the day according
!  to the TimPak clock settings.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!       hour: hour of the day (0-23)                     [I4]
!     minute: minutes of the hour (0-59)                 [I4]
!        sec: seconds of the minute                      [R4]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
INTEGER(I4K) hour,minute
REAL(R4K) sec
INTENT(out) :: hour,minute,sec

  IF (TP_time(1).LT.0) THEN
    !Initialize to current time of day
    CALL TP_SETSYS
  ENDIF

  HOUR=TP_TIME(1)
  MINUTE=TP_TIME(2)
  SEC=TP_TIME(3) + TP_TIME(4)*0.001D+00

RETURN
END

SUBROUTINE TP_TOD6(hour,minute,sec)
!  Return the hour, minutes, seconds of the day according
!  to the TimPak clock settings.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!       hour: hour of the day (0-23)                     [I4]
!     minute: minutes of the hour (0-59)                 [I4]
!        sec: seconds of the minute                      [I4]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
INTEGER(I4K) hour,minute,sec
REAL(R8K) msec
INTENT(out) :: hour,minute,sec

  IF (TP_time(1).LT.0) THEN
    !Initialize to current time of day
    CALL TP_SETSYS
  ENDIF

  HOUR=TP_TIME(1)
  MINUTE=TP_TIME(2)
  MSEC=TP_TIME(4)*0.001D+00 + HALF
  SEC=TP_TIME(3) + MSEC

RETURN
END

SUBROUTINE TP_AddDat(stamp)
!  Format and append the date in active (TimPak) memory to
!  stamp.  Dates are formated, for example, as 01-JAN-2000.
!  Use the TP_SETDAT or TP_SETSYS routines to set the date
!  prior to calling TP_AddDat.  If neither the TP_SETDAT or
!  TP_SETSYS routine is explicitly called prior to using the
!  TP_AddDat routine for the first time in a program, then
!  the TP_SETSYS is automatically invoked.
!
!     INPUT/OUTPUT
!        STAMP: character variable                     [A*]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE STRPAK
USE TIMPAK_WorkSpace
CHARACTER(LEN=*) STAMP
CHARACTER(LEN=3) :: MONTH(12)= (/'JAN','FEB','MAR','APR','MAY','JUN',  &
                                 'JUL','AUG','SEP','OCT','NOV','DEC' /)
INTENT(inout)  :: STAMP

IF (TP_DATE(1).LT.0) THEN
  !Set to current system date/time
  CALL TP_SETSYS
ENDIF

CALL STR_ADD(TP_DATSTR(7:8)//'-',STAMP)
CALL STR_ADD(MONTH(TP_DATE(2))//'-',STAMP)
CALL STR_ADD(TP_DATSTR(1:4),STAMP)

RETURN
END

SUBROUTINE TP_AddTim(stamp)
!  Format and append the time in active (TimPak) memory to
!  stamp. Use the TP_SETTIM or TP_SETSYS routines to set the
!  time prior to calling TP_AddTim.  If neither the TP_SETDAT or
!  TP_SETSYS routine is explicitly called prior to using the
!  TP_AddDat routine for the first time in a program, then
!  the TP_SETSYS is automatically invoked.
!
!     INPUT/OUTPUT
!        STAMP: character variable                     [A*]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE STRPAK
USE TIMPAK_WorkSpace
CHARACTER(LEN=*) STAMP
CHARACTER(LEN=2) HOUR,AMPM
INTEGER JFLG
INTENT(inout)  :: STAMP

IF (TP_TIME(1).LT.0) THEN
  !Set to current system date/time
  CALL TP_SETSYS
ENDIF

!DEFAULT FORMAT HH:MM:SS AM or PM
IF (TP_TIME(1).LT.13) THEN
  WRITE(HOUR,fmt='(I2)',iostat=jflg) TP_TIME(1)
  AMPM='AM'
ELSE
  WRITE(HOUR,fmt='(I2)',iostat=jflg) TP_TIME(1)-12
  AMPM='PM'
ENDIF
CALL CHR_SAR(' ','0', HOUR)
CALL STR_ADD(HOUR//':',STAMP)                  !hour of day
CALL STR_ADD(TP_TIMSTR(3:4)//':',STAMP)        !minutes
CALL STR_ADD(TP_TIMSTR(5:6)//' '//AMPM,STAMP)  !seconds

RETURN
END

SUBROUTINE TP_CPUTIM(tim)
!  Returns a processor-dependent time index that can be used to compute
!  the amount of CPU time (s) used between subsequent calls.  This routine
!  is a stub for the Fortran 95 intrinsic CPU_TIME.  Sections of a
!  program can be timed as follows.
!
!       CALL TP_CPUTIM(t1)
!         ... do some stuff [optional]
!       CALL TP_CPUTIM(t2)
!       WRITE(*,*) 'Elapsed time (s): ',t2-t1
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!        TIM: elapsed CPU time (s)                   [R8]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: If the Fortran '95 intrinsic function CPU_TIME
!               is not available with your computer system, replace
!                  CALL CPU_TIME(tim)
!               with
!                  CALL TP_ETIME(tim)
!
!               The TP_ETIME routine is a Fortran 90 based TimPak
!               routine that provides equivalent functionality.
USE TIMPAK
REAL(R8K) tim
INTENT(out)  :: tim

  !CALL CPU_TIME(tim)  !Fortran 95 intrinsic
  CALL TP_ETIME(tim) !Fortran 90 based TimPak routine

RETURN
END

SUBROUTINE TP_ETIME(tim)
!  Returns the number of seconds that have elapsed since
!  midnight on the current date (or the date when TP_ETIME
!  was first called in a program).  Sections of a program
!  can be timed as follows.
!
!       CALL TP_ETIME(t1)
!         ... do some stuff [optional]
!       CALL TP_ETIME(t2)
!       WRITE(*,*) 'Elapsed time (s): ',t2-t1
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!        TIM: elapsed time (s)                             [R8]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: The TP_ETIME routine is only provided in
!               in double precision (KIND=8) to minimize
!               round-off errors associated with adding
!               and subtracting large elapsed times.
!
USE TIMPAK_WorkSpace
REAL(R8K) tim,jdn,Tod
REAL(R8K) :: jdr = -1.0D+00
INTEGER(I4K) DATE_TIME(8),iflg
CHARACTER(LEN=12) REAL_CLOCK(3)
INTENT(out)  :: tim

  !Fortran 90 instrinsic
  CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3), Date_Time)

  !Get julian day number
  CALL TP_JDN(Date_Time(1),Date_Time(2),Date_Time(3), jdn,iflg)
  IF (jdr.LT.ZERO) THEN
    !initialize jdr to current julian day number
    jdr=jdn
  ENDIF

  !Compute time of day
  TOD =   DATE_Time(5)                                   &    !hours
        + (DATE_Time(6)                                  &    !minutes
        +  (DATE_Time(7)                                 &    !seconds
        +   0.001D+00*DATE_Time(8))/60.0D+00)/60.0D+00        !milliseconds
  TIM = ((jdn-jdr)*24.0D+00 + TOD)*3600.0D+00

RETURN
END

SUBROUTINE TP_GETDAT(year,month,day)
!  Return the day, month, and year in active TimPak memory.
!  If the date has not been explicitly set using the TP_SETDAT
!  routine prior to calling TP_GETDAT, the date is returned as
!  the current system date.
!
!     INPUT:
!        N/A
!
!     OUTPUT:
!        day: day of month (1-31)                        [I4]
!      month: month of the year (1-12)                   [I4]
!       year: 4-digit year                               [I4]
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE TIMPAK_WorkSpace
INTEGER(I4K) Day,Month,Year
INTENT(out)  :: Day, Month, Year

  IF (TP_DATE(1).LT.0) THEN
    !Set to current system date/time
    CALL TP_SETSYS
  ENDIF

  Day=TP_DATE(1)
  Month=TP_DATE(2)
  Year=TP_DATE(3)

RETURN
END

LOGICAL FUNCTION TP_IsLeapYear(year)
!  Return TP_IsLeapYear as true if YEAR is a number that is
!  exactly divisible by 4, except for century years which
!  must also be divisible by 400.
!
!     INPUT:
!         year: 4-digit number                  [I4]
!
!     OUTPUT:
!        TP_IsLeapYear: logical value           [L}
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE TIMPAK_WorkSpace
INTEGER(I4K) year
LOGICAL Leap
INTENT(in) :: year

IF (MOD(year,4).EQ.0) THEN
  !Possible leap year
  IF (MOD(year,100).EQ.0) THEN
    !Test for special case
    IF (MOD(year,400).EQ.0) THEN
      Leap=yes
    ELSE
      Leap=no
    ENDIF
  ELSE
    Leap=yes
  ENDIF
ELSE
  Leap=no
ENDIF
TP_IsLeapYear=Leap

RETURN
END

SUBROUTINE TP_JDN(year,month,day, jdn,iflg)
!  Return a julian day number for a specified Gregorian
!  calendar date.  This julian day number is the number
!  of days that has elapsed since noon on January 1 of
!  the year -4712 (in the Julian proleptic calendar proposed
!  by Joseph Scalizer in 1582).  This implementation is based
!  on version 1.003 (5-Aug-1998) of the algorithm described
!  by Peter Baum on the website
!
!        http://www.capecod.net/~pbaum/
!
!     INPUT:
!         year: 4-digit number                  [I4]
!        month: month of the year (1-12)        [I4]
!          day: day of month (1-31)             [I4]
!
!     OUTPUT:
!          jdn: julian day number               [R8]
!         iflg: status flag                     [I4]
!               0: no problems encountered
!              -1: invalid month
!              -2: invalid day
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!       20-APR-1999 RD STEWART
!         Replaced AINT function with FLOOR function to correct
!         bug in initial version of TP_JDN
!
!     COMMENTS:
!       The Julian calendar is the calendar established by Julius
!       Caesar in 46 B.C.  In this calendar, the year is fixed at
!       365 days except for every 4th year when it has 366 days
!       (i.e., a leap year).  Alternatively, the year is considered
!       to have 365.25 days but fractions are ignored until they
!       accumulate to make a full day.
!
USE TIMPAK_WorkSpace
INTEGER :: TLU(12)=(/306,337,0,31,61,92,122,153,184,214,245,275/)
!TLU = AINT ((153*M-457)/5) where M = the month (1-12)
INTEGER(I4K) year,month,day, iflg
INTEGER yyyy
REAL(R8K) jdn
INTENT(in) year,month,day
INTENT(out) jdn,iflg

IFLG=0
IF ((month.LT.1).OR.(month.GT.12)) THEN
  !Trap bad months
  IFLG=-1
  RETURN
ELSEIF ((day.LT.1).OR.(day.GT.31)) THEN
  !Trap bad days
  IFLG=-2
  RETURN
ENDIF
yyyy=year
IF (month.LT.3) THEN
  yyyy=yyyy-1
ENDIF
jdn = Day + TLU(month)  + 365.0D+00*yyyy            &
                        + FLOOR(0.25D+00*yyyy)      &
                        - FLOOR(0.01D+00*yyyy)      &
                        + FLOOR(0.0025D+00*yyyy)    &
                        + 1721118.5D+00

RETURN
END

SUBROUTINE TP_SETDAT(year,month,day, iflg)
!  Set the TimPak day, month, and year.
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     INPUT:
!        day: day of month (1-31)                        [I4]
!      month: month of the year (1-12)                   [I4]
!       year: 4-digit year (> 0)                         [I4]
!
!     OUTPUT:
!         IFLG: Task Status Flag                         [I4]
!               1: WARNING. Year less than 1st year
!                  of gregorian calendar.
!               0: no problemo
!              -1: day is outside allowed range (1-31)
!              -2: month is outside allowed range (1-12)
!              -3: month is outside allowed range (1-9999)
!              -4: unable to updated TP_DATSTR
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
USE STRPAK
INTEGER(I4K) Day,Month,Year
INTEGER(I4K) iflg,jflg
INTENT(in)  :: Day, Month, Year
INTENT(out) :: iflg

  IF ((day.LT.1).OR.(day.GT.31)) THEN
    IFLG=-1
    RETURN
  ELSEIF ((month.LT.1).OR.(month.GT.12)) THEN
    IFLG=-2
    RETURN
  ELSEIF ((year.LT.1).OR.(year.GT.9999)) THEN
    IFLG=-3
    RETURN
  ELSEIF (year.LT.1582) THEN
    !Warning: year less than first year of gegorian calendar
    IFLG=1
  ELSE
    IFLG=0
  ENDIF

  TP_DATE(1)=day
  TP_DATE(2)=month
  TP_DATE(3)=year

  !Update date-string
  TP_DATSTR=' '
  WRITE(TP_DATSTR,fmt='(I4,I2,I2)',iostat=jflg) YEAR,MONTH,DAY
  IF (jflg.NE.0) THEN
    iflg=-4
  ENDIF
  CALL CHR_SAR(' ','0', TP_DATSTR)  !Replace blanks with zero.

RETURN
END

SUBROUTINE TP_SETTIM(hour,minute,sec,msec, iflg)
!  Set the hour, minute, second, and milliseconds that have
!  elapsed since midnight (hour=0,minute=0,sec=0,msec=0) of
!  the current day.
!
!     INPUT:
!       hour: hour of the day (0-23)                     [I4]
!     minute: minutes of the hour (0-59)                 [I4]
!        sec: seconds of the minute (0-59)               [I4]
!       msec: milliseconds of the second (0-999)         [I4]
!
!     OUTPUT:
!         IFLG: Task Status Flag                         [I4]
!               0: no problemo
!              -1: hour is outside allowed range (0-23)
!              -2: minute is outside allowed range (0-59)
!              -3: second is outside allowed range (0-59)
!              -4: msec is outside allowed range (0-999)
!              -5: task failed while generating TP_TIMSTR
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
!
USE TIMPAK_WorkSpace
INTEGER(I4K) hour,minute,sec,msec
INTEGER(I4K) iflg,jflg
INTENT(in)  :: hour,minute,sec,msec
INTENT(out) :: iflg

  IF ((hour.LT.0).OR.(hour.GT.23)) THEN
    IFLG=-1
    RETURN
  ELSEIF ((minute.LT.0).OR.(minute.GT.59)) THEN
    IFLG=-2
    RETURN
  ELSEIF ((sec.LT.0).OR.(sec.GT.59)) THEN
    IFLG=-3
    RETURN
  ELSEIF ((msec.LT.0).OR.(msec.GT.999)) THEN
    !Warning: year less than first year of gegorian calendar
    IFLG=-4
  ELSE
    IFLG=0
  ENDIF

  TP_TIME(1)=hour
  TP_TIME(2)=minute
  TP_TIME(3)=sec
  TP_TIME(4)=msec

  !Update time-string
  TP_TIMSTR=' '
  WRITE(TP_TIMSTR,fmt='(I2,I2,I2,A,I3)',iostat=jflg) &
      HOUR,MINUTE,SEC,'.',MSEC
  IF (jflg.NE.0) THEN
    iflg=-5
  ENDIF
  CALL CHR_SAR(' ','0', TP_TIMSTR)  !Replace blanks with zero.

RETURN
END

SUBROUTINE TP_SETSYS
!  Set the TimPak clock to the date/time indicated by the
!  system clock.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@pnl.gov
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!       Uses Fortran 90 DATE_AND_TIME intrinsic.
!
USE TIMPAK_WorkSpace
INTEGER DATE_TIME(8)
CHARACTER(LEN=12) REAL_CLOCK(3)

  !Fortran 90 instrinsic
  CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3), Date_Time)
  TP_DATSTR=REAL_CLOCK(1)(:8)
  TP_TIMSTR=REAL_CLOCK(2)(:10)

  TP_DATE(1)=Date_Time(1)  !4-digit year
  TP_DATE(2)=Date_Time(2)  !Month
  TP_DATE(3)=Date_Time(3)  !Day

  TP_TIME(1)=Date_Time(5)  !hour of day
  TP_TIME(2)=Date_Time(6)  !minute of hour
  TP_TIME(3)=Date_Time(7)  !second of hour
  TP_TIME(4)=Date_Time(8)  !milliseconds

RETURN
END

SUBROUTINE PP_BLANK
!  Print a blank line to the active PrnPak output device(s).
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!        STR: character variable                          [A*]
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE PRNPAK_GLOBAL
INTEGER(i4k) jflg

  WRITE(unit=u4o,fmt='(/)',advance='NO',iostat=jflg) 
  IF (echo.AND.(u4e.NE.u4o)) THEN
    WRITE(unit=u4e,fmt='(/)',advance='NO',iostat=jflg)
  ENDIF

RETURN
END

SUBROUTINE PP_GETMAR(lmar,rmar)
!  Return the current PrnPak left and right page margins.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!        LMAR: left margin                          [I4]
!        RMAR: right margin                         [I4]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 15-JUN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PRNPAK_GLOBAL
INTEGER(i4k) lmar,rmar
INTENT(out) :: lmar,rmar

  lmar=Margin_Left
  rmar=Margin_Right

RETURN
END

SUBROUTINE PP_GETU4E(io)
!  Return as IO the active PrnPak echo-output device number.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!
!     OUTPUT:
!         IO: output device number                           [I4]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PRNPAK_GLOBAL
INTEGER(i4k) io
INTENT(out)  :: io

  io=u4e

RETURN
END

SUBROUTINE PP_GETU4O(io)
!  Return as IO the active PrnPak output device number.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!
!     OUTPUT:
!         IO: output device number                           [I4]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PRNPAK_GLOBAL
INTEGER(i4k) io
INTENT(out)  :: io

  io=u4o

RETURN
END

SUBROUTINE PP_SETMAR(lmar,rmar, iflg)
!  Set the PrnPak left and right page margins.  Printing begins in
!  in column (lmar+1) and ends in column (rmar-1). 
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!        LMAR: left margin                                [I4]
!        RMAR: right margin                               [I4]
!
!     OUTPUT:
!       IFLG: status flag                                 [I4]
!              0: no problemo
!             -1: left margin is negative
!             -2: right margin less than left margin
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 09-JUN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: By default, the PrnPak margins are set to 
!               lmar=0 and rmar=80.
!
USE PRNPAK_GLOBAL
USE STRPAK
INTEGER(i4k) lmar,rmar,iflg
INTENT(in) :: lmar,rmar
INTENT(out) :: iflg

  IF (lmar < 0) THEN
    iflg=-1
    RETURN
  ELSEIF (rmar < lmar) THEN
    iflg=-2
    RETURN
  ENDIF
  Margin_Left=lmar
  Margin_Right=rmar

  !GENERATE NEW PAGE FORMAT
  Page_Fmat='(T   '
  CALL SP_CATNUM(Margin_Left+1,'I','', Page_Fmat)
  CALL STR_ADD(',A)',Page_Fmat)

RETURN
END

SUBROUTINE PP_SETU4E(io, iflg)
!  Set PrnPak echo-output device number to unit=io.  If io < 0,
!  then PrnPak echo-output is turned off.  If io = 0, then echo-output
!  is set to the default output device (typically the screen = 6).
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         IO: integer number                              [I4]
!
!     OUTPUT:
!       IFLG: status flag                                 [I4]
!              0: no problemo
!             -1: specified output device not ready for
!                 for i/o.  Echo turned off.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 09-JUN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PRNPAK_GLOBAL
INTEGER(i4k) io,iflg,jflg
LOGICAL ok
INTENT(in)  :: io
INTENT(out) :: iflg

  iflg=0
  IF (io < 0) THEN
    Echo=no
    u4e=0
  ELSEIF (io.EQ.0) THEN
    u4e=u4s
    Echo=yes
  ELSEIF (io.EQ.u4s) THEN
    u4e=u4s
    Echo=yes
  ELSE
    INQUIRE(unit=io,opened=ok,iostat=jflg)
    IF (ok.AND.(jflg.EQ.0)) THEN
      Echo=yes
      u4e=io
    ELSE
      Echo=no
      u4e=0
      iflg=-1
    ENDIF
  ENDIF

RETURN
END

SUBROUTINE PP_SETU4O(io, iflg)
!  Set PrnPak output device number to unit=io.  If IO=0, then the
!  output device number is set to the default output device (typically
!  the screen = 6) and PrnPak auto-echo feature is turned off.  The
!  auto-echo feature can be turned back on by calling the PP_SETU4E
!  routine after calling the PP_SETU4O routine.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         IO: integer number                                      [I4]
!
!     OUTPUT:
!       IFLG: status flag                                         [I4]
!              1: WARNING: Output device number same as 
!                 auto-echo device
!              0: no problemo
!             -1: specified output device not ready for
!                 for i/o.
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 09-JUN-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PRNPAK_GLOBAL
INTEGER(i4k) io,iflg,jflg
LOGICAL ok
INTENT(in)  :: io
INTENT(out) :: iflg

  iflg=0
  IF (io < 0) THEN
    iflg=-1
    u4o=0
    RETURN
  ELSEIF (io.EQ.0) THEN
    u4o=u4s
    Echo=no
    u4e=0
  ELSEIF (io.EQ.u4s) THEN
    u4o=io
  ELSE
    INQUIRE(unit=io,opened=ok,iostat=jflg)
    IF (ok.AND.(jflg.EQ.0)) THEN
      u4o=io
    ELSE
      u4o=0
      iflg=-1
      RETURN
    ENDIF
  ENDIF
  IF (u4e.EQ.u4o) THEN
    IFLG=1
  ELSE
    IFLG=0
  ENDIF

RETURN
END

SUBROUTINE PP_TextWrap
!  Toggle the PrnPak text wrap feature on/off.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!        N/A
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE PRNPAK_GLOBAL
  
  TextWrap=(.not.TextWrap)

RETURN
END

SUBROUTINE PP_PRNSTR0(str,adv)
!  Print the contents of STR to the PrnPak output device(s) and
!  optionally advance to the next line.  Trailing blanks in
!  STR are ignored.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!        STR: character variable                          [A*]
!        ADV: line feed?                                  [L4]
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: When output is sent to the screen (unit=u4s) and the
!               right margin is set at 81, the DVF 6.0 compiler does
!               not suppress the line feed (even with the advance='NO'
!               option specified).  I am not sure if this is done by
!               *design* or is merely a flaw in the compiler.  Regardless,
!               it is very inconvenient.  For now, the best work around
!               seems to be to use a right margin of 80 columns.
!
USE PRNPAK_GLOBAL
USE STRPAK
CHARACTER(LEN=*) STR
LOGICAL(L4K) adv,again,more
INTEGER(I4K) ilen,jflg
INTEGER(I4K) j1,j2,dm,cnt
INTENT(in)  :: STR,ADV

  j1=1
  ilen=MAX(1,LENTRIM(str))       !Length of STR minus trailing blanks
  dm=Margin_Right-Margin_Left-2  !Characters per line - 2
  cnt=0
  again=(dm > 0)
  DO WHILE (again)
    j2=MIN(j1+dm,ilen)

    !LEFT JUSTIFY TEXT WHEN WRAPPING TO NEXT LINE
    cnt=cnt+1
    more=((cnt>1).AND.(str(j1:j1).EQ.' '))
    DO WHILE(more)
      j1=j1+1
      j2=j2+1
      more=((j2<ilen).AND.(str(j1:j1).EQ.' '))
    ENDDO

    !PRINT TEXT TO MAIN OUTPUT DEVICE
    WRITE(unit=u4o,fmt=Page_Fmat,advance='NO',iostat=jflg) str(j1:j2)
    again=((j2+1 < ilen).AND.TextWrap)
    IF (echo) THEN
      !Print text to echo-output device
      IF (u4e.EQ.u4o) WRITE(unit=u4o,fmt='(/)',advance='NO',iostat=jflg)
      WRITE(unit=u4e,fmt=Page_Fmat,advance='NO',iostat=jflg) str(j1:j2)
    ENDIF

    j1=j2+1
    IF (again) CALL PP_BLANK
  ENDDO
  IF (adv) CALL PP_BLANK

RETURN
END

SUBROUTINE PP_PRNSTR1(str)
!  Print the contents of STR to PrnPak output device and
!  the advance to the next line.  Trailing blanks in STR
!  are ignored.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!        STR: character variable                          [A*]
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE PRNPAK_GLOBAL
CHARACTER(LEN=*) STR
INTENT(in) :: str

  CALL PP_PRNSTR0(str,yes)

RETURN
END

SUBROUTINE PP_PRNSTR2(io,str,adv)
!  Print the contents of STR to unit=io and then optionally
!  advance to the next line.  Trailing blanks in STR are ignored.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         IO: valid output device number                  [I4]
!        STR: character variable                          [A*]
!        ADV: line feed?                                  [L4]
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS: 
!
USE PRNPAK_GLOBAL
USE STRPAK
CHARACTER(LEN=*) STR
LOGICAL(L4K) adv,again,more
INTEGER(I4K) io,ilen,jflg
INTEGER(I4K) j1,j2,dm,cnt
INTENT(in)  :: IO,STR,ADV

  IF (io.NE.u4s) THEN
    !Check output device
    INQUIRE(unit=io,opened=again)
    IF (.not.again) RETURN 
  ENDIF

  j1=1
  ilen=MAX(1,LENTRIM(str))       !Length of STR minus trailing blanks
  dm=Margin_Right-Margin_Left-2  !Characters per line - 2
  cnt=0
  again=(dm > 0)
  DO WHILE (again)
    j2=MIN(j1+dm,ilen)

    !LEFT JUSTIFY TEXT WHEN WRAPPING TO NEXT LINE
    cnt=cnt+1
    more=((cnt>1).AND.(str(j1:j1).EQ.' '))
    DO WHILE(more)
      j1=j1+1
      j2=j2+1
      more=((j2<ilen).AND.(str(j1:j1).EQ.' '))
    ENDDO

    !PRINT TEXT TO MAIN OUTPUT DEVICE
    WRITE(unit=io,fmt=Page_Fmat,advance='NO',iostat=jflg) str(j1:j2)
    again=((j2+1 < ilen).AND.TextWrap)

    j1=j2+1
    IF (again) WRITE(unit=io,fmt='(/)',advance='NO',iostat=jflg) 
  ENDDO
  IF (adv) WRITE(unit=io,fmt='(/)',advance='NO',iostat=jflg) 

RETURN
END

SUBROUTINE PP_PRNSTR3(io,str)
!  Print the contents of STR to PrnPak output device and
!  the advance to the next line.  Trailing blanks in STR
!  are ignored.
!
!      LIBRARY: [PRNPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         IO: valid output device number                  [I4]
!        STR: character variable                          [A*]
!
!     OUTPUT:
!        N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 18-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
USE PRNPAK_GLOBAL
CHARACTER(LEN=*) STR
INTEGER(I4K) io
INTENT(in) :: io,str

  CALL PP_PRNSTR2(io,str,yes)

RETURN
END

SUBROUTINE EP_Blank
!  Print a blank line to the active ErrPak output device(s).
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!         N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PrnPak
USE ERRPAK_GLOBAL
INTEGER(I4K) io

  IF (EP_ect.EQ.Diagnostic) THEN
    CALL PP_GETU4E(io)
    IF (io > 0) THEN
      CALL PP_PRNSTR(io,' ',yes)
    ELSE
      CALL PP_BLANK
    ENDIF
  ELSE
    CALL PP_BLANK
  ENDIF

RETURN
END

SUBROUTINE EP_PrnStr(str)
!  Print a character string to the standard ErrPak
!  output device(s).
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     INPUT:
!         STR: ErrPak message or text               [A*]
!
!     OUTPUT:
!         N/A
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) STR
INTEGER(I4K) io,lmar,rmar,jflg
INTENT(in) :: str

  !Save PrnPak Page Margins and then change to ErrPak default
  CALL PP_GetMar(lmar,rmar)
  IF (EP_ect.EQ.Diagnostic) THEN
    CALL PP_SetMar(2,80, jflg)
  ELSEIF (EP_ect.EQ.Warning) THEN
    CALL PP_SetMar(2,80, jflg)
  ELSE
    CALL PP_SetMar(10,80, jflg)
  ENDIF

  IF (EP_ect.EQ.Diagnostic) THEN
    CALL PP_GETU4E(io)
    IF (io > 0) THEN
      CALL PP_PRNSTR(io,str,yes)
    ELSE
      CALL PP_PRNSTR(str,yes)
    ENDIF
  ELSE
    CALL PP_PRNSTR(str,yes)
  ENDIF

  !Restore original PrnPak Page Margins
  CALL PP_SetMar(lmar,rmar, jflg)

RETURN
END


SUBROUTINE EP_NewErr0(nor,eid,ect)
!  Initiate a new ErrPak traceback sequence.
!
!     INPUT:
!        NOR: Name of routine generating error            [A*]
!        EID: Message identification number               [I4]
!        ECT: Error Code Type                             [I4]
!             >= 1: Warning
!                0: Code diagnostic
!            <= -1: Unrecoverable (fatal) Error
!
!     OUTPUT:
!        N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL

CHARACTER(LEN=*) NOR
CHARACTER(LEN=MXML) STR
INTEGER(I4K) eid,ect,io
INTEGER(I4K) lmar,rmar,jflg
INTENT(in)  :: nor,eid,ect

  IF (ActiveError) THEN
    CALL EP_EndErr0(' ',0)
  ENDIF

  !Decode and save type of error message.
  IF (ect > 0) THEN
    EP_ect = Warning
  ELSEIF (ect < 0) THEN
    EP_ect = Fatal
  ELSE
    EP_ect = Diagnostic
  ENDIF
  ActiveError=yes  !Start new error sequence

  !Save PrnPak Page Margins and then change to ErrPak default
  CALL PP_GetMar(lmar,rmar)
  CALL PP_SetMar(0,80, jflg)

  !CONSTRUCT INITIAL TRACEBACK MESSAGE
  STR=' '
  IF (ECT < 0) THEN
    !FATAL ERROR
    STR='UNRECOVERABLE ERROR IN'
    CALL STR_ADD(' '//Nor, STR)
    CALL PP_PRNSTR(STR,yes)
    STR='  TRACEBACK (reverse calling sequence)'
    CALL PP_PRNSTR(STR,yes)
    CALL EP_NewErr1(nor,eid)

  ELSEIF (ECT > 0) THEN
    !WARNING MESSAGE
    CALL STR_ADD(Nor, STR)
    CALL STR_ADD(' '//'WARNING ', STR)
    CALL SP_CATNUM(EID,'I',' (', STR)
    CALL STR_ADD(')', STR)
    CALL PP_PRNSTR(STR,yes)

  ELSE
    !CODE DIAGNOSTIC.  Diagnostic information is printed to
    !the PrnPak echo-output device unit but not the main
    !output device number.
    CALL STR_ADD(Nor, STR)
    CALL STR_ADD(' '//'DIAGNOSTIC', STR)
    CALL SP_CATNUM(EID,'I',' (', STR)
    CALL STR_ADD(')', STR)
    CALL PP_GETU4E(io)
    IF (io > 0) THEN
      CALL PP_PRNSTR(io,STR,yes)
    ELSE
      CALL PP_PRNSTR(STR,yes)
    ENDIF
  ENDIF

  !Restore original PrnPak Page Margins
  CALL PP_SetMar(lmar,rmar, jflg)

RETURN
END


SUBROUTINE EP_NewErr1(nor,eid)
!  Print an error/warning message to the default PrnPak
!  output device(s).
!
!     INPUT:
!        NOR: Name of routine generating error            [A*]
!        EID: Error code identification number            [I4]
!
!     OUTPUT:
!       IFLG: status flag                                 [I4]
!              0: no problemo
!             -1: specified output device not ready for
!                 for i/o.  Echo turned off.
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) NOR
CHARACTER(LEN=MXML) STR
INTEGER(I4K) eid,ilen
INTEGER(I4K) lmar,rmar,jflg
INTENT(in)  :: nor,eid

  IF (.not.ActiveError) THEN
    CALL EP_NewErr0(nor,eid,Diagnostic)  !Default to diagnostic
  ENDIF
  IF (LineFeed) THEN
    CALL EP_BLANK
    LineFeed=no
  ENDIF

  !Save PrnPak Page Margins and then change to ErrPak default
  CALL PP_GetMar(lmar,rmar)
  IF (EP_ect.EQ.Diagnostic) THEN
    CALL PP_SetMar(2,80, jflg)
  ELSEIF (EP_ect.EQ.Warning) THEN
    CALL PP_SetMar(2,80, jflg)
  ELSE
    CALL PP_SetMar(6,80, jflg)
  ENDIF

  STR=' '
  CALL STR_ADD(Nor, STR)
  ilen=LENTRIM(STR)+2
  CALL SP_CATNUM(EID,'I','(', STR(ilen:))
  CALL STR_ADD(')', STR)
  CALL PP_PRNSTR(STR,yes)

  !Restore original PrnPak Page Margins
  CALL PP_SetMar(lmar,rmar, jflg)

RETURN
END


SUBROUTINE EP_EndErr0(nor,iopt)
!  End the active Errpak traceback sequence and optionally
!  terminate program execution.
!
!     INPUT:
!        NOR: Name of routine generating error               [A*]
!       IOPT: task options                                   [I4]
!             0: continue program execution (no message)
!            -1: terminate program execution immediately
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL

INTEGER(I4K) iopt
CHARACTER(LEN=*) NOR
CHARACTER(LEN=MXML) STR
INTEGER(I4K) lmar,rmar,jflg
INTENT(in)  :: iopt,nor

  IF (.not.ActiveError) RETURN  !Nothing to do
  IF (LineFeed) THEN
    CALL EP_Blank
    LineFeed=no
  ENDIF

  !Save PrnPak Page Margins and then change to ErrPak default
  CALL PP_GetMar(lmar,rmar)
  CALL PP_SetMar(0,80, jflg)

  ActiveError=no  !Reset error trap
  IF (iopt < 0) THEN
    !KILL PROGRAM
    STR='PROGRAM EXECUTION TERMINATED'
    IF (LENTRIM(nor)>0) THEN
      CALL STR_ADD(' IN '//Nor, STR)
    ENDIF
    CALL PP_PRNSTR(STR,yes)
    STOP
  ENDIF

  IF (EP_ect.EQ.Diagnostic) THEN
    CALL EP_Blank
  ELSE
    IF (EP_ect.EQ.Fatal) THEN
      STR='ERROR CONDITION TRAPPED'
      IF (LENTRIM(nor) > 0) THEN
        CALL STR_ADD(' IN '//Nor, STR)
      ENDIF
      CALL PP_PRNSTR(STR,yes)
    ENDIF
    CALL PP_BLANK
  ENDIF

  !Restore original PrnPak Page Margins
  CALL PP_SetMar(lmar,rmar, jflg)

RETURN
END

SUBROUTINE EP_EndErr1(nor)
!  End the active Errpak traceback sequence and continue
!  execution of the program.
!
!     INPUT:
!        NOR: Name of routine generating error         [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
CHARACTER(LEN=*) NOR
INTENT(in)  :: nor

  CALL EP_EndErr0(nor,0)  !continue program execution

RETURN
END

SUBROUTINE EP_EndErr2(iopt)
!  End the active Errpak traceback sequence and optionally
!  terminate program execution.
!
!     INPUT:
!       IOPT: task options                                [I4]
!             0: continue program execution (default)
!            -1: terminate program execution immediately
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL

INTEGER(I4K) iopt
INTENT(in)  :: iopt

  CALL EP_EndErr0(' ',iopt)

RETURN
END

SUBROUTINE EP_EndErr3
!  End the active Errpak traceback sequence and continue
!  execution of the program.
!
!     INPUT:
!         N/A
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:

  CALL EP_EndErr0(' ',0)  !continue program execution

RETURN
END


SUBROUTINE EP_AddTxt(str)
!  Print (additional) traceback information to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
INTENT(in)  :: str

  IF (.not.ActiveError) THEN
    !Default to unidentified diagnostic
    CALL EP_NewErr0('UNKNOWN',999,Diagnostic)
  ENDIF
  CALL EP_PrnStr(str)

RETURN
END

SUBROUTINE EP_AddR8N0(str,r8n,fmat)
!  Print a text message followed by a double precision number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        R8N: a number                         [R8]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
CHARACTER(LEN=MXML) BUF
REAL(R8K) r8n
INTENT(in)  :: str,r8n,fmat

  IF (.not.ActiveError) THEN
    !Default to unidentified diagnostic
    CALL EP_NewErr0('UNKNOWN',999,Diagnostic)
  ENDIF
  BUF=' '
  CALL STR_Add(str,buf)
  CALL SP_CATNUM(r8n,fmat,' ', buf)
  CALL EP_PrnStr(buf)

RETURN
END

SUBROUTINE EP_AddR4N0(str,r4n,fmat)
!  Print a text message followed by a single precision number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        R4N: a number                         [R4]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
REAL(R4K) r4n
REAL(R8K) r8n
INTENT(in)  :: str,r4n,fmat

  r8n=r4n
  CALL EP_AddR8N0(str,r8n,fmat)

RETURN
END

SUBROUTINE EP_AddI4N0(str,i4n,fmat)
!  Print a text message followed by an integer number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        I4N: number                           [I4]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
INTEGER(I4K) i4n
REAL(R8K) r8n
INTENT(in)  :: str,i4n,fmat

  r8n=i4n
  CALL EP_AddR8N0(str,r8n,fmat)

RETURN
END


SUBROUTINE EP_AddR8N1(str,r8n)
!  Print a text message followed by a real number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        R8N: a number                         [R8]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
REAL(R8K) r8n
INTENT(in)  :: str,r8n

  CALL EP_AddR8N0(str,r8n,'E15')  !Default format

RETURN
END

SUBROUTINE EP_AddR4N1(str,r4n)
!  Print a text message followed by a real number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        R4N: a number                         [R4]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
REAL(R4K) R4N
INTENT(in)  :: str,r4n

  CALL EP_AddR4N0(str,r4n,'E7')  !Default format

RETURN
END

SUBROUTINE EP_AddI4N1(str,i4n)
!  Print a text message followed by an integer number
!  to the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!        I4N: a number                         [I4]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
INTEGER(I4K) i4n
INTENT(in)  :: str,i4n

  CALL EP_AddI4N0(str,i4n,'I')  !Default format

RETURN
END


SUBROUTINE EP_AddR8A0(str,xa,n,fmat)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [R8]
!          N: logical size of array            [I4]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
CHARACTER(LEN=MXML) BUF
INTEGER(I4K) i,n,ilen
INTEGER(I4K) lmar,rmar
REAL(R8K) xa(*)
INTENT(in)  :: str,xa,n,fmat

  IF (.not.ActiveError) THEN
    !Default to unidentified diagnostic
    CALL EP_NewErr0('UNKNOWN',999,Diagnostic)
  ENDIF
  CALL PP_GetMar(lmar,rmar)
  BUF=' '
  CALL STR_Add(str,buf)
  ilen=LENTRIM(buf)
  DO i=1,n
    CALL SP_CATNUM(xa(i),fmat,' ', buf(ilen:))
    IF (LENTRIM(buf) > rmar-16) THEN
      CALL EP_PrnStr(buf)
      BUF=' '
    ENDIF
  ENDDO
  IF (LENTRIM(buf) > 0) THEN
    CALL EP_PrnStr(buf)
  ENDIF

RETURN
END

SUBROUTINE EP_AddR4A0(str,xa,n,fmat)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [R4]
!          N: logical size of array            [I4]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
CHARACTER(LEN=MXML) BUF
INTEGER(I4K) i,n,ilen
INTEGER(I4K) lmar,rmar
REAL(R4K) xa(*)
INTENT(in)  :: str,xa,n,fmat

  IF (.not.ActiveError) THEN
    !Default to unidentified diagnostic
    CALL EP_NewErr0('UNKNOWN',999,Diagnostic)
  ENDIF
  CALL PP_GetMar(lmar,rmar)
  BUF=' '
  CALL STR_Add(str,buf)
  ilen=LENTRIM(buf)
  DO i=1,n
    CALL SP_CATNUM(xa(i),fmat,' ', buf(ilen:))
    IF (LENTRIM(buf) > rmar-8) THEN
      CALL EP_PrnStr(buf)
      BUF=' '
    ENDIF
  ENDDO
  IF (LENTRIM(buf) > 0) THEN
    CALL EP_PrnStr(buf)
  ENDIF

RETURN
END

SUBROUTINE EP_AddI4A0(str,xa,n,fmat)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [I4]
!          N: logical size of array            [I4]
!       FMAT: short-hand format specifier      [A*]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE StrPak
USE PrnPak
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str,fmat
CHARACTER(LEN=MXML) BUF
INTEGER(I4K) i,n,ilen
INTEGER(I4K) lmar,rmar
INTEGER(I4K) xa(*)
INTENT(in)  :: str,xa,n,fmat

  IF (.not.ActiveError) THEN
    !Default to unidentified diagnostic
    CALL EP_NewErr0('UNKNOWN',999,Diagnostic)
  ENDIF
  CALL PP_GetMar(lmar,rmar)
  BUF=' '
  CALL STR_Add(str,buf)
  ilen=LENTRIM(buf)
  DO i=1,n
    CALL SP_CATNUM(xa(i),fmat,' ', buf(ilen:))
    IF (LENTRIM(buf) > rmar-9) THEN
      CALL EP_PrnStr(buf)
      BUF=' '
    ENDIF
  ENDDO
  IF (LENTRIM(buf) > 0) THEN
    CALL EP_PrnStr(buf)
  ENDIF

RETURN
END


SUBROUTINE EP_AddR8A1(str,xa,n)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [R8]
!          N: logical size of array            [I4]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
INTEGER(I4K) n
REAL(R8K) xa(*)
INTENT(in)  :: str,xa,n

  CALL EP_AddR8A0(str,xa,n,'E15')  !Default number format.

RETURN
END


SUBROUTINE EP_AddR4A1(str,xa,n)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [R4]
!          N: logical size of array            [I4]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
INTEGER(I4K) n
REAL(R4K) xa(*)
INTENT(in)  :: str,xa,n

  CALL EP_AddR4A0(str,xa,n,'E7')  !Default number format.

RETURN
END

SUBROUTINE EP_AddI4A1(str,xa,n)
!  Print a text message followed by an array of numbers to
!  the standard ErrPak output device(s).
!
!     INPUT:
!        STR: message                          [A*]
!         XA: array of numbers                 [I4]
!          N: logical size of array            [I4]
!
!     OUTPUT:
!         N/A
!
!      LIBRARY: [ERRPAK]
!       STATUS: [ALPHA]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 20-JUL-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
USE ERRPAK_GLOBAL
CHARACTER(LEN=*) str
INTEGER(I4K) n
INTEGER(I4K) xa(*)
INTENT(in)  :: str,xa,n

  CALL EP_AddI4A0(str,xa,n,'I')  !Default number format.

RETURN
END

