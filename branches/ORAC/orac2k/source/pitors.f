      SUBROUTINE pitors(itorsl,ntors,o1,type,ttorsp,ptorsp,ntorsp,o2,
     x                 ptorsx,p2,bug,nbug)

************************************************************************
*                                                                      *
*     Check improper torsions with potential and generate parameter    *
*     list for improper torsions.                                      *
*                                                                      *
*     ARGUMENTS:                                                       *
*                                                                      *
*     ITORSL  : List with all p-torsions found from sr. SITORS.        *
*               >> integer ITORSL(4,N1) <<                             *
*               ITORSL(1..4,I) contain the numbers of atom 1-2-3-4     *
*               in torsion no. I, I = 1,..,NTORS.             (INPUT)  *
*                                                                      *
*     NTORS   : Number of torsions.                           (INPUT)  *
*               >> integer NTORS <<                                    *
*     N1      : Physical column dim. ITORSL/PTORSX.           (INPUT)  *
*               >> integer N1 <<                                       *
*     TYPE    : Types for each atom in the torsions list.     (INPUT)  *
*               >> character*7 TYPE(*) <<                              *
*     TTORSP  : List of torsions in the potential model.      (INPUT)  *
*               >> character*7 TTORSP(4,N2) <<                         *
*               TTORSP(1..4,I) contains the types of atom 1-2-3-4      *
*               in potential no. I, I = 1,..,NTORSP.                   *
*     PTORSP  : List with the potential parameters            (INPUT)  *
*               >> real PTORSP(N2,M2) <<                               *
*     NTORSP  : Number of potential parameters                (INPUT)  *
*               >> integer NTORSP <<                                   *
*     N2      : Physical row dimension of PTORSP and column            *
*               dimension of TTORSP                                    *
*               >> integer N2 <<                                       *
*     PTORSX  : List with the potential parameters            (OUTPUT) *
*               corresponding to the torsions in ITORSL.               *
*               >> real PTORSK(M2,N1) <<                               *
*     M2      : Physical column dimension of PTORSP/PTORSX    (INPUT)  *
*               >> integer M2 <<                                       *
*     WORK    : Work array                                    (INPUT)  *
*               >> integer WORK(4,N1) <<                               *
*               In the calling program must have the same              *
*               physical column dimension as ITORSL                    *
*     IRET    : Return code.                                  (OUTPUT) *
*               >> integer IRET <<                                     *
*               IRET = 0,..,1. For each value of IRET a                *
*               corresponding error message is assigned to             *
*               ERRMSG.                                                *
*     ERRMSG  : Error message.                                (OUTPUT) *
*               >> CHARACTER*80 ERRMSG <<                              *
*                                                                      *
*---- LAST UPDATE: 05/24/89 -------------------------------------------*
*                                                                      *
*     Modified version of a subroutine written by                      *
*     Gerald Kneller Dept 48B, IBM Kingston 1989                       *
*                                                                      *
*     EXTERNALS:                                                       *
*                                                                      *
*     - NONE                                                           *
*                                                                      *
************************************************************************

*==== DECLARATIONS: ===================================================*

      IMPLICIT none

*---- ARGUMENTS: ------------------------------------------------------*

      INTEGER      ntors,ntorsp,o1,o2,p2,iret,itorsl(4,o1),bug(*),nbug
      REAL*8       ptorsx(o1,p2),ptorsp(o2,p2)
      CHARACTER*80 errmsg
      CHARACTER*7  match1(4),type(*),ttorsp(4,o2)
      CHARACTER*1 wild

*---- LOCAL VARIABLES: ------------------------------------------------*

      INCLUDE 'parst.h'
      INTEGER na,work
      PARAMETER (na=m4)
      COMMON /rag1/ work(4,na)
      INTEGER      itors,itorsp,ntorsk,m,j,j1,j2
      INTEGER      ip,iperm(6,4)
      LOGICAL err
      DATA iperm     /1,1,1,1,1,1,
     X                2,3,2,3,4,4,
     X                3,2,4,4,2,3,
     X                4,4,3,2,3,2/
      DATA wild/'x'/

*==== EXECUTABLE STATEMENTS: ==========================================*


*---- INITIALIZATION:  ------------------------------------------------*

      DO 10 itors = 1,o1
          DO 20 j2 = 1,p2
              ptorsx(itors,j2) = 0.0D0
20        CONTINUE
10    CONTINUE

*=======================================================================

*------- Look for all the matching torsions ----------------------------

      nbug = 0
      ntorsk = 0
      DO 30 itors = 1,ntors

*------- Keep matching improper torsions -------------------------------

          DO 40 itorsp = 1,ntorsp
              DO 50 ip=1,2
                  DO 60 j=1,4
                      match1(j)=ttorsp(j,itorsp)
60                CONTINUE
                  IF( type(itorsl(iperm(ip,1),itors)).EQ.match1(1).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(4).OR.
     x                type(itorsl(iperm(ip,1),itors)).EQ.match1(4).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(1)) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF
                      DO 70 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
70                    CONTINUE
                      DO 80 m=1,4
                          work(m,ntorsk) = itorsl(m,itors)
80                    CONTINUE
                      GOTO 100
                  END IF
50            CONTINUE
40        CONTINUE
          DO 90 itorsp = 1,ntorsp
              DO 110 ip=1,2
                  DO 120 j=1,4
                      IF(ttorsp(j,itorsp).EQ.wild) THEN
                          match1(j)=type(itorsl(iperm(ip,j),itors))
                      ELSE
                          match1(j)=ttorsp(j,itorsp)
                      END IF
120               CONTINUE
                  IF( type(itorsl(iperm(ip,1),itors)).EQ.match1(1).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(4).OR.
     x                type(itorsl(iperm(ip,1),itors)).EQ.match1(4).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(1)) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na*3) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF
                      DO 130 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
130                   CONTINUE
                      DO 140 m=1,4
                          work(m,ntorsk) = itorsl(m,itors)
140                   CONTINUE
                      GOTO 100
                  END IF
110           CONTINUE
90        CONTINUE
          DO 150 itorsp = 1,ntorsp
              DO 160 ip=3,6
                  DO 170 j=1,4
                      match1(j)=ttorsp(j,itorsp)
170               CONTINUE
                  IF( type(itorsl(iperm(ip,1),itors)).EQ.match1(1).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(4).OR.
     x                type(itorsl(iperm(ip,1),itors)).EQ.match1(4).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(1)) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na*3) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF

                      DO 180 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
180                   CONTINUE
                      DO 190 m=1,4
                          work(m,ntorsk) = itorsl(iperm(ip,m),itors)
190                   CONTINUE
                      GOTO 100
                  END IF
160           CONTINUE
150       CONTINUE
          DO 200 itorsp = 1,ntorsp
              DO 210 ip=3,6
                  DO 220 j=1,4
                      IF(ttorsp(j,itorsp).EQ.wild) THEN
                          match1(j)=type(itorsl(iperm(ip,j),itors))
                      ELSE
                          match1(j)=ttorsp(j,itorsp)
                      END IF
220               CONTINUE
                  IF( type(itorsl(iperm(ip,1),itors)).EQ.match1(1).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(4).OR.
     x                type(itorsl(iperm(ip,1),itors)).EQ.match1(4).AND.
     x                type(itorsl(iperm(ip,2),itors)).EQ.match1(2).AND.
     x                type(itorsl(iperm(ip,3),itors)).EQ.match1(3).AND.
     x                type(itorsl(iperm(ip,4),itors)).EQ.match1(1)) THEN
                      ntorsk = ntorsk + 1
                      IF(ntorsk .GT. na*3) THEN
                          errmsg='In pitors: Physical dimensions of the'
     x                    //' work array are insufficient. '
                          iret=1
                          RETURN
                      END IF
                      DO 230 m=1,p2
                          ptorsx(ntorsk,m) = ptorsp(itorsp,m)
230                   CONTINUE
                      DO 240 m=1,4
                          work(m,ntorsk) = itorsl(iperm(ip,m),itors)
240                   CONTINUE
                      GOTO 100
                  END IF
210           CONTINUE
200       CONTINUE
          nbug=nbug+1
          bug(nbug)=itors
c---      one thousend Bugs: we had enough 
          if(nbug.gt.1000) RETURN 
100       CONTINUE
30    CONTINUE

      if(nbug.ge.1) RETURN

*---------- Change the list of p-torsions ------------------------------

      ntors=ntorsk
      DO 250 j=1,ntors
          DO 260 j1=1,4
              itorsl(j1,j)=work(j1,j)
260       CONTINUE
250   CONTINUE

*---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
