      SUBROUTINE kinetic(nstart,nend,ss_index,co,nato_slv,nmol
     &     ,cnstpp_slv,ntap,cnstpp,x0,y0,z0,vpx,vpy,vpz,vpxcm,vpycm
     &     ,vpzcm,nprot,protl,mass,massp,wmtp,ucek,pucek,temp,tempt
     &     ,tempr,temppr,tcm,rcm,stresstk,cpress,isostress,vco,masspp
     &     ,ucepr,temppra)

************************************************************************
*                                                                      *
*   KINETIC computes all Kinetic energies and temperatures of both     *
*   solvent and solute, c.o.m velocities and rms of solvent molecules. *
*   KINETIC is called only when NRESPA > 1                             *
*                                                                      *
*     SOLVEN  :  Solvent flag.                                    (I)  *
*     PROTEI  :  Macromolecule flag.                              (I)  *
*     OC      :  Transform the coordinates to simulation box      (I)  *
*                frame.                                                *
*                >> real*8  OC(3,3)  <<                                *
*     CO      :  Transform the coordinates to orthogonal          (I)  *
*                crystallographic frame.                               *
*                >> real*8  OC(3,3)  <<                                *
*     XP0     :  Solute coordinates, packed by molecule.          (I)  *
*     YP0        >> real*8  XP0(*), ... <<                             *
*     ZP0                                                              *
*     VPX     :  Solute velocities , packed by molecule.          (I)  *
*     VPY        >> real*8  XP0(*), ... <<                             *
*     VPZ                                                              *
*     VCX     :  Solvent velocities, packed by sites.             (I)  *
*     VCY        >> real*8  XP0(*), ... <<                             *
*     VCZ                                                              *
*     WMSS    :  List of the solvent atomic masses, packed by     (I)  *
*                site.                                                 *
*                >> real*8  MASSA(*) <<                                *
*     WTMOL   :  Total mass of the solvent molecule.              (I)  *
*     MASS    :  List of solute atomic masses, packed by molecule.(I)  *
*                >> real*8  MASSB(*) <<                                *
*     WMTP    :  Total mass of the solute molecules.              (I)  *
*     NATOW   :  Number of atoms per solvent molecule.            (I)  *
*     NMOL    :  Number of solvent molecules.                     (I)  *
*     NTAP    :  Number of atomic sites on the solute.            (I)  *
*                                                                      *
*     XSMOVE  :  Solvent rms                                     (I/O) *
*     YSMOVE  :  >> real*8  XSMOVE(*), ... <<                    (I/O) *
*     ZSMOVE  :                                                  (I/O) *
*                                                                      *
*     UCEK    :  Kinetic energy of the solvent molecules.         (O)  *
*     PUCEK   :  Kinetic energy of the solute.                    (O)  *
*     TEMP    :  Total temperature of the system                  (O)  *
*     TEMPT   :  Translational temperature of the solvent.        (O)  *
*     TEMPR   :  Rotational + intra temperature of the solvent.   (O)  *
*     TEMPPR  :  Temperature of the solute                        (O)  *
*     TCM     :  Translational temperature of the solute.         (O)  *
*     RCM     :  Rotational temperature of the solute.            (O)  *
*     V0X     :  Solute c.o.m. velocities                         (O)  *
*     V0Y        >> real*8  XP0(*), ... <<                        (0)  *
*     V0Z                                                         (0)  *
*                                                                      *
*     NA .GE. NMOL*NATOW, NB .GE. NTAP, NC .GE. NMOL                   *
*     If this conditions are not satisfied the program aborts.         *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      USE Module_Stress
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER   nato_slv,nmol,ntap,nprot,protl(*),ss_index(*),cnstpp_slv
     &     ,cnstpp
      REAL*8    mass(*),wmtp,co(3,3),vco(3,3),ucepr,masspp(*),temppra
      REAL*8    x0(*),y0(*),z0(*),vpx(*),vpy(*),vpz(*)
      REAL*8    ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk(3,3)
      REAL*8  vpxcm(*),vpycm(*),vpzcm(*),massp(*)
      LOGICAL cpress,isostress

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,nf,nfslv,nfslt,type,iii,jjj,count,j,nstart,nend
      REAL*8  vxyz(3)
      REAL*8  tto(2),ttcm(2),tucek,mtot,iner(3,3),i_iner(3,3),omega(3)
     &     ,erot
      INCLUDE 'parst.h'
      INTEGER na,nb,offset
      PARAMETER (na=npm,nb=m1)


*==================== EXECUTABLE STATEMENTS ============================
      
      rcm=0.0D0
      DO i=1,2
         tto(i)=0.d0
      END DO
      DO i=nstart,nend
         type=ss_index(i)
         tto(type)=tto(type)+(0.50D+0)*mass(i)*(vpx(i)**2+vpy(i)**2
     &        +vpz(i)**2)
      end do   

      DO iii=1,3
         DO jjj=1,3
            stresstk(iii,jjj) = 0.0D0
         END DO
      END DO
      DO i=1,2
         ttcm(i)=0.0D0
      END DO
      IF(nstart .EQ. 1 .AND. nend .EQ. ntap) THEN
         count=0
         rcm=0.0D0
         DO i=1,nprot
            offset=protl(count+1)
            CALL comp_ang_vel(.TRUE.,protl(count+1),mass,vpx,vpy,vpz,x0
     &           ,y0,z0,omega,iner,i_iner,erot)
            rcm=rcm+erot
            count=count+1+offset
         END DO

         count=0
         DO i=1,nprot
            offset=protl(count+1)
            type=ss_index(protl(count+1+1))
            vxyz(1)=co(1,1)*vpxcm(i)+co(1,2)*vpycm(i)+co(1,3)*vpzcm(i)
            vxyz(2)=co(2,1)*vpxcm(i)+co(2,2)*vpycm(i)+co(2,3)*vpzcm(i)
            vxyz(3)=co(3,1)*vpxcm(i)+co(3,2)*vpycm(i)+co(3,3)*vpzcm(i)
            ttcm(type)=ttcm(type)+0.5D0*massp(i)*(vxyz(1)**2+vxyz(2)**2
     &           +vxyz(3)**2)
            DO iii=1,3
               DO jjj=1,3
                  stresstk(iii,jjj) = stresstk(iii,jjj) + 
     x                 massp(i)*vxyz(iii)*vxyz(jjj)
               END DO
            END DO
            count=count+offset+1
         END DO
      END IF

      nf=3*ntap-cnstpp

      tucek=(tto(1)+ttcm(1)+tto(2)+ttcm(2))*efact
      pucek=(tto(1)+ttcm(1))*efact  
      ucek=(tto(2)+ttcm(2))*efact
      nfslv=3*nato_slv*nmol-cnstpp_slv
      nfslt=nf-nfslv

      temp=2.0D0*tucek/(gascon*nf)
      IF(nfslt .NE. 0) THEN
         temppr=2.0D0*pucek/(gascon*nfslt)
      ELSE
         temppr=0.0D0
      END IF

      tcm=ttcm(1)+ttcm(2)
      tcm=2.0D0*tcm*efact/(gascon*3.0D0*DFLOAT(nprot))
      IF(nstart .EQ. 1 .AND. nend .EQ. ntap) THEN
         IF(nf-3*nprot .NE. 0) THEN
            rcm=2.0D0*rcm*efact/(gascon*DFLOAT(nf-3*nprot))
         END IF
      END IF
      tempt=tcm
      tempr=rcm

      IF(cpress) THEN
         ucepr=0.0D0
         DO i=1,3
            DO j=i,3
               ucepr=ucepr+masspp(i)*vco(i,j)*vco(i,j)
            END DO
         END DO
         ucepr=0.5D0*ucepr*efact
         IF(isostress) THEN
            temppra=2.0D0*ucepr/gascon
         ELSE IF(FixedAngles_Stress) THEN
            temppra=(2.0D0/3.0D0)*ucepr/gascon
         ELSE IF(Orthogonal_Stress) THEN
            temppra=(2.0D0/3.0D0)*ucepr/gascon
         ELSE
            temppra=(2.0D0/6.0D0)*ucepr/gascon
         END IF
      END IF

      RETURN
      END
