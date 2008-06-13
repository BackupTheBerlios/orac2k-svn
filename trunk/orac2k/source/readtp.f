      SUBROUTINE readtp(kdump,mapnl,nmapnl)

************************************************************************
*                                                                      *
*                                                                      *
*     READTP will dump all data necessary to produce a restart         *
*     to the data file.                                                *
*                                                                      *
*                                                                      *
*----- Last update 06/5/90 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS none.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

      include 'parst.h'
      include 'cpropar.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER kdump,nmapnl,mapnl(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,nsevere
      CHARACTER*80 errmsg
      CHARACTER*8 fmt
      REAL*8 dummy

*==================== EXECUTABLE STATEMENTS ============================

      REWIND kdump
      nsevere=0

      WRITE(6,1000) 

      READ(kdump,END=100) llacc,lldon,nbone,ntap,lbond
     &     ,lbend,ltors,litor,lpnbd,int14p,int13p,nspw,nunit,lphyd,ngrp
     &     ,nrigg,nrigat,nmapnl,nmol,nato_slv,adihed,slv_exist,slt_exist


*-----------------------------------------------------------------------
*--- Check dimensions --------------------------------------------------
*-----------------------------------------------------------------------

      IF(ntap .GT. m1) THEN
         nsevere=nsevere+1
         errmsg=
     &'Number of atoms'
     &//' exceed dimensions. _SIT_SOLU_ should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(ngrp .GT. m11) THEN
         nsevere=nsevere+1
         errmsg=
     &'Number of groups'
     &//' exceed dimensions. _TGROUP_ should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(nsevere .GT. 0 .AND. nsevere .LT. 99) THEN
         CALL int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE READING TPGPRM FILE'
         CALL xerror(errmsg,80,1,2)
      END IF

      IF(llacc .GT. m1) THEN
         nsevere=nsevere+1
         errmsg=
     &'lacc dimensions'
     &//' insufficient. _SIT_SOLU_ should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(lldon .GT. m1) THEN
         nsevere=nsevere+1
         errmsg=
     &'ldon dimensions'
     &//' insufficient. _SIT_SOLU_ should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(lbond .GT. m9) THEN
         nsevere=nsevere+1
         errmsg=
     &'Bonds exceed dimensions.'
     &//' M9 in parst.inc should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(lbend .GT. m2) THEN
         nsevere=nsevere+1
         errmsg=
     &'Bends exceed dimensions.'
     &//' M2 in parst.inc should be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(ltors .GT. m3) THEN
         nsevere=nsevere+1
         errmsg=
     &'Torsions exceed dimensions.'
     &//' M3 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(litor .GT. m4) THEN
         nsevere=nsevere+1
         errmsg=
     &'I-torsions exceed dimensions.'
     &//' M4 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(lpnbd .GT. m6) THEN
         nsevere=nsevere+1
         errmsg=
     &'Non-bonded pot. exceed dimensions.'
     &//' M6 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(int14p .GT. m3) THEN
         nsevere=nsevere+1
         errmsg=
     &'1-4 interacts. exceed dimensions.'
     &//' M2 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(int13p .GT. m2) THEN
         nsevere=nsevere+1
         errmsg=
     &'1-3 interacts. exceed dimensions.'
     &//' M2 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(lpnbd*(lpnbd+1)/2 .GT. m5) THEN
         nsevere=nsevere+1
         errmsg=
     &'Non-bonded pot. exceed dimensions.'
     &//' M5 in parst.inc must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(nunit .GT. n5) THEN
         nsevere=nsevere+1
         errmsg=
     &'Residues exceed dimensions.'
     &//'_TOT_RES_ must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(nmapnl .GT. m8) THEN
         nsevere=nsevere+1
         errmsg=
     &'Bonded map exceeds dimensions.'
     &//' M8 in parst.inc  must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(nsevere .GT. 0 .AND. nsevere .LT. 99) THEN
         CALL int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE READING TPGPRM FILE'
         CALL xerror(errmsg,80,1,2)
      END IF
      
      READ(kdump,END=100) ss_point(1,1),ss_point(1,2)
      READ(kdump,END=100) (ss_point(i+1,1),i=1,ss_point(1,1))
     &     ,(ss_point(i+1,2),i=1,ss_point(1,2))
      READ(kdump,END=100) (ss_index(i),i=1,ntap)
      READ(kdump,END=100) (lacc(1,i),lacc(2,i),i=1,llacc),
     x             (ldon(1,i),ldon(2,i),i=1,lldon)
      READ(kdump,END=100) (nbtype(i),i=1,ntap),(nres(i,1),i=1,ntap)
     &     ,(nres(i,2),i=1,ntap),(mback(i),i=1,ntap),(grppt(1,i),i=1
     &     ,ngrp),(grppt(2,i),i=1,ngrp)
      READ(kdump,END=100) (lbnd(1,i),i=1,lbond),(lbnd(2,i),i=1,lbond),
     &     (concta(i,1),(concta(i,j),j=2,concta(i,1)+1),i=1,ntap),
     &             ((lbndg(j,i),i=1,lbend),j=1,3),
     &             ((ltor(j,i),i=1,ltors),j=1,4),
     &             ((litr(j,i),i=1,litor),j=1,4)
      READ(kdump,END=100) ((int14(j,i),j=1,2),i=1,int14p),
     x            (type14(i),i=1,int14p),
     x             ((int13(j,i),j=1,2),i=1,int13p),
     x             ((nhtype(i,j),j=1,lpnbd),i=1,lpnbd)
      READ(kdump,END=100) wmtp,(chrge(i),i=1,ntap),(mass(i),i=1,ntap)
      READ(kdump,END=100) (mapnl(i),i=1,nmapnl)
      IF(amphi) THEN
          READ(kdump,END=100) nbun,iz,
     x                 ((potbo(i,j),i=1,lbond),j=1,2),
     x                 ((potbe(i,j),i=1,lbend),j=1,4),
     x                 ((ptorj(i,j),i=1,ltors),j=1,4),
     x                 ((potit(i,j),i=1,litor),j=1,3),
     x                 (pnbd1(i),i=1,lpnbd),(pnbd2(i),i=1,lpnbd),
     x                 (pnbd3(i),i=1,lpnbd),(pnbd4(i),i=1,lpnbd),
     &                 (sjorg(i),i=1,lpnbd),
     x                 (ejorg(i),i=1,lpnbd),(phyd1(i),i=1,lphyd),
     x                 (phyd2(i),i=1,lphyd)
          IF(iz .EQ. 1) THEN
             READ(kdump) 
     &                 (ecc6(i),i=1,lpnbd*(lpnbd+1)/2),
     &                 (ecc12(i),i=1,lpnbd*(lpnbd+1)/2)
          END IF
      ELSE
         IF(old_tpg) THEN
            READ(kdump,END=100) nbun,iz,
     &           ((potbo(i,j),i=1,lbond),j=1,2),
     &           ((potbe(i,j),i=1,lbend),j=1,4),
     &           ((dummy,i=1,lbend),j=1,3),
     &           ((potto(i,j),i=1,ltors),j=1,2),
     &           ((potit(i,j),i=1,litor),j=1,3),
     &           (pnbd1(i),i=1,lpnbd),(pnbd2(i),i=1,lpnbd),
     &           (pnbd3(i),i=1,lpnbd),(pnbd4(i),i=1,lpnbd),
     &           (phyd1(i),i=1,lphyd),
     &           (phyd2(i),i=1,lphyd)
         ELSE
            READ(kdump,END=100) nbun,iz,
     &           ((potbo(i,j),i=1,lbond),j=1,2),
     &           ((potbe(i,j),i=1,lbend),j=1,4),
     &           ((potto(i,j),i=1,ltors),j=1,2),
     &           ((potit(i,j),i=1,litor),j=1,3),
     &           (pnbd1(i),i=1,lpnbd),(pnbd2(i),i=1,lpnbd),
     &           (pnbd3(i),i=1,lpnbd),(pnbd4(i),i=1,lpnbd),
     &           (phyd1(i),i=1,lphyd),
     &           (phyd2(i),i=1,lphyd)
         END IF
            IF(iz .EQ. 1) THEN
             READ(kdump) 
     &                 (ecc6(i),i=1,lpnbd*(lpnbd+1)/2),
     &                 (ecc12(i),i=1,lpnbd*(lpnbd+1)/2)
          END IF
      END IF

      IF(nbun .GT. nores) THEN
         nsevere=nsevere+1
         errmsg=
     &'Units exceed dimensions.'
     &//' _NRES_ must be increased.'
         CALL xerror(errmsg,80,1,30)
      END IF

      IF(nsevere .GT. 0 .AND. nsevere .LT. 99) THEN
         CALL int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE READING TPGPRM FILE'
         CALL xerror(errmsg,80,1,2)
      END IF

      READ(kdump,END=100) (nasitp(i),i=1,nrigg), (nbsitp(i),i=1,nrigg)
      READ(kdump,END=100) ((asitp(j,i),j=1,nasitp(i)),i=1,nrigg)
      READ(kdump,END=100) ((bsitp(j,i),j=1,nbsitp(i)),i=1,nrigg)
      READ(kdump,END=100) (beta(i),i=1,ntap),(betb(i),i=1,ntap),
     x     (alnbd(i),i=1,lpnbd)
      READ(kdump,END=100) (prsymb(i),i=1,nunit)
      READ(kdump,END=200) (mend(i),i=1,nbun)

!===== Make the file compatible with older versions 

200   CONTINUE
      RETURN

100   CONTINUE
      WRITE(6,*)
      errmsg= ' End of file encountered. Wrong version or corrupted '// 
     x              ' topology file. Abort.'
      WRITE(6,'(a80)') errmsg

      WRITE(6,*)
1000  FORMAT(//5x,'<------ Reading Topology-Parameter File ------->'//)
       
*================= END OF EXECUTABLE STATEMENTS ========================

      STOP
      END
