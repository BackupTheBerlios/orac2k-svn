      SUBROUTINE P_pme_assign_node(node,pme,nprocs,ictxt,nodex,nodey
     &     ,nodez,npy,npz)

************************************************************************
*   Time-stamp: <1999-10-29 16:58:12 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb  7 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ictxt,nodex,nodey,nodez,npy,npz,nprocs,node
      LOGICAL pme

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined T3E
      INTEGER*8 mpx,mpy,mpz,one
      INTEGER*8 ictxta,nodexa,nodeya,nodeza,npya,npza
#else
      INTEGER mpx,mpy,mpz,one
#endif
      CHARACTER*80 errmsg
      DATA one/1/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      npy=1
      npz=1
      nodex=0
      nodey=0
      nodez=0
#if defined _FFT_T3E_
      IF(pme .AND. nprocs .NE. 1) THEN
         IF(nprocs .EQ. 4) THEN
            npy=2
            npz=2
         ELSE IF(nprocs .EQ. 8) THEN
            npy=4
            npz=2
         ELSE IF(nprocs .EQ. 16) THEN
            npy=4
            npz=4
         ELSE IF(nprocs .EQ. 32) THEN
            npy=8
            npz=4
         ELSE IF(nprocs .EQ. 64) THEN
            npy=8
            npz=8
         ELSE
            WRITE(errmsg,'('' Cannot run parallel PME with '',i5,
     &           '' PEs. Must be power of 2. '')') nprocs
            npy=1
            npz=1
            nodex=0
            nodey=0
            nodez=0
            CALL xerror(errmsg,80,1,2)
         END IF
         npya=npy
         npza=npz
         nodexa=nodex
         nodeya=nodey
         nodeza=nodez
         ictxta=ictxt
         CALL gridinit3d(ictxta,one,npya,npza)
         CALL gridinfo3d(ictxta,mpx,mpy,mpz,nodexa,nodeya,nodeza)
         ictxt=ictxta
         nodex=nodexa
         nodey=nodeya
         nodez=nodeza
      END IF
#elif defined _GPFFT_
      IF(pme .AND. nprocs .NE. 1) THEN
         npy=1
         npz=nprocs
         nodex=0
         nodey=0
         nodez=node
      END IF
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
