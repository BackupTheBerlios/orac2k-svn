c*-*
c*-* Copyright (C) 2005 Massimo Marchi and Piero Procacci
c*-* Contact for info M. Marchi, CEA,  Gif Sur Yvette 91191 (FRANCE) 
c*-* Email:Massimo.Marchi.NOSPAM@NOSPAM.cea.fr
c*-* 
c*-* This program is free software; you can redistribute it and/or
c*-* modify it under the terms of the GNU General Public License
c*-* as published by the Free Software Foundation; either version 2
c*-* of the License, or (at your option) any later version.
c*-* 
c*-* This program is distributed in the hope that it will be useful,
c*-* but WITHOUT ANY WARRANTY; without even the implied warranty of
c*-* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c*-* GNU General Public License for more details.
c*-* 
c*-* You should have received a copy of the GNU General Public License
c*-* along with this program; if not, write to the Free Software
c*-* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
      PROGRAM orac

************************************************************************
*                                                                      *
*                       ORAC of BLAKES 7                               *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*                              O R A C                                 *
*                                                                      *
*                           (Beta Version)                             *
*                                                                      *
*          "A Molecular Dynamics Program to Simulate Complex Systems"  *
*                                                                      *
*              Copyright(C) 1989 - 2005 by Massimo Marchi and          *
*                          Piero Procacci                              *
*                                                                      *
*                        All Right Reserved                            *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
* ORAC is provided "as is" and without any warranty express or implied.*
* The user assumes all risks of using ORAC.                            *
*                                                                      *
* The user may make copies of ORAC for his/her own use, and modify     *
* those copies. The user MAY NOT distribute any copy of the origina    *
* or modified source code to any users at any sites other than his/her *
* own.                                                                 *
*                                                                      *
*                                                                      *
*                                                                      *
*                      Written by  Massimo Marchi                      *
*                                                                      *
*                      Centre  d'Etudes Saclay                         *
*                       Gif sur Yvette, FRANCE                         *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*======================= DECLARATIONS ==================================

#if defined PARALLEL
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'
      INCLUDE 'orac.h'
      INCLUDE 'parallel.h'

*-------------- DEFINITION OF THE SCRATCH COMMON BLOCK -----------------


*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER i,iret,ierr,TimeToGo,TRemain
      REAL*8 vfcp,tfcp,elapse

*==================== EXECUTABLE STATEMENTS ============================

#if defined PARALLEL

*=======================================================================
*----- Initialization of the Parallel interface ------------------------
*=======================================================================

      CALL P_whoami(node,nprocs,ncube,nbyte)
      CALL P_open_input(node,nprocs,ncube,nbyte)
#else
      iret=0
      errmsg='none'
      node=0
      nprocs=1
      nodex=0
      nodey=0
      nodez=0
      npy=1
      npz=1
      ncube=1
#endif

*=======================================================================
*----- Read form input files and initialiase topology arrays -----------
*=======================================================================

      CALL starta(mapnl,mapnl_slv,mapp,mapp_slv,node)

#if defined PARALLEL
*=======================================================================
*----- If pme initialize 3D grid ---------------------------------------
*=======================================================================
      
      CALL P_pme_assign_node(node,pme,nprocs,ictxt,nodex,nodey,nodez
     &     ,npy,npz)      
#endif
*=======================================================================
*----- Build the simulation box and, if needed, insert a macromolecule -
*----- in an equilibrated configuration of solvent molecules -----------
*=======================================================================

      CALL bldbox(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,iret,errmsg
     &     ,node,nprocs,ncube)

#ifdef PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif

*----------------- If iret .EQ. 1 STOP ** ------------------------------
*                                 ====

      IF(iret .EQ. 2) CALL xerror(errmsg,80,1,20)
      IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
      IF(.NOT. stoprun) THEN


*=======================================================================
*----- Verify some variables -------------------------------------------
*=======================================================================

         CALL verify_variables(iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*----- Run an MD simulation or analysis --------------------------------
*=======================================================================

         IF(analys) THEN
            CALL drive_analysis(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm
     &           ,ypcm,zpcm,node,nodex,nodey,nodez,ictxt,npy,npz,nprocs
     &           ,ncube)
         ELSE IF(mdsim) THEN
            CALL mtsmd(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm,ypcm,zpcm)
         ELSE IF(minimize .OR. frequencies) THEN
            IF(nprocs .GT. 1) THEN
               errmsg=
     & 'Minimization must be performed with no more than'/
     & /' *one* CPU.'
               CALL xerror(errmsg,80,1,20)
            END IF
            CALL run_minimize(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,gh,xpcm
     &           ,ypcm,zpcm)
         END IF
      END IF

      WRITE(kprint,100)
      WRITE(kprint,200)
      CALL close_files
#if defined PARALLEL
      CALL MPI_FINALIZE(ierr)
#endif

*================= END OF EXECUTABLE STATEMENTS ========================

100   FORMAT(/ /
     &     '     *****************************************************',
     &     '*******************'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *                        PROGRAM COMPLETED           ',
     &     '                  *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *****************************************************',
     &     '*******************'/)

200   FORMAT(/ /10x,
     &'|--------------------------------------------------------------|' 
     &,/10x,
     &"|    ``They're unfriendly.  Which is fortunate really;         |"
     &,/10x,
     &"|      they'd be difficult to love.''                          |"
     &,/10x,
     &"|                              . . . Avon    Blake's Seven     |"
     &,/10x,
     &"|                                                              |"
     &,/10x,
     &"|--------------------------------------------------------------|"
     &/ /)
      STOP
      END
