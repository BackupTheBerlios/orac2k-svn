/*   Configuration file */

/* Effective number of molecules of the solute.  */

#define  _N_PROT_MAX_    30000

/* Sites of the solvent molecules */

#define  _SIT_SOLV_	100

/* Total sites of the solute molecules */

#define  _SIT_SOLU_    100000

/* Number of types of the solvent molecule*/

#define  _TYP_SOLV_	2

/* Number of types of the solute molecules */

#define  _TYP_SOLU_	250

/* Total Number of residues in the molecule */ 

#define  _NRES_		28700

/* Total Number of residues in the solvent molecule */ 

#define  _MOL_NRES_     2

/* Total number of groups */

#define  _TGROUP_	34000

#define  _SOLV_GROUP_	100

/* Maximum bin dimension of the G(r) */

#define  _MAX_BIN_      200

/* Maximum bin dimension for Hydrogen bond histogram */

#define  _MAX_HHISTO_BIN_     1000


/* Dimension of neighbor lists. Defines expected number of neighbors 
for respectively, solvent-solvent, solvent-solute and solute-solute 
neighbor list. Here, neighbors are counted for molecules and groups */

#define  _NAT_PP_    190   /* Dimension inter-neighbor list _NAT_PP_*_TGROUP_*/


/* EWALD summation */
/* 
*=================== Definition of the sum boundaries =================*
*                                                                      *
*     The sum over the reciprocal space is carried out in the 3        *
*     dimensions from -LMAX, 0, LMAX  -MMAX, 0, MMAX  -NMAX, 0, NMAX   *
*                                                                      *
*=================== Definition of the sum boundaries =================*
*/

#define  _LMAX_		2
#define  _MMAX_		2
#define  _NMAX_		2


/* Number of Hoover chains  */

#define _HOOVER_	50

/* Total number of residue available, i.e. how many residue types are
used.                                                             */

#define _TOT_RES_	100

/* Total number of connection per atom in the connection table */

#define _CONN_	7 

/*
=================================================================
==========  Analysis Dimensions Start here ======================
=================================================================
*/

/****************************************************************
*********** Dimension  of the bond, bending, and torsions *******
*********** arrays used by drive_analysis ***********************
*****************************************************************/

#define  _ANL_TOPO_	500

/****************************************************************
*********** Types of solvent atoms used by gofr *****************
*****************************************************************/

/* Number of types of the solvent molecule*/

#define  _TYP_SOLV_GOFR_	4

/****************************************************************
****************** Voronoi **************************************
*****************************************************************/

/*    Maximum of neighbors for atom to be considered            */

/*
#define _MAX_NEIGH_	120
*/

#define _MAX_NEIGH_	 590

/*    Maximum of bonded neighbors                               */

/*
                        40
*/
#define _MAX_IG_NNL_	 40

/*    Maximum allowed planes in the Voronoi construcion         */

/*
#define   MAXPLA	40
*/

#define   MAXPLA        80

/*    Maximum allowed verteces in the Voronoi construcion       */

/*
#define   MAXVER	20
*/

#define   MAXVER       120

#define   _CONTACT_VORONOI_	3

/*     Cavities calculation                                     */

/*
#define   _MAX_CAVITIES_       2400
#define   _MAX_CAVITIES_ATOM_  30
#define   _MAX_CAVITIES_BIN_   310
#define   _MAX_CAVITIES_NRES_  100
*/
#define   _MAX_CAVITIES_       2
#define   _MAX_CAVITIES_ATOM_  2
#define   _MAX_CAVITIES_BIN_   2
#define   _MAX_CAVITIES_NRES_  2

/****************************************************************
****************** Correlations *********************************
*****************************************************************/

/* Parameters for the calculation of the residence times */
/* Number of groups for which to count water neighbors */

/* #define  _DLST_	3500 */
#define  _DLST_	1

/* Total number of molecules counted for all the protein */

/* #define  _VIC_	8000 */
#define  _VIC_	1

/* Number of timesteps for which to accumulate correlation time */
/* Use for residence times                                      */

/* #define  _DTM_	1000 */
#define  _DTM_	1

/****************************************************************
****************** Velocity AutoCorrelations ********************
*****************************************************************/
/* Atom of the solvent */

#define _NAT_SOLV_	2

/* Atom of the solute */

#define _NAT_SOLU_	1

/****************************************************************
********************* MTS VACF parameters ***********************
****************************************************************/

#define _MTS_NUM_CORR_  2
#define _MTS_NUM_SPEC_  2
#define _MTS_SUBUN_  2


/* Number of correlation timestep to accumulate during the run */

#define _NUM_CORR_	1
/* Total number of subunits. To be used in analysis */

#define _SUBUN_	 _N_PROT_MAX_

/****************************************************************
****************** Type of boundary conditions ******************
****************** define the appropriate choice ****************
*****************************************************************/

#undef MAGIC /* Use magic number for cubic boundary conditions */

#define  ANINT /* Use ANINT for cubic boundary conditions */

/****************************************************************
****************** ABMD crystalize parameters *******************
*****************************************************************/

#define _ABMD_CRYST_    20

/****************************************************************
****************** More exotic parameters ***********************
******************* Never need to change*************************
*****************************************************************/


/* Respectively, number of fields, excluded fields and maximum 
number of atoms for each field */

#define _F1_	5
#define _F2_	1
#define _F3_	5


/****************************************************************
****************** Assign dimension to the buffer ***************
********** Total buffer is 3*_BUFFER_ATOMS_*_BUFFER_TIME_ *******
*****************************************************************/

/*
#define _BUFFER_ATOMS_   100
#define _BUFFER_TIME_    4096
#define _BUFFER_FFT_     2*_BUFFER_TIME_
*/
#define _BUFFER_ATOMS_   2
#define _BUFFER_TIME_    2
#define _BUFFER_FFT_     2*_BUFFER_TIME_

/****************************************************************
********** Dimensions of the dynamic matrix *********************
*****************************************************************/

#define _DYNAMIC_DIM_      2
#define _CHEB_ORDER_       2

/* Grid points for PME */
#define _MAX_GRID_  144
#define _FFT1_      144
#define _FFT2_      144
#define _FFT3_      144

/* B-spline order for PME */
#define _MORD_ 7

/* B-spline Dimensions for Direct sum*/

#define _DIR_SPLINE_      6000

/* Maximum number of cell for linked cell neighbor list */ 

/*
#define _CELLMAX_    125000  *  = 50 x 50 x 50 *
*/

#define _CELLMAX_    125000  /*  = 50 x 50 x 50 */  


/* Maximum number of cells surrounding the target cell in linked cell list */ 

/*
#define _INDMAX_  5000
*/
#define _INDMAX_  5000

/* Define length of a word */

#define _NBYTE_ 4

/* Define length of a real */

#define _RBYTE_ 8
