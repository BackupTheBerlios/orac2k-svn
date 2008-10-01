#/---------------------------------------------------------------------\
#   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
#                                                                      |
#    This software is a computer program named oracDD whose            |
#    purpose is to simulate and model complex molecular systems.       |
#    The code is written in fortran 95 compliant with Technical        |
#    Report TR 15581, and uses MPI-1 routines for parallel             |
#    coding.                                                           |
#                                                                      |
#    This software is governed by the CeCILL license under             |
#    French law and abiding by the rules of distribution of            |
#    free software.  You can  use, modify and/ or redistribute         |
#    the software under the terms of the CeCILL icense as              |
#    circulated by CEA, CNRS and INRIA at the following URL            |
#    "http://www.cecill.info".                                         |
#                                                                      |
#    As a counterpart to the access to the source code and rights      |
#    to copy, modify and redistribute granted by the license,          |
#    users are provided only with a limited warranty and the           |
#    software's author, the holder of the economic rights, and         |
#    the successive licensors have only limited liability.             |
#                                                                      |
#    The fact that you are presently reading this means that you       |
#    have had knowledge of the CeCILL license and that you accept      |
#    its terms.                                                        |
#                                                                      |
#    You should have received a copy of the CeCILL license along       |
#    with this program; if not, you can collect copies on the URL's    |
#    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
#    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
#                                                                      |
#----------------------------------------------------------------------/
# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header: /cvsroot/orac2k/oracDD/src/Makefile,v 1.20 2007/07/20 14:07:43 mmarchi Exp $
# $Log: Makefile,v $
# Revision 1.20  2007/07/20 14:07:43  mmarchi
# *** empty log message ***
#
# Revision 1.19  2007/02/15 16:17:17  mmarchi
# *** empty log message ***
#
# Revision 1.18  2007/01/31 17:10:51  mmarchi
# *** empty log message ***
#
# Revision 1.16  2007/01/26 06:58:17  mmarchi
# *** empty log message ***
#
# Revision 1.15  2007/01/25 15:21:59  mmarchi
# *** empty log message ***
#
# Revision 1.14  2007/01/24 16:51:53  mmarchi
# *** empty log message ***
#
# Revision 1.13  2007/01/23 18:07:40  mmarchi
#
# Added License statement in the header of each module or routine
#
# Revision 1.12  2007/01/23 16:43:47  mmarchi
# Read from topology/parameter binary file!
#
# Revision 1.11  2007/01/22 14:17:15  mmarchi
# *** empty log message ***
#
# Revision 1.10  2007/01/18 16:44:49  mmarchi
# *** empty log message ***
#
# Revision 1.9  2007/01/17 16:57:32  mmarchi
# *** empty log message ***
#
# Revision 1.8  2007/01/17 09:50:42  mmarchi
# *** empty log message ***
#
# Revision 1.7  2007/01/16 15:06:01  mmarchi
# *** empty log message ***
#
# Revision 1.6  2007/01/14 16:39:05  mmarchi
# *** empty log message ***
#
# Revision 1.5  2007/01/13 21:25:30  mmarchi
# *** empty log message ***
#
# Revision 1.4  2007/01/12 19:41:26  mmarchi
# *** empty log message ***
#
# Revision 1.3  2007/01/12 09:59:49  mmarchi
# *** empty log message ***
#
# Revision 1.2  2007/01/12 09:57:52  mmarchi
# *** empty log message ***
#
# Revision 1.1.1.1  2007/01/11 18:48:57  mmarchi
# Imported sources
#
# Revision 1.6  2006/12/20 21:11:57  mmarchi
# *** empty log message ***
#
# Revision 1.5  2006/12/20 08:53:19  mmarchi
# *** empty log message ***
#
# Revision 1.4  2006/12/19 13:42:51  mmarchi
# *** empty log message ***
#
# Revision 1.3  2006/12/18 22:23:44  mmarchi
# *** empty log message ***
#
# Revision 1.2  2006/12/18 20:49:30  mmarchi
# *** empty log message ***
#
# Revision 1.1.1.1  2006/12/16 20:50:03  mmarchi
# initial import into CVS
#


include Make.what

#==========================================================================
# Define default rules
#==========================================================================

.f.o :
	$(FC) $(FFLAGS)  -c $<

.SUFFIXES: .o .f90
.f90.o :
	$(F90) $(FFLAGS) -c $<

.f90.mod :
	$(F90) $(FFLAGS) -c $<



OBJS = Constants.o Types.o Parameters_Globals.o Errors.o flib.o Node.o Strings.o\
	Tree.o Inputs.o Grammars.o Myparse.o Keyword.o Resid.o Parameters.o \
	Numerics.o Units.o PDB.o Setups.o Process.o Hash_Tops.o Tops.o SecondarySeq.o\
	IndPatch.o AtomCnt.o IndSequence.o ResidueTpg.o  SystemTpg.o TypesPrm.o \
	SystemPrm.o LennardJones.o PrmUtilities.o BondedPrm.o SimulationBox.o\
	ReadStore.o AddHydrogens_.o AtomBox.o Neighbors_S.o Neighbors.o Utilities.o Groups.o \
	Banner.o Atom.o Inout.o MDRun.o FactorizeNo.o PI_Decompose.o PI_.o PI_Neighbors.o \
	PI_Communicate.o NeighCells.o Potential.o Ewald.o PME.o PME_Routines.o Parallel.o\
	Print_Defs.o Rfft3d.o IndBox.o Forces.o Direct.o oracDD.o 

OBJS_f = AddH.o duni.o

SRC = $(patsubst %.o, %.f90, $(OBJS)) $(patsubst %.o, %.f, $(OBJS_f)) 

all: oracDD.x

oracDD.x: $(OBJS) $(OBJS_f)
	$(LD) -o oracDD.x $(OBJS) $(OBJS_f) -L$(FFTW_LIB) -L$(LIBS) $(LLFLAGS)\
	$(LOCFLAGS)

clean:
	rm -f *.o *.x *.mod *~

copy: 
	for j in $(SRC);\
	do cp $$j oppo/$$j;done;


# Dependencies
io.h: $(shell ../bin/makedepf90 $(SRC) > io.h)
include io.h
$(shell rm io.h)
