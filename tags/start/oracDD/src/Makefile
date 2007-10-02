# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header$
# $Log$
# Revision 1.1  2007/01/11 18:48:57  mmarchi
# Initial revision
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

#FC = gfortran
#F90 = gfortran
#CC = gcc
#LD = gfortran

FC = g95
F90 = g95
CC = gcc
LD = g95

#FC = f90
#F90 = f90
#CC = gcc
#LD = f90

#FC = pgf90
#F90 = pgf90
#CC = gcc
#LD = pgf90

#FC = ifort
#F90 = ifort
#CC = gcc
#LD = ifort

CFLAGS = -O
SHELL = /bin/sh
DESCRIPTION = Merge version MAX-PIE. Unit of time is fs.
STATUS = beta version 
RELEASE = 1.0
CPP = /lib/cpp -traditional
CC = cc
LIB = nolib
CURDIR = `pwd`
DIST_DIR = tmp_dir/orac/source
LIBRARIES = NONE
#FFLAGS = -g 
#FFLAGS = -g -fbounds-check -ftrace=frame 
FFLAGS = -O3
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
	SystemPrm.o LennardJones.o PrmUtilities.o BondedPrm.o oracS.o 

SRC = $(patsubst %.o, %.f90, $(OBJS)) 

all: oracS.x

oracS.x: $(OBJS) 
	$(LD) -o oracS.x $(OBJS) $(LIBS)

clean:
	rm -f *.o *.x *.mod *~

copy: 
	for j in $(SRC);\
	do cp $$j oppo/$$j;done;


# Dependencies
AtomCnt.o : AtomCnt.f90 AtomCnt__GetConnections.f90 Strings.o Node.o flib.o Tops.o SecondarySeq.o Myparse.o Errors.o IndSequence.o Constants.o IndPatch.o 
AtomCnt__GetConnections.o : AtomCnt__GetConnections.f90 
Atoms.o : Atoms.f90 AtomCnt.o Constants.o 
BondedPrm.o : BondedPrm.f90 Node.o Strings.o PrmUtilities.o SystemTpg.o Resid.o Errors.o 
Constants.o : Constants.f90 
Errors.o : Errors.f90 Types.o Constants.o 
Globals.o : Globals.f90 Parameters_Globals.o Types.o 
Grammars.o : Grammars.f90 flib.o Strings.o Tree.o Types.o Errors.o 
Hash_Tops.o : Hash_Tops.f90 Errors.o Constants.o 
IndPatch.o : IndPatch.f90 Parameters_Globals.o SecondarySeq.o Errors.o 
IndSequence.o : IndSequence.f90 Tops.o IndPatch.o SecondarySeq.o Constants.o Errors.o Parameters_Globals.o 
Inputs.o : Inputs.f90 Strings.o Types.o flib.o Errors.o 
Keyword.o : Keyword.f90 flib.o Errors.o Node.o Myparse.o Constants.o 
LennardJones.o : LennardJones.f90 Errors.o Myparse.o Strings.o Resid.o Constants.o TypesPrm.o Units.o 
Myparse.o : Myparse.f90 flib.o Constants.o 
Newresidues_.o : Newresidues_.f90 
Node.o : Node.f90 Constants.o 
Numerics.o : Numerics.f90 dminv.f 
PDB.o : PDB.f90 Node.o Strings.o Errors.o Constants.o 
Parameters.o : Parameters.f90 Resid.o flib.o Myparse.o Strings.o Errors.o Tree.o Parameters_Globals.o Constants.o SecondarySeq.o 
Parameters_Globals.o : Parameters_Globals.f90 Types.o 
PrmUtilities.o : PrmUtilities.f90 flib.o Errors.o TypesPrm.o Node.o Myparse.o AtomCnt.o Constants.o 
Process.o : Process.f90 Tree.o Errors.o Grammars.o Setups.o Parameters.o 
Resid.o : Resid.f90 flib.o Strings.o Myparse.o Errors.o Keyword.o Node.o Constants.o 
ResidueTpg.o : ResidueTpg.f90 Tops.o Node.o Myparse.o Errors.o Parameters_Globals.o flib.o Strings.o AtomCnt.o 
SecondarySeq.o : SecondarySeq.f90 Node.o flib.o Myparse.o Parameters_Globals.o Errors.o 
Setups.o : Setups.f90 Parameters_Globals.o Myparse.o Strings.o Tree.o PDB.o Units.o Numerics.o flib.o Errors.o Constants.o 
Strings.o : Strings.f90 flib.o Types.o Constants.o Errors.o 
SystemPrm.o : SystemPrm.f90 Errors.o BondedPrm.o BondedPrm.o BondedPrm.o BondedPrm.o TypesPrm.o LennardJones.o PrmUtilities.o 
SystemTpg.o : SystemTpg.f90 Tops.o Node.o Constants.o ResidueTpg.o AtomCnt.o Errors.o Parameters_Globals.o SecondarySeq.o IndSequence.o 
Tops.o : Tops.f90 Newresidues_.f90 Tops_Charmm.f90 Tops_Orac.f90 Errors.o Parameters_Globals.o flib.o Node.o Resid.o Hash_Tops.o Strings.o Myparse.o Types.o SecondarySeq.o 
Tops_Charmm.o : Tops_Charmm.f90 
Tops_Orac.o : Tops_Orac.f90 
Tree.o : Tree.f90 Errors.o Strings.o Constants.o 
Types.o : Types.f90 Constants.o 
TypesPrm.o : TypesPrm.f90 Errors.o Strings.o Myparse.o Resid.o 
TypesTpg.o : TypesTpg.f90 
Units.o : Units.f90 
flib.o : flib.f90 
oracS.o : oracS.f90 SystemPrm.o SystemTpg.o AtomCnt.o Tops.o Units.o Process.o Grammars.o Inputs.o Tree.o 
