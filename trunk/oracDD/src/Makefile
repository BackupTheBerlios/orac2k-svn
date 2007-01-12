# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header$
# $Log$
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
io.h: $(shell ../bin/makedepf90 $(SRC) > io.h)
include io.h
$(shell rm io.h)
