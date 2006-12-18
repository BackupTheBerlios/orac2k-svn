# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header$
# $Log$
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
FFLAGS = -O
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



OBJS = Class_Tree.o Xerror_Mod.o GRAMMAR_Mod.o flib.o  STRINGS_Mod.o  NODES_Mod.o\
	ORAC_INPUT_Mod.o TYPES.o ORAC_PROCESS_Mod.o MYPARSE_Mod.o ERROR_Mod.o\
	TOPS_Mod.o PARAMETERS_Mod.o PARAMETERS_GLOBALS.o GLOBALS.o\
	CONSTANTS.o TOPPAR_STORE_Mod.o TOPOLOGY_Mod.o TOPOLOGY_GLOBALS.o \
	Class_AtomCnt.o Class_SystemTpg.o TYPES_Utils.o LINKED_LIST.o UNITS.o\
	Class_SystemPrm.o Class_ResidueTpg.o oracS.o 

SRC = $(patsubst %.o, %.f90, $(OBJS)) 

all: oracS.x

oracS.x: $(OBJS) 
	$(LD) -o oracS.x $(OBJS) $(LIBS)

clean:
	rm -f *.o *.x *.mod *~

copy: 
	for j in $(SRC);\
	do cp $$j oppo/$$j;done;


oracS.o: Class_Tree.o ORAC_INPUT_Mod.o GRAMMAR_Mod.o ORAC_PROCESS_Mod.o TOPOLOGY_Mod.o
GRAMMAR_Mod.o: Class_Tree.o STRINGS_Mod.o flib.o TYPES.o ERROR_Mod.o
STRINGS_Mod.o: TYPES.o flib.o 
TOPPAR_STORE_Mod.o: CONSTANTS.o NODES_Mod.o flib.o PARAMETERS_GLOBALS.o
NODES_Mod.o: CONSTANTS.o
ORAC_INPUT_Mod.o: flib.o TYPES.o STRINGS_Mod.o Xerror_Mod.o
ORAC_PROCESS_Mod.o: Class_Tree.o MYPARSE_Mod.o STRINGS_Mod.o ERROR_Mod.o PARAMETERS_Mod.o
Class_Tree.o: flib.o Xerror_Mod.o CONSTANTS.o STRINGS_Mod.o
TOPS_Mod.o: TOPS_Mod.f90 TOPS_Charmm.f90 TOPS_Orac.f90 TYPES.o Xerror_Mod.o MYPARSE_Mod.o\
	STRINGS_Mod.o PARAMETERS_GLOBALS.o
PARAMETERS_Mod.o: CONSTANTS.o PARAMETERS_GLOBALS.o Class_Tree.o ERROR_Mod.o MYPARSE_Mod.o\
	 TOPPAR_STORE_Mod.o TOPS_Mod.o Store_TpgPar.f90
GLOBALS.o: PARAMETERS_GLOBALS.o 
PARAMETERS_GLOBALS.o: TYPES.o
TYPES.o: CONSTANTS.o
MYPARSE_Mod.o: CONSTANTS.o flib.o
ERROR_Mod.o: CONSTANTS.o
TOPOLOGY_Mod.o: PARAMETERS_GLOBALS.o TOPOLOGY_GLOBALS.o Class_AtomCnt.o\
	TYPES_Utils.o STRINGS_Mod.o MYPARSE_Mod.o ERROR_Mod.o Class_SystemTpg.o\
	Class_SystemPrm.o New_Residues.f90 
Class_AtomCnt.o: ERROR_Mod.o PARAMETERS_GLOBALS.o TOPOLOGY_GLOBALS.o flib.o \
	MYPARSE_Mod.o STRINGS_Mod.o LINKED_LIST.o
Class_SystemTpg.o: ERROR_Mod.o CONSTANTS.o Class_AtomCnt.o PARAMETERS_GLOBALS.o\
	TOPOLOGY_GLOBALS.o Class_ResidueTpg.o
TYPES_Utils.o: TYPES.o

Class_SystemPrm.o: UNITS.o Class_AtomCnt.o ERROR_Mod.o PARAMETERS_GLOBALS.o TOPOLOGY_GLOBALS.o flib.o \
	MYPARSE_Mod.o STRINGS_Mod.o LINKED_LIST.o

Class_ResidueTpg.o: ERROR_Mod.o PARAMETERS_GLOBALS.o TOPOLOGY_GLOBALS.o \
	MYPARSE_Mod.o LINKED_LIST.o flib.o
