# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header$
# $Log$
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


# Dependencies
CONSTANTS.o : CONSTANTS.f90 
Class_AtomCnt.o : Class_AtomCnt.f90 LINKED_LIST.o MYPARSE_Mod.o STRINGS_Mod.o flib.o\
	ERROR_Mod.o TOPOLOGY_GLOBALS.o PARAMETERS_GLOBALS.o 
Class_ResidueTpg.o : Class_ResidueTpg.f90 LINKED_LIST.o MYPARSE_Mod.o\
	ERROR_Mod.o TOPOLOGY_GLOBALS.o PARAMETERS_GLOBALS.o flib.o STRINGS_Mod.o 
Class_SystemPrm.o : Class_SystemPrm.f90 Class_SystemTpg.o LINKED_LIST.o\
	flib.o UNITS.o CONSTANTS.o STRINGS_Mod.o MYPARSE_Mod.o ERROR_Mod.o\
	TOPOLOGY_GLOBALS.o PARAMETERS_GLOBALS.o 
Class_SystemTpg.o : Class_SystemTpg.f90 LINKED_LIST.o CONSTANTS.o Class_ResidueTpg.o\
	Class_AtomCnt.o ERROR_Mod.o TOPOLOGY_GLOBALS.o PARAMETERS_GLOBALS.o 
Class_Tree.o : Class_Tree.f90 Xerror_Mod.o STRINGS_Mod.o CONSTANTS.o 
ERROR_Mod.o : ERROR_Mod.f90 TYPES.o CONSTANTS.o 
GLOBALS.o : GLOBALS.f90 PARAMETERS_GLOBALS.o TYPES.o 
GRAMMAR_Mod.o : GRAMMAR_Mod.f90 flib.o STRINGS_Mod.o Class_Tree.o TYPES.o ERROR_Mod.o 
LINKED_LIST.o : LINKED_LIST.f90 CONSTANTS.o 
MYPARSE_Mod.o : MYPARSE_Mod.f90 flib.o CONSTANTS.o 
NODES_Mod.o : NODES_Mod.f90 CONSTANTS.o 
New_Residues.o : New_Residues.f90 
ORAC_INPUT_Mod.o : ORAC_INPUT_Mod.f90 STRINGS_Mod.o TYPES.o flib.o Xerror_Mod.o 
ORAC_PROCESS_Mod.o : ORAC_PROCESS_Mod.f90 Class_Tree.o MYPARSE_Mod.o\
	ERROR_Mod.o PARAMETERS_Mod.o 
PARAMETERS_GLOBALS.o : PARAMETERS_GLOBALS.f90 TYPES.o 
PARAMETERS_Mod.o : PARAMETERS_Mod.f90 Store_TpgPar.f90 TOPPAR_STORE_Mod.o\
	TOPS_Mod.o MYPARSE_Mod.o STRINGS_Mod.o ERROR_Mod.o Class_Tree.o \
	PARAMETERS_GLOBALS.o CONSTANTS.o 
STRINGS_Mod.o : STRINGS_Mod.f90 flib.o TYPES.o CONSTANTS.o 
Store_TpgPar.o : Store_TpgPar.f90 
TOPOLOGY_GLOBALS.o : TOPOLOGY_GLOBALS.f90 TYPES.o 
TOPOLOGY_Mod.o : TOPOLOGY_Mod.f90 New_Residues.f90 TOPOLOGY_GLOBALS.o\
	PARAMETERS_GLOBALS.o ERROR_Mod.o TYPES_Utils.o MYPARSE_Mod.o STRINGS_Mod.o\
	Class_SystemPrm.o Class_SystemTpg.o 
TOPPAR_STORE_Mod.o : TOPPAR_STORE_Mod.f90 flib.o NODES_Mod.o CONSTANTS.o\
	PARAMETERS_GLOBALS.o 
TOPS_Charmm.o : TOPS_Charmm.f90 
TOPS_Mod.o : TOPS_Mod.f90 TOPS_Charmm.f90 TOPS_Orac.f90 STRINGS_Mod.o MYPARSE_Mod.o\
	Xerror_Mod.o PARAMETERS_GLOBALS.o TYPES.o 
TOPS_Orac.o : TOPS_Orac.f90 
TYPES.o : TYPES.f90 CONSTANTS.o 
TYPES_Utils.o : TYPES_Utils.f90 TYPES.o 
UNITS.o : UNITS.f90 
Xerror_Mod.o : Xerror_Mod.f90 
flib.o : flib.f90 
oracS.o : oracS.f90 TOPOLOGY_Mod.o ORAC_PROCESS_Mod.o GRAMMAR_Mod.o\
	ORAC_INPUT_Mod.o Class_Tree.o