#-------------------------------------------------------------------\
#                                                                   | 
#  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>   |
#                                                                   |
#     This program is free software;  you  can  redistribute  it    |
#     and/or modify it under the terms of the GNU General Public    |
#     License version 2 as published  by  the  Free  Software       |
#     Foundation.                                                   |
#                                                                   |
#     This program is distributed in the hope that  it  will  be    |
#     useful, but WITHOUT ANY WARRANTY; without even the implied    |
#     warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR    |
#     PURPOSE.   See  the GNU General  Public License for more      |
#     details.                                                      |
#                                                                   |
#     You should have received a copy of the GNU General  Public    |
#     License along with this program; if not, write to the Free    |
#     Software Foundation, Inc., 59  Temple  Place,  Suite  330,    |
#     Boston, MA 02111-1307  USA                                    |
#                                                                   |
# ------------------------------------------------------------------/
# Makefile for 
# Copyright (C) 2006 Massimo Marchi <marchi@job.saclay.cea.fr>
# $Header$
# $Log$
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
	SystemPrm.o LennardJones.o PrmUtilities.o BondedPrm.o SimulationBox.o\
	ReadStore.o AddHydrogens_.o oracDD.o 

OBJS_f = AddH.o duni.o

SRC = $(patsubst %.o, %.f90, $(OBJS)) $(patsubst %.o, %.f, $(OBJS_f)) 

all: oracDD.x

oracDD.x: $(OBJS) $(OBJS_f)
	$(LD) -o oracDD.x $(OBJS) $(OBJS_f) $(LIBS)

clean:
	rm -f *.o *.x *.mod *~

copy: 
	for j in $(SRC);\
	do cp $$j oppo/$$j;done;


# Dependencies
io.h: $(shell ../bin/makedepf90 $(SRC) > io.h)
include io.h
$(shell rm io.h)
