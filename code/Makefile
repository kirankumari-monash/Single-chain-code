# Time-stamp: <Makefile 18:44, 14 Oct 2007 by P Sunthar>
# $Log: Makefile,v $
# Revision 1.1  2008/03/15 01:22:19  sunthar
# gfortran
#

CC = cc
FC90 = gfortran
FC = g77

#CC = icc
#FC90 = ifort
#MPIF90 = mpif90
#FC = ifort


## optimisation level
OPT = -O3 #-fast
#OPT = -fast
#OPT = -cm -O3  -xW -w95 -ipo -tpp7
#OPT = -cm -w95
#db = -g
#DEBUG = -mcmodel=large  ## will increase the execution time
#PRECISION = -r8
# for lc at apac
#LAMHOME=/opt/lam-6.5.8/Intel
#INCLUDES = -I$(LAMHOME)/include

# for mpich at brecca.vpac
#MPIHOME = /usr/local/mpich-1.2.5-i-gm
#INCLUDES = -I$(MPIHOME)/include

FFLAGS = -fbounds-check -O3 #-fdefault-real-8
## Compiler flags
#FFLAGS = -O3 -fast -CB -r8 #-mcmodel=large
#FFLAGS = $(db) $(INCLUDES) $(OPT) $(PRECISION) $(DEBUG)\
# -warn unused -warn nousage -hpfwarn hpf \
# -warn declarations -warn argument_checking


FC90FLAGS = $(FFLAGS) 

# profiling while linking
#LDFLAGS = -p

# if ipo was used during compilation
#LDFLAGS = -ipo

MODS = globals.o interfaces.o

MODS = modules.o
LIBS = modules.o pcsi.o utils.o properties.o 


# for ifort in lc
#SYSLIBS = -L$(MKL)/lib/32 \
#       -lmkl_solver -lmkl_lapack -lmkl_ia32 \
#       -lguide -lpthread -lm 

#SYSLIBS =  -llapack -lblas -lm 
SYSLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread


CHAINLIB = libchain.a

sens: modules.o gsipc.o utils.o properties.o sensemble.o 
	$(FC90) -o $@ $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

cens: modules.f90 gsipc.f90 utils.f90 properties.f90 censemb.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

CSstrin: modules.f90 gsipc.f90 utils.f90 properties.f90 Sinit.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

dens: modules.f90 gsipc.f90 utils.f90 properties.f90 debens.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

lodbal: lbal.f90
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)



ecstep: modules.f90 gsipc.f90 utils.f90 properties.f90 ensbl_cstep.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

wfb: fbnk.f90 
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

wb2fb: fbl0.f90 utils.f90
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@

textra: textra.o hcoutils.o
	$(FC90) $(LDFLAGS) -o textra textra.o hcoutils.o 
	strip $@

fbspring: fb.f90 utils.f90
	$(FC90) -o $@ $(FC90FLAGS) $(LDFLAGS) $+ $(SYSLIBS)
	strip $@


.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .F .c

clean:
	rm *.o *.mod 

# Fixed source, if any:
.f.o:
	$(FC) -c $(FFLAGS) $<

.F.o:
	$(FC) -c $(CCOPTS) $<

.f.mod:
	$(FC) -c $(FC90FLAGS) $<

# Free form source:
.f90.o:
	$(FC90) -c $(FC90FLAGS) $<

.f90.mod:
	$(FC90) -c $(FC90FLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<
