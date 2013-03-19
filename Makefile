PROCESSOR := $(shell uname -m)

#ifeq ($(PROCESSOR),ia64)
  F90=gfortran
  FFLAGS=-g -C -O0 -ffree-form -fcheck=all -I/opt/local/include 
  FFLAGS2=$(FFLAGS)
  LDFLAGS=-L/opt/local/lib -lnetcdff
#else
#  ifeq ($(PROCESSOR),x86_64)
#    F90=/usr/local/pgi/linux86-64/5.2/bin/pgf90
#    FFLAGS=-Kieee -fastsse #-Kieee # 
#    LDFLAGS=-lnetcdf
#  else
#    F90=lf95
#    FFLAGS=-g --chk aesux --trap -I/usr/local/include #--o2  -I/usr/local/include #
#    LDFLAGS=-lnetcdf
#  endif
#endif

SOURCES= DGmod.f90 dgsweep.f90
OBJECTS=$(SOURCES:.f90=.o)

all: $(SOURCES) nodal_test

nodal_test: $(OBJECTS) nodal_test.f90
	$(F90) $(FFLAGS) $(OBJECTS) nodal_test.f90 -o $@ $(LDFLAGS) 

.PHONY:

clean:
	rm -f *.o nodal_test

%.o : %.f90
	$(F90) -c $(FFLAGS) $<
