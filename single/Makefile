COMPILER =  gfortran

COMPFLAGS = 
FLAGS = -lfftw3 -lm -lnetcdff -lnetcdf
INCLUDE = /usr/include/
#LIB = /usr/local/lib
LIB = /usr/lib
	   
SOURCES = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsbox.f90 
SOURCESTEST = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 test.f90 

ALL: $(SOURCES)
	 ${COMPILER} ${COMPFLAGS} -o dns.x $(SOURCES) -I$(INCLUDE) $(FLAGS) -L$(LIB)

test: $(SOURCESTEST)
	  ${COMPILER} ${COMPFLAGS} -o test.x $(SOURCESTEST) -I$(INCLUDE) $(FLAGS) -L$(LIB)

clean:
	rm -f *.o
	rm -f *.mod
