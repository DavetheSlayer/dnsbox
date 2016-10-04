COMPILER =  gfortran

FLAGS = -lfftw3 
INCLUDE = /usr/local/include/
LIB = /usr/local/lib
	   
SOURCES = parameters.f90 variables.f90 state.f90 dnsbox.f90 

ALL: $(SOURCES)
	 ${COMPILER} -o dns.x $(SOURCES) -I$(INCLUDE) $(FLAGS) -L$(LIB)

clean:
	rm -f *.o
	rm -f *.mod
