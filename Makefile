COMPILER =  /cluster/home/nbudanur/local2/bin/mpif90 
COMPFLAGS	= -ffree-line-length-none -x f95-cpp-input -c -O4 \
			  -DGNU -DMEASURE -DSTRIDE1 -DFFTW \
			  -I/cluster/home/nbudanur/local2/include/

MODSOBJ = parameters.o variables.o state.o rhs.o io.o

LIBS = /cluster/home/nbudanur/local2/lib/libp3dfft.a \
	   /cluster/home/nbudanur/local2/lib/libfftw3.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5hl_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_hl.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5.a -lz -ldl
	   
SOURCES = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsboxheun.f90 
SOURCESHEUN = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsboxheun.f90 
SOURCESBAND = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsboxband.f90 
SOURCESIMP = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsboximplicit.f90 
SOURCESADJ = parameters.f90 variablesAdj.f90 state.f90 rhsAdj.f90 io.f90 dnsboxAdj.f90 
SOURCESPOST = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 post.f90 

all: $(MODSOBJ) dnsboxheun.f90
	$(COMPILER) $(COMPFLAGS) dnsboxheun.f90
	$(COMPILER) -o ./dns.x dnsboxheun.o $(MODSOBJ) $(FLAGS) $(LIBS)

#ALL: $(SOURCES)
#	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCES) $(LIBS)

test: $(MODSOBJ) test.f90
	$(COMPILER) $(COMPFLAGS) test.f90
	$(COMPILER) -o ./test.x test.o $(MODSOBJ) $(FLAGS) $(LIBS)
 

clean:
	rm -f *.o
	rm -f *.mod

adj: $(SOURCESADJ)
	 ${COMPILER} -o dnsadj.x $(FLAGS) $(SOURCESADJ) $(LIBS)

post: $(SOURCESPOST)
	  ${COMPILER} -o post.x $(FLAGS) $(SOURCESPOST) $(LIBS)

heun: $(SOURCESHEUN)
	  ${COMPILER} -o dns.x $(FLAGS) $(SOURCESHEUN) $(LIBS)

band: $(SOURCESBAND)
	  ${COMPILER} -o dns.x $(FLAGS) $(SOURCESBAND) $(LIBS)

imp: $(SOURCESIMP)
	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCESIMP) $(LIBS)

#------------------------------------------------------------------------------
parameters.o : parameters.f90 
	$(COMPILER) $(COMPFLAGS) parameters.f90

variables.o : variables.f90 parameters.o
	$(COMPILER) $(COMPFLAGS) variables.f90

state.o: state.f90 variables.o
	$(COMPILER) $(COMPFLAGS) state.f90

rhs.o: rhs.f90 variables.o
	$(COMPILER) $(COMPFLAGS) rhs.f90

io.o : io.f90 state.o rhs.o
	$(COMPILER) $(COMPFLAGS) io.f90
