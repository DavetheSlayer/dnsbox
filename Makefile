COMPILER =  /cluster/home/nbudanur/local2/bin/mpif90 

FLAGS = -DGNU -DMEASURE -DSTRIDE1 -DFFTW -I/cluster/home/nbudanur/local2/include/ -g -O2

LIBS = /cluster/home/nbudanur/local2/lib/libp3dfft.a \
	   /cluster/home/nbudanur/local2/lib/libfftw3.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5hl_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_hl.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5.a -lz -ldl
	   
SOURCES = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsbox.f90 
SOURCESADJ = parameters.f90 variablesAdj.f90 state.f90 rhsAdj.f90 io.f90 dnsboxAdj.f90 
SOURCESRK = parameters.f90 variables.f90 io.f90 rhs.f90 nsboxRK4.f90 
SOURCESETDRK = parameters.f90 variables.f90 io.f90 rhs.f90 nsboxETDRK4.f90 
SOURCESETD = parameters.f90 variables.f90 io.f90 rhs.f90 nsboxETD.f90 
SOURCESEULER = parameters.f90 variables.f90 io.f90 rhs.f90 nsboxEuler.f90 
SOURCESVAR = parameters.f90 variables.f90 testvar.f90
SOURCESRHS = parameters.f90 variables.f90 rhs.f90 testrhs.f90
SOURCESSPEC = parameters.f90 variables.f90 io.f90 rhs.f90 testspec.f90 
SOURCESDEAL = parameters.f90 variables.f90 io.f90 rhs.f90 testdealias.f90 
SOURCESEKIN = parameters.f90 variables.f90 rhs.f90 io.f90 testekin.f90 
SOURCESSTATE = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 teststate.f90 

ALL: $(SOURCES)
	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCES) $(LIBS)

adj: $(SOURCESADJ)
	 ${COMPILER} -o dnsadj.x $(FLAGS) $(SOURCESADJ) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod

var: $(SOURCESVAR)
	 ${COMPILER} -o testvar.x $(FLAGS) $(SOURCESVAR) $(LIBS)

rhs: $(SOURCESRHS)
	 ${COMPILER} -o testrhs.x $(FLAGS) $(SOURCESRHS) $(LIBS)

rk: $(SOURCESRK)
	${COMPILER} -o testRK.x $(FLAGS) $(SOURCESRK) $(LIBS)

etdrk: $(SOURCESETDRK)
	${COMPILER} -o testETDRK.x $(FLAGS) $(SOURCESETDRK) $(LIBS)

etd: $(SOURCESETD)
	${COMPILER} -o testETD.x $(FLAGS) $(SOURCESETD) $(LIBS)

euler: $(SOURCESEULER)
	${COMPILER} -o testEuler.x $(FLAGS) $(SOURCESEULER) $(LIBS)

spec: $(SOURCESSPEC)
	 ${COMPILER} -o testSpec.x $(FLAGS) $(SOURCESSPEC) $(LIBS)

deal: $(SOURCESDEAL)
	 ${COMPILER} -o testDeal.x $(FLAGS) $(SOURCESDEAL) $(LIBS)
	 
ekin: $(SOURCESEKIN)
	 ${COMPILER} -o testEkin.x $(FLAGS) $(SOURCESEKIN) $(LIBS)

state: $(SOURCESSTATE)
	 ${COMPILER} -o testSkin.x $(FLAGS) $(SOURCESSTATE) $(LIBS)
