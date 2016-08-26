COMPILER =  /cluster/home/nbudanur/local2/bin/mpif90 

FLAGS = -DGNU -DMEASURE -DSTRIDE1 -DFFTW -I/cluster/home/nbudanur/local2/include/ -g -O2

LIBS = /cluster/home/nbudanur/local2/lib/libp3dfft.a \
	   /cluster/home/nbudanur/local2/lib/libfftw3.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5hl_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_hl.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5.a -lz -ldl
	   
SOURCES = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsbox.f90 
SOURCESHEUN = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 dnsboxheun.f90 
SOURCESADJ = parameters.f90 variablesAdj.f90 state.f90 rhsAdj.f90 io.f90 dnsboxAdj.f90 
SOURCESPOST = parameters.f90 variables.f90 state.f90 rhs.f90 io.f90 post.f90 

ALL: $(SOURCES)
	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod

adj: $(SOURCESADJ)
	 ${COMPILER} -o dnsadj.x $(FLAGS) $(SOURCESADJ) $(LIBS)

post: $(SOURCESPOST)
	  ${COMPILER} -o post.x $(FLAGS) $(SOURCESPOST) $(LIBS)

heun: $(SOURCESHEUN)
	  ${COMPILER} -o dns.x $(FLAGS) $(SOURCESHEUN) $(LIBS)
