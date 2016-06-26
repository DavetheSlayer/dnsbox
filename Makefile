COMPILER =  mpif90 

FLAGS = -DGNU -DMEASURE -DSTRIDE1 -DFFTW -I/cluster/home/nbudanur/local/include/  -g -O2

LIBS = /cluster/home/nbudanur/local/lib/libp3dfft.a /cluster/home/nbudanur/local/lib/libfftw3.a

SOURCES = parameters.f90 variables.f90 rhs.f90 nsbox.f90 

ALL: $(SOURCES)
	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
