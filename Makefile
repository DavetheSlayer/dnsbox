COMPILER =  /cluster/home/nbudanur/local2/bin/mpif90 

FLAGS = -DGNU -DMEASURE -DSTRIDE1 -DFFTW -I/cluster/home/nbudanur/local2/include/ -g -O2

LIBS = /cluster/home/nbudanur/local2/lib/libp3dfft.a \
	   /cluster/home/nbudanur/local2/lib/libfftw3.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5hl_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_hl.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5_fortran.a \
	   /clusterhome/nbudanur/local2/lib/libhdf5.a -lz -ldl
	   
SOURCES = parameters.f90 variables.f90 io.f90 rhs.f90 nsbox.f90 

ALL: $(SOURCES)
	 ${COMPILER} -o dns.x $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
