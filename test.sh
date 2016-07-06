make clean
make
cp dns.x test/1015/ 
cd test/1015/
/clusterhome/nbudanur/local2/bin/mpirun -np 4 ./dns.x &> test.out &
cd ..
cd ..
