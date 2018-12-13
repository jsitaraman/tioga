cd grid;
p2=`expr $1 \* 2`
echo $p2
mpirun -np $p2 ../../build/gridGen/buildGrid
cd ..
mpirun -np $1 ../build/driver/tioga.exe
