cd grid;
p2=`expr $1 \* 2`
echo $p2
mpirun -np $p2 ../../build4/gridGen/buildGrid
cd ..
mpirun -np $1 ../build4/driver/tioga.exe
