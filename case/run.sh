cd grid;
mpirun -np $1 ../../gridGen/buildGrid
cd ..
mpirun -np $1 ../driver/tioga.exe
