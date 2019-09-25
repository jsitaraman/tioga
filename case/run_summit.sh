cd grid;
p2=`expr $1 \* 2`
echo $p2
jsrun -n $p2 ../../build/gridGen/buildGrid
cd ..
jsrun -n $1 ../build/driver/tioga.exe
