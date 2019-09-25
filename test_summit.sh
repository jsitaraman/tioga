module load gcc/7.4.0 cmake/3.15.2
mkdir build
cd build
cmake -DTIOGA_HAS_NODEGID:BOOL=off -DTIOGA_ENABLE_TIMERS:BOOL=on ..
make
cd ..
cd case/;./run_summit.sh 8;cd -
