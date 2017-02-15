mkdir build
cd build
cmake ../src
make -j8
cd ..
cd driver;make;cd -
cd gridGen;make;cd -
cd case/;./run.sh 16;cd -
