cd src;make clean;cd ..
cd driver;make clean;cd ..
cd gridGen;make clean;cd ..
cd case/grid;rm -rf *.plt 2>/dev/null;cd -
rm -rf case/flow*.dat 2> /dev/null
