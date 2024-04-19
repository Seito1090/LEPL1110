#! /bin/bash
#
cd ProjectPreProcessor/
cmake CMakeLists.txt -B build/
cd build/
make
./myFem
cd ..
rm -rf build
cd ../Project
cmake CMakeLists.txt -B build
cd build 
make 
./myFem
cd ..
rm -rf build
cd ../ProjectPostProcessor
cmake CMakeLists.txt -B build
cd build 
make
./myFem
cd ..
rm -rf build
