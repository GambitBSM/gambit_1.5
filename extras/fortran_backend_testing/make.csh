# testing libfcode.so

#gfortran -c -fPIC fcode.f90
#gfortran -shared -o libfcode.so fcode.o
#g++ -I../../modules/Utils/include/ fcode_test.cpp -ldl -o fcode_test -std=c++0x

ifort -c -fPIC fcode.f90
ifort -shared -o libfcode.so fcode.o
icpc -I../../modules/Utils/include/ fcode_test.cpp -ldl -o fcode_test -std=c++0x
