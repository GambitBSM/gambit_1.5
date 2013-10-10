gfortran -shared -Wl,-soname,lib_fortran.so lib_fortran.f90 -o lib_fortran.so -fPIC
g++ -I../../modules/Utils/include/ -I../../modules/Backends/include/ backend_tester.cpp -o backend_tester.out -ldl -std=c++0x

