icpc -c -fPIC libfirst.cpp
icpc -shared -o libfirst.so libfirst.o  
icpc factory_test.cpp -o factory_test

#g++ -c -fPIC libfirst.cpp
#g++ -shared -o libfirst.so libfirst.o  
#g++ factory_test.cpp -ldl -o factory_test
