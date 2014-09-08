cp output/*.hpp minimal/
cp output/*.cpp minimal/
g++ -shared minimal/*.cpp -fPIC -o minimal/libminimal.so
