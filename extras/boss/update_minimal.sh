cp output/additions_for_external_code/*.hpp minimal/
cp output/additions_for_external_code/*.cpp minimal/
g++ -shared minimal/*.cpp -fPIC -o minimal/libminimal.so
