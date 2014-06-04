echo
echo "Compiling..."
echo 
g++ -c -pedantic classes.cpp functions.cpp factory_U.cpp factory_T.cpp factory_X.cpp factory_Container.cpp -fPIC
echo "...done."
echo

echo "Linking..."
echo
g++ -shared classes.o functions.o factory_U.o factory_T.o factory_X.o factory_Container.o -o lib.so
echo "...done."
echo

echo "Cleaning..."
echo
rm classes.o functions.o factory_U.o factory_T.o factory_X.o factory_Container.o
echo "...done."
echo
