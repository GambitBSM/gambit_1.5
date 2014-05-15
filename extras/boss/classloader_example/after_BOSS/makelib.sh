echo
echo "Compiling..."
echo 
g++ -c classes.cpp functions.cpp factories.cpp -fPIC
echo "...done."
echo

echo "Linking..."
echo
g++ -shared classes.o functions.o factories.o -o lib.so
echo "...done."
echo

echo "Cleaning..."
echo
rm classes.o factories.o
echo "...done."
echo
