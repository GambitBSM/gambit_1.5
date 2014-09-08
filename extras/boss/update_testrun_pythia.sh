
cp pythia8186/lib/libpythia8.so testrun_pythia/

cp output/*.h testrun_pythia/loader/
#cp output/GAMBIT_wrapper_functions.cpp testrun_pythia/loader

rm testrun_pythia/loader/Basics.h
rm testrun_pythia/loader/Event.h
rm testrun_pythia/loader/Info.h
rm testrun_pythia/loader/Pythia.h

