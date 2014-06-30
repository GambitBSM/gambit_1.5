PY8_CXXFLAGS=`pythia8-config --cxxflags`
PY8_LDFLAGS=`pythia8-config --libs`

ROOT_CXXFLAGS=`root-config --cflags`
ROOT_LDFLAGS=`root-config --libs`

FJ_CXXFLAGS=`fastjet-config --cxxflags`
FJ_LDFLAGS=`fastjet-config --libs`

SIMPLE_HEP_INC=/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/

FASTSIM_INC=/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include 

#g++ -c -fopenmp -std=c++0x -Wall -g -O3 -pedantic -fPIC -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include FastSim.cpp

#g++ -c -fopenmp -std=c++0x -Wall -g -O3 -pedantic -fPIC -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include DetectorResponse.cpp

#g++ -c -fopenmp -std=c++0x -Wall -g -O3 -pedantic -fPIC -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include fastsim_interface.cpp

#g++ -c -fopenmp -std=c++0x -Wall -g -O3 -pedantic -fPIC -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include FastSim_Reader.cpp

#g++ -c -fopenmp -std=c++0x -Wall -g -O3 -pedantic -fPIC -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include jsoncpp.cc

#g++ -shared FastSim.o DetectorResponse.o fastsim_interface.o FastSim_Reader.o -Wl,-whole-archive -L/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/lib -l fastjet.a -Wl,-no-whole-archive -o libfastsim.so

#g++ -Wall -std=c++0x -g -O3 -pedantic -fPIC  -I${SIMPLE_HEP_INC} -I${FASTSIM_INC} ${PY8_CXXFLAGS} ${ROOT_CXXFLAGS} ${FJ_CXXFLAGS} -L. -lfastsim ${PY8_LDFLAGS} ${ROOT_LDFLAGS} ${FJ_LDFLAGS} example-fastsim-py8.cpp -o fastsimpythia  

#g++ -fopenmp -Wall -std=c++0x -g -O3 -pedantic -fPIC  -I${SIMPLE_HEP_INC} -I${FASTSIM_INC} ${PY8_CXXFLAGS} ${ROOT_CXXFLAGS} ${FJ_CXXFLAGS} -L. -lfastsim ${PY8_LDFLAGS} ${ROOT_LDFLAGS} ${FJ_LDFLAGS} example-fastsim-py8-parallel.cpp -o fastsimpythia-parallel  

#g++ -shared -Wall -fPIC -g -I/import/sydpp2/atlas/Aldos_Data/gambit/external/fastjet/include/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include -o libAldo.so fastsim_interface.cpp
g++ -shared -Wall -fPIC -g -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/modules/contrib/hep_simple_lib/ -I/import/sydpp2/atlas/Aldos_Data/gambit/hepforge2/extras/fast_sim/include -o libAldo.so fastsim_interface.cpp
