rm -rf pythia8186
cp -r pythia8186_original pythia8186

cp headers_by_hand/identification_pythia.hpp pythia_BOSS_output/backend_types/BOSSedPythia_1_0/identification.hpp
cp headers_by_hand/loaded_types_pythia.hpp pythia_BOSS_output/backend_types/BOSSedPythia_1_0/loaded_types.hpp

cp pythia_BOSS_output/*.h pythia8186/include/
cp pythia_BOSS_output/*.hpp pythia8186/include/

mv pythia8186/include/Pythia.h pythia8186/include/Pythia8/
mv pythia8186/include/Basics.h pythia8186/include/Pythia8/
mv pythia8186/include/Event.h pythia8186/include/Pythia8/
mv pythia8186/include/Info.h pythia8186/include/Pythia8/

cp -r pythia_BOSS_output/backend_types pythia8186/include/

cp -r pythia_BOSS_output/*.cc pythia8186/src/

