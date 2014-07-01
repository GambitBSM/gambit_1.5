# to compile fastsim from scratch and run it using the zee pythia 8 command file and an ideal version of the atlas detector

gmake
gmake example
./example -cmd ../pythia8/cmndpythia_zee_1000.cmnd -j ../detectors/ideal_atlas.json -o tree_zee_ideal.root -d 3

