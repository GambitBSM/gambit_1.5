# just compiles the ReadDelphesTRee code
  cmd="g++ ReadDelphesTree.cc"
  cmd+=" -I$(root-config --incdir)"
  cmd+=" $(root-config --libs)"
  cmd+=" -I$DELPHESLOCATION"
  cmd+=" -L$DELPHESLOCATION"
  cmd+=" -lDelphes"
  cmd+=" -o ReadDelphesTree"

  # Compile the script:
  echo "$cmd"
  $cmd

