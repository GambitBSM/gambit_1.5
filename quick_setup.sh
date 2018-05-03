source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/6.2.0/x86_64-slc6/setup.sh
export PYTHONPATH="$PWD/PyYAML-3.11/lib:$PYTHONPATH"

lsetup "gsl 2.1-x86_64-slc6-gcc62-opt"
lsetup "hdf5 1.10.0.patch1-slc6"

CMAKENAME=cmake-3.3.2-Linux-x86_64
if test ! -d $CMAKENAME; then
  wget https://cmake.org/files/v3.3/$CMAKENAME.tar.gz
  tar xf $CMAKENAME.tar.gz
  rm -f $CMAKENAME.tar.gz
fi
export PATH=$PWD/$CMAKENAME/bin:$PATH

PYYAMLNAME=PyYAML-3.11
if test ! -d $PYYAMLNAME; then
  wget http://pyyaml.org/download/pyyaml/$PYYAMLNAME.tar.gz
  tar xf $PYYAMLNAME.tar.gz
  rm -f $PYYAMLNAME.tar.gz
fi
export PYTHONPATH="$PWD/$PYYAMLNAME/lib:$PYTHONPATH"

ename=3.3.3
EIGENNAME=eigen-$ename
if [[ ! -e $EIGENNAME ]]; then
  wget http://bitbucket.org/eigen/eigen/get/$ename.tar.gz
  tar xvf $ename.tar.gz
  rm -f $ename.tar.gz
  ename=$(ls -d eigen-eigen-*)
  ln -sf $ename $EIGENNAME
fi

#LAPACKNAME=lapack-3.7.0
#if test ! -d $LAPACKNAME; then
#    wget http://www.netlib.org/lapack/lapack-3.7.0.tgz
#    tar xvf $LAPACKNAME.tgz
#    rm -f $LAPACKNAME.tgz
#fi
#export PATH=$PWD/$LAPACKNAME/bin:$PATH

#only used when first compiled or remove dependencies after failed build

export OMP_NUM_THREADS=1
