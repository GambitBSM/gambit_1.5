# Package information
PKGNAME            := FlexibleSUSY
FLEXIBLESUSY_VERSION := 2.0.0-dev
FLEXIBLESUSY_MAJOR := 2
FLEXIBLESUSY_MINOR := 0
FLEXIBLESUSY_PATCH := 0
FLEXIBLESUSY_EXTRA := -dev
FLEXIBLESUSY_PKG   := $(PKGNAME)-$(FLEXIBLESUSY_VERSION)
FLEXIBLESUSY_TAG   := v$(FLEXIBLESUSY_VERSION)
ABSBASEDIR         := /mnt/xrootdd/gambit/rose/gambit/contrib/MassSpectra/flexiblesusy
INSTALL_DIR        := 

# Makefile switches
ENABLE_COMPILE     := yes
ENABLE_FFLITE      := no
ENABLE_HIMALAYA    := no
ENABLE_LAPACK      := yes
ENABLE_LIBRARYLINK := no
ENABLE_LOOPTOOLS   := no
ENABLE_META        := no
ENABLE_SHARED_LIBS := no
ENABLE_SQLITE      := yes
ENABLE_STATIC      := no
ENABLE_THREADS     := yes
ENABLE_TSIL        := no

# C/C++ preprocessor defines
ENABLE_COLORS         := no
ENABLE_DEBUG          := 
ENABLE_CHECK_EIGENVALUE_ERROR := no
ENABLE_ILP64MKL_WORKAROUND := no
ENABLE_SILENT         := no
ENABLE_VERBOSE        := no

# Makefile modules
MODELS             :=  models/CMSSM models/MSSM models/MSSMatMGUT models/MSSM_mAmu models/MSSMatMSUSY_mAmu models/MSSMatMGUT_mAmu models/SingletDMZ3 models/SingletDM
ADDONS             := 
MODULES            := config src fflite slhaea doc
MODULES            += $(MODELS) $(ADDONS)

# Variables for Mathematica meta code
MATH               := math
SOLVERS            := two_scale semi_analytic

# Variables for compilation
CXX                := /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/6.2.0/x86_64-slc6/bin/g++
CPPFLAGS           :=  $(patsubst %,-I%,$(MODULES)) -I.
CXXFLAGS           :=  -fPIC -std=c++14 -fopenmp -Wall -Wextra -Wno-misleading-indentation -O2 -DNDEBUG -Wno-deprecated-declarations -Wno-unused-parameter -Wno-unused-variable
FC                 := /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/6.2.0/x86_64-slc6/bin/gfortran
FFLAGS             :=  -fPIC -ffree-line-length-none -ffixed-line-length-none -cpp -fopenmp -O2
FLIBS              := 
MODULE_LIBEXT      := .a
MODULE_MAKE_LIB_CMD:= ar cru
LDFLAGS            := -fopenmp
LDLIBS             := 
LIBLNK_LIBEXT      := .so
LIBLNK_MAKE_LIB_CMD:= g++ -shared -o
BLASLIBS           := /usr/lib64/liblapack.so /usr/lib64/libblas.so
BOOSTTESTLIBS      := -lboost_unit_test_framework
BOOSTTHREADLIBS    := 
BOOSTFLAGS         :=  -I/cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/boost-1_62
EIGENFLAGS         := -I/mnt/xrootdd/gambit/rose/gambit/build/../eigen-3.3.3
GSLLIBS            := -L/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/gsl/2.1-x86_64-slc6-gcc62-opt/lib -lgsl -lgslcblas -lm
GSLFLAGS           := -I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/gsl/2.1-x86_64-slc6-gcc62-opt/include
LAPACKLIBS         := /usr/lib64/liblapack.so /usr/lib64/libblas.so
LOOPFUNCFLAGS      := 
LOOPFUNCLIBS       := 
HIMALAYAFLAGS      := 
HIMALAYALIBS       := 
SQLITEFLAGS        := 
SQLITELIBS         := -lsqlite3
THREADLIBS         := -lpthread
TSILFLAGS          := 
TSILLIBS           := 

# the modules add their dependency files to this variable
ALLDEP   :=
# the modules add soucre files to be created to this variable
ALLSRC   :=
# the modules add their libraries to this variable
ALLLIB   :=
# the modules add executables to this variable
ALLEXE   :=
# the modules add LibraryLink libraries to this variable
ALLLL    :=
# the modules add test executables to this variable
ALLTST   :=
# the addons add libraries to this variable to be linked against executables
ADDONLIBS:=

# flexiblesusy-config script
FSCONFIG := ./flexiblesusy-config

# configure script
CONFIGURE_SCRIPT := configure

# dependency generator
DEPGEN = $(DEPGEN_EXE)

# script which lists the SARAH model files
SARAH_DEP_GEN := ./config/list_sarah_model_files.sh

# script which installs a file without export markers
INSTALL_STRIPPED := config/install_stripped.sh

# script which replaces path\like\this with path/like/this
# thereby enabling FlexibleSUSY to run on Cygwin
CONVERT_DOS_PATHS := config/convert_dos_paths.sh

# script which finds Mathematica's include paths
MATH_INC_PATHS := config/mathematica_include_paths.sh

# package information files
README_FILE    := README
AUTHORS_FILE   := AUTHORS
COPYING_FILE   := COPYING
ChangeLog_FILE := ChangeLog

.PHONY:  all allsrc allexec alllib allll alltest clean clean-dep \
	 clean-executables clean-generated clean-lib clean-obj \
	 depend distclean librarylink release-tag release-head \
	 showbuild tag

all:
ifeq ($(ENABLE_META),yes)
all:     allsrc
endif
ifeq ($(ENABLE_COMPILE),yes)
all:     alllib allexec
endif
ifeq ($(ENABLE_LIBRARYLINK),yes)
all:     allll
endif

ifeq ($(INSTALL_DIR),)
install-src::
	$(error Installation directory is not set!  To set in please run: ./configure --with-install-dir=/path/to/install/dir)
else
install-src::
	install -d $(INSTALL_DIR)
	$(INSTALL_STRIPPED) $(CONFIGURE_SCRIPT) $(INSTALL_DIR) -m u=rwx,g=r,o=r
	$(INSTALL_STRIPPED) $(README_FILE) $(INSTALL_DIR) -m u=rw,g=r,o=r
	$(INSTALL_STRIPPED) $(AUTHORS_FILE) $(INSTALL_DIR) -m u=rw,g=r,o=r
	$(INSTALL_STRIPPED) $(COPYING_FILE) $(INSTALL_DIR) -m u=rw,g=r,o=r
	$(INSTALL_STRIPPED) $(ChangeLog_FILE) $(INSTALL_DIR) -m u=rw,g=r,o=r

ifeq ($(ENABLE_META),yes)
install-src:: allsrc
endif
endif

include config/abspathx.mk
include $(patsubst %, %/module.mk, $(MODULES))

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),install-src)
ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),doc)
ifeq ($(ENABLE_COMPILE),yes)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
ifeq ($(findstring pack-,$(MAKECMDGOALS)),)
ifeq ($(findstring print-,$(MAKECMDGOALS)),)
ifeq ($(findstring release-,$(MAKECMDGOALS)),)
-include $(ALLDEP)
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

allsrc:   $(ALLSRC)
allexec:  $(ALLEXE)
alllib:   $(ALLLIB)
allll:    $(ALLLL)
alltest:  $(ALLTST)

clean-dep:
	-rm -f $(ALLDEP)

clean-executables:
	-rm -f $(ALLEXE)

clean-lib:
	-rm -f $(ALLLIB) $(ALLLL)

depend:  clean-dep
depend:  $(ALLDEP)

%.d: %.cpp | $(DEPGEN)
	$(DEPGEN) $(CPPFLAGS) -MM -o '$@' -MT '$*.o' $^

%.d: %.c | $(DEPGEN)
	$(DEPGEN) $(CPPFLAGS) -MM -o '$@' -MT '$*.o' $^

%.d: %.f | $(DEPGEN)
	$(DEPGEN) $(CPPFLAGS) -MM -o '$@' -MT '$*.o' $^

%.d: %.F | $(DEPGEN)
	$(DEPGEN) $(CPPFLAGS) -MM -o '$@' -MT '$*.o' $^

distclean::
	-rm -f Makefile
	-rm -f $(FSCONFIG)
	-rm -f $(SARAH_DEP_GEN)
	-rm -f config.boost config.log config.math config.sarah config.status

librarylink: allll
	@true

showbuild:
	@echo "PKGNAME            = $(PKGNAME)"
	@echo "VERSION            = $(FLEXIBLESUSY_VERSION)"
	@echo "ABSBASEDIR         = $(ABSBASEDIR)"
	@echo "INSTALL_DIR        = $(INSTALL_DIR)"
	@echo ""
	@echo "MATH               = $(MATH)"
	@echo "MODELS             = $(MODELS)"
	@echo "ADDONS             = $(ADDONS)"
	@echo "SOLVERS            = $(SOLVERS)"
	@echo ""
	@echo "CXX                = $(CXX)"
	@echo "CPPFLAGS           = $(CPPFLAGS)"
	@echo "CXXFLAGS           = $(CXXFLAGS)"
	@echo "FC                 = $(FC)"
	@echo "FFLAGS             = $(FFLAGS)"
	@echo "FLIBS              = $(FLIBS)"
	@echo "LDFLAGS            = $(LDFLAGS)"
	@echo "LDLIBS             = $(LDLIBS)"
	@echo "BLASLIBS           = $(BLASLIBS)"
	@echo "BOOSTTESTLIBS      = $(BOOSTTESTLIBS)"
	@echo "BOOSTTHREADLIBS    = $(BOOSTTHREADLIBS)"
	@echo "BOOSTFLAGS         = $(BOOSTFLAGS)"
	@echo "EIGENFLAGS         = $(EIGENFLAGS)"
	@echo "GSLLIBS            = $(GSLLIBS)"
	@echo "GSLFLAGS           = $(GSLFLAGS)"
	@echo "HIMALAYAFLAGS      = $(HIMALAYAFLAGS)"
	@echo "HIMALAYALIBS       = $(HIMALAYALIBS)"
	@echo "LAPACKLIBS         = $(LAPACKLIBS)"
	@echo "LIBLNK_MAKE_LIB_CMD= $(LIBLNK_MAKE_LIB_CMD)"
	@echo "LIBLNK_LIBEXT      = $(LIBLNK_LIBEXT)"
	@echo "LOOPFUNCFLAGS      = $(LOOPFUNCFLAGS)"
	@echo "LOOPFUNCLIBS       = $(LOOPFUNCLIBS)"
	@echo "MODULE_MAKE_LIB_CMD= $(MODULE_MAKE_LIB_CMD)"
	@echo "MODULE_LIBEXT      = $(MODULE_LIBEXT)"
	@echo "SQLITELIBS         = $(SQLITELIBS)"
	@echo "SQLITEFLAGS        = $(SQLITEFLAGS)"
	@echo "THREADLIBS         = $(THREADLIBS)"
	@echo "TSILLIBS           = $(TSILLIBS)"
	@echo "TSILFLAGS          = $(TSILFLAGS)"
	@echo ""
	@echo "ENABLE_COLORS      = $(ENABLE_COLORS)"
	@echo "ENABLE_COMPILE     = $(ENABLE_COMPILE)"
	@echo "ENABLE_DEBUG       = $(ENABLE_DEBUG)"
	@echo "ENABLE_CHECK_EIGENVALUE_ERROR = $(ENABLE_CHECK_EIGENVALUE_ERROR)"
	@echo "ENABLE_FFLITE      = $(ENABLE_FFLITE)"
	@echo "ENABLE_HIMALAYA    = $(ENABLE_HIMALAYA)"
	@echo "ENABLE_LAPACK      = $(ENABLE_LAPACK)"
	@echo "ENABLE_LIBRARYLINK = $(ENABLE_LIBRARYLINK)"
	@echo "ENABLE_ILP64MKL_WORKAROUND = $(ENABLE_ILP64MKL_WORKAROUND)"
	@echo "ENABLE_LOOPTOOLS   = $(ENABLE_LOOPTOOLS)"
	@echo "ENABLE_META        = $(ENABLE_META)"
	@echo "ENABLE_SHARED_LIBS = $(ENABLE_SHARED_LIBS)"
	@echo "ENABLE_SILENT      = $(ENABLE_SILENT)"
	@echo "ENABLE_SQLITE      = $(ENABLE_SQLITE)"
	@echo "ENABLE_STATIC      = $(ENABLE_STATIC)"
	@echo "ENABLE_THREADS     = $(ENABLE_THREADS)"
	@echo "ENABLE_TSIL        = $(ENABLE_TSIL)"
	@echo "ENABLE_VERBOSE     = $(ENABLE_VERBOSE)"
	@echo ""
	@echo "The list of modules to be built:"
	@echo "--------------------------------"
	@echo "$(MODULES)"

print-% : ; @echo $* = $($*)

tag:
	git tag $(FLEXIBLESUSY_TAG)

release-tag:
	git archive --worktree-attributes \
		--prefix=$(FLEXIBLESUSY_PKG)/ \
		--output=$(FLEXIBLESUSY_PKG).tar.gz $(FLEXIBLESUSY_TAG)
	md5sum $(FLEXIBLESUSY_PKG).tar.gz \
		> $(FLEXIBLESUSY_PKG).tar.gz.md5

release-head:
	$(eval GIT_HEAD_DESCR := $(shell git describe --tags HEAD))
	$(eval FLEXIBLESUSY_HEAD_PKG := $(PKGNAME)-$(GIT_HEAD_DESCR))
	git archive --worktree-attributes \
		--prefix=$(FLEXIBLESUSY_HEAD_PKG)/ \
		--output=$(FLEXIBLESUSY_HEAD_PKG).tar.gz HEAD
	md5sum $(FLEXIBLESUSY_HEAD_PKG).tar.gz \
		> $(FLEXIBLESUSY_HEAD_PKG).tar.gz.md5
