DIR          := models/NSM
MODNAME      := NSM
SARAH_MODEL  := NSM
WITH_$(MODNAME) := yes

NSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NSM_MK     := \
		$(DIR)/module.mk

NSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NSM_INCLUDE_MK := \
		$(NSM_SUSY_BETAS_MK) \
		$(NSM_SOFT_BETAS_MK)

NSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NSM_generated \
		$(DIR)/LesHouches.in.NSM

NSM_GNUPLOT := \
		$(DIR)/NSM_plot_rgflow.gnuplot \
		$(DIR)/NSM_plot_spectrum.gnuplot

NSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNSM_SRC := \
		$(DIR)/NSM_a_muon.cpp \
		$(DIR)/NSM_edm.cpp \
		$(DIR)/NSM_effective_couplings.cpp \
		$(DIR)/NSM_info.cpp \
		$(DIR)/NSM_input_parameters.cpp \
		$(DIR)/NSM_mass_eigenstates.cpp \
		$(DIR)/NSM_observables.cpp \
		$(DIR)/NSM_physical.cpp \
		$(DIR)/NSM_slha_io.cpp \
		$(DIR)/NSM_soft_parameters.cpp \
		$(DIR)/NSM_susy_parameters.cpp \
		$(DIR)/NSM_utilities.cpp \
		$(DIR)/NSM_weinberg_angle.cpp

EXENSM_SRC := \
		$(DIR)/run_NSM.cpp \
		$(DIR)/run_cmd_line_NSM.cpp \
		$(DIR)/scan_NSM.cpp
LLNSM_LIB  :=
LLNSM_OBJ  :=
LLNSM_SRC  := \
		$(DIR)/NSM_librarylink.cpp

LLNSM_MMA  := \
		$(DIR)/NSM_librarylink.m \
		$(DIR)/run_NSM.m

LIBNSM_HDR := \
		$(DIR)/NSM_cxx_diagrams.hpp \
		$(DIR)/NSM_a_muon.hpp \
		$(DIR)/NSM_convergence_tester.hpp \
		$(DIR)/NSM_edm.hpp \
		$(DIR)/NSM_effective_couplings.hpp \
		$(DIR)/NSM_ewsb_solver.hpp \
		$(DIR)/NSM_ewsb_solver_interface.hpp \
		$(DIR)/NSM_high_scale_constraint.hpp \
		$(DIR)/NSM_info.hpp \
		$(DIR)/NSM_initial_guesser.hpp \
		$(DIR)/NSM_input_parameters.hpp \
		$(DIR)/NSM_low_scale_constraint.hpp \
		$(DIR)/NSM_mass_eigenstates.hpp \
		$(DIR)/NSM_model.hpp \
		$(DIR)/NSM_model_slha.hpp \
		$(DIR)/NSM_observables.hpp \
		$(DIR)/NSM_physical.hpp \
		$(DIR)/NSM_slha_io.hpp \
		$(DIR)/NSM_spectrum_generator.hpp \
		$(DIR)/NSM_spectrum_generator_interface.hpp \
		$(DIR)/NSM_soft_parameters.hpp \
		$(DIR)/NSM_susy_parameters.hpp \
		$(DIR)/NSM_susy_scale_constraint.hpp \
		$(DIR)/NSM_utilities.hpp \
		$(DIR)/NSM_weinberg_angle.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
-include $(DIR)/two_scale.mk
endif
ifneq ($(findstring lattice,$(SOLVERS)),)
-include $(DIR)/lattice.mk
endif
ifneq ($(findstring semi_analytic,$(SOLVERS)),)
-include $(DIR)/semi_analytic.mk
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NSM_SUSY_BETAS_MK)
-include $(NSM_SOFT_BETAS_MK)
-include $(NSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
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

# remove duplicates in case all solvers are used
LIBNSM_SRC := $(sort $(LIBNSM_SRC))
EXENSM_SRC := $(sort $(EXENSM_SRC))

LIBNSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNSM_SRC)))

EXENSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENSM_SRC)))

EXENSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENSM_SRC)))

LIBNSM_DEP := \
		$(LIBNSM_OBJ:.o=.d)

EXENSM_DEP := \
		$(EXENSM_OBJ:.o=.d)

LLNSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNSM_SRC)))

LLNSM_OBJ  := $(LLNSM_SRC:.cpp=.o)
LLNSM_LIB  := $(LLNSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNSM) $(EXENSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNSM_SRC) $(NSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNSM_HDR) $(NSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENSM_SRC) $(NSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNSM_SRC) $(NSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNSM_MMA) $(NSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NSM_MK) $(NSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NSM_INCLUDE_MK) $(NSM_INSTALL_DIR)
ifneq ($(NSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NSM_SLHA_INPUT) $(NSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NSM_GNUPLOT) $(NSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNSM_DEP)
		-rm -f $(EXENSM_DEP)
		-rm -f $(LLNSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNSM)
		-rm -f $(LLNSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNSM_OBJ)
		-rm -f $(EXENSM_OBJ)
		-rm -f $(LLNSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NSM_TARBALL) \
		$(LIBNSM_SRC) $(LIBNSM_HDR) \
		$(EXENSM_SRC) \
		$(LLNSM_SRC) $(LLNSM_MMA) \
		$(NSM_MK) $(NSM_INCLUDE_MK) \
		$(NSM_SLHA_INPUT) $(NSM_GNUPLOT)

$(LIBNSM_SRC) $(LIBNSM_HDR) $(EXENSM_SRC) $(LLNSM_SRC) $(LLNSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_NSM)"
		@echo "Note: to regenerate NSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NSM):
		@true
endif

$(LIBNSM_DEP) $(EXENSM_DEP) $(LLNSM_DEP) $(LIBNSM_OBJ) $(EXENSM_OBJ) $(LLNSM_OBJ) $(LLNSM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNSM_DEP) $(EXENSM_DEP) $(LLNSM_DEP) $(LIBNSM_OBJ) $(EXENSM_OBJ) $(LLNSM_OBJ) $(LLNSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNSM_OBJ) $(LLNSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBNSM): $(LIBNSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNSM_LIB): $(LLNSM_OBJ) $(LIBNSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBNSM_DEP) $(EXENSM_DEP)
ALLSRC += $(LIBNSM_SRC) $(EXENSM_SRC)
ALLLIB += $(LIBNSM)
ALLEXE += $(EXENSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNSM_DEP)
ALLSRC += $(LLNSM_SRC)
ALLLL  += $(LLNSM_LIB)
endif
