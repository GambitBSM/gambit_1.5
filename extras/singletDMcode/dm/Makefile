
SPACE:=
SPACE+=
ESCAPED:=\$(SPACE)
DM_ROOT=$(subst $(SPACE),$(ESCAPED),$(PWD)/)

dm : config_paths dm_lib install_tables


all : dm

config_paths :
	perl $(DM_ROOT)scr/config.pl $(DM_ROOT)

dm_lib dmlib :
	$(MAKE) -C src

##############################


install_tables :
	cd $(DM_ROOT)dat; make all


clean :
	cd $(DM_ROOT)src; $(MAKE) clean; \
	rm -f $(DM_ROOT)lib/sim* $(DM_ROOT)lib/*.bin; \
	rm -f $(DM_ROOT)dat/ascii2bin	




.PHONY: all config_paths dm_lib dmlib install_tables \
	clean
