# Makefile for the DarkSUSY library
# Technical questions about this makefile should be directed to
# Joakim Edsjo, edsjo@physto.se
# Date: August 29, 2000.

# You shouldn't have to change anything in this file.
# You should only change in the makefile in the DarkSUSY root directory.

###############################################
#### DO NOT CHANGE ANYTHING BELOW THIS LINE ###
###############################################

# Here we set where the library goes and where include files are.
# NOTE. These need to be specified with relative pathnames to be
# used properly by the makefiles in the underlying directories.
LIB=../lib
INC=../inc


vpath %.a $(LIB)

DIRS = an ge ha hm hr ib ini xcern xcmlib interf 

DOBJ = $(DIRS:%=%_obj)

OBJWDIR = */*.o

all : libdm.a

libdm.a : $(DOBJ)
	ar rv $(LIB)/libdm.a $(OBJWDIR) 
	ranlib $(LIB)/libdm.a

$(DOBJ) : 
	cd $(@:_obj=) ; make all

.PHONY: clean $(DOBJ)

clean : 
	rm -f */*.o; rm -f  $(LIB)/libdm.a


