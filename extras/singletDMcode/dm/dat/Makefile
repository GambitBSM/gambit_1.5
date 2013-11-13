# Makefile for dat directory.
# Technical questions about this makefile should be directed to
# Joakim Edsjo, edsjo@physto.se
# Date: August 29, 2000

# This is where the tables will be installed
LIB=../lib

# These are tables that are just copied
TAB_A = simdiff151.dat simdiff152.dat simdiff153.dat simdiff154.dat	\
simdiff155.dat simdiff156.dat simdiff157.dat simdiff158.dat	\
simdiff171.dat simdiff172.dat simdiff173.dat simint51.dat simint52.dat	\
simint53.dat simint54.dat simint55.dat simint56.dat simint57.dat \
simint58.dat simint71.dat simint72.dat simint73.dat

# These are tables that will be converted to binary tables
#TAB_B = simicediff101.dat.gz simicediff102.dat.gz simicediff103.dat.gz	\
simiceint1.dat.gz simiceint2.dat.gz simiceint3.dat.gz

TAB_B = we-res-b-full-fl21-f1.dat.gz \
we-res-b-full-fl21-f2.dat.gz \
we-res-b-full-fl21-f3.dat.gz \
we-res-b-full-fl21-f4.dat.gz \
we-res-b-full-fl21-f5.dat.gz \
we-res-b-full-fl21-f6.dat.gz \
we-res-b-full-fl21-f7.dat.gz \
we-res-b-full-fl21-f8.dat.gz \
we-res-b-full-fl21-f9.dat.gz \
we-res-b-full-fl21-f10.dat.gz \
we-res-b-full-fl21-f11.dat.gz \
we-res-b-full-fl21-f12.dat.gz \
we-res-b-full-fl21-f13.dat.gz \
we-res-b-full-fl21-f14.dat.gz \
we-res-b-full-fl21-f15.dat.gz \
we-res-b-full-fl21-f16.dat.gz \
we-res-b-full-fl21-f17.dat.gz \
we-res-b-full-fl21-f18.dat.gz \
we-res-b-full-fl21-f19.dat.gz \
we-res-b-full-fl21-f20.dat.gz \
we-res-b-full-fl21-f21.dat.gz \
we-res-b-full-fl21-f22.dat.gz \
we-res-b-full-fl21-f23.dat.gz \
we-res-b-full-fl21-f24.dat.gz \
we-res-b-full-fl21-f25.dat.gz \
we-res-b-full-fl21-f26.dat.gz \
we-res-b-full-fl22-f1.dat.gz \
we-res-b-full-fl22-f2.dat.gz \
we-res-b-full-fl22-f3.dat.gz \
we-res-b-full-fl22-f4.dat.gz \
we-res-b-full-fl22-f5.dat.gz \
we-res-b-full-fl22-f6.dat.gz \
we-res-b-full-fl22-f7.dat.gz \
we-res-b-full-fl22-f8.dat.gz \
we-res-b-full-fl22-f9.dat.gz \
we-res-b-full-fl22-f10.dat.gz \
we-res-b-full-fl22-f11.dat.gz \
we-res-b-full-fl22-f12.dat.gz \
we-res-b-full-fl22-f13.dat.gz \
we-res-b-full-fl22-f14.dat.gz \
we-res-b-full-fl22-f15.dat.gz \
we-res-b-full-fl22-f16.dat.gz \
we-res-b-full-fl22-f17.dat.gz \
we-res-b-full-fl22-f18.dat.gz \
we-res-b-full-fl22-f19.dat.gz \
we-res-b-full-fl22-f20.dat.gz \
we-res-b-full-fl22-f21.dat.gz \
we-res-b-full-fl22-f22.dat.gz \
we-res-b-full-fl22-f23.dat.gz \
we-res-b-full-fl22-f24.dat.gz \
we-res-b-full-fl22-f25.dat.gz \
we-res-b-full-fl22-f26.dat.gz


TAB_OUT_A = $(TAB_A:%=$(LIB)/%)
TAB_OUT_B = $(TAB_B:%.dat.gz=$(LIB)/%.bin)

all : ascii2bin $(TAB_OUT_A) $(TAB_OUT_B)

$(TAB_OUT_A) : ascii2bin $(TAB_A)
	perl tab_install.pl $(@F) $(LIB)

$(TAB_OUT_B) : ascii2bin $(TAB_B)
	perl tab_install.pl $(@F:.bin=.dat.gz) $(LIB)

ascii2bin : ascii2bin.f
	$(FC) -o ascii2bin ascii2bin.f

clean :
	rm -f tmp.dat tmp.dat.gz core
