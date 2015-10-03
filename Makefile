include Makefile.inc

BINS	= neighcorr config-test
OBJLIBS	= 
OBJS	= trajectory.o configuration.o

all: $(BINS) $(OBJLIBS)

neighcorr: neighcorr.o
	$(LINK) neighcorr neighcorr.o $(OBJS)
neighcorr.o: neighcorr.cpp optionparser.h trajectory.o
	$(COMP) neighcorr.cpp

config-test: config-test.o
	$(LINK) config-test config-test.o $(OBJS)
config-test.o: config-test.cc optionparser.h trajectory.o
	$(COMP) config-test.cc

trajectory.o: trajectory.cpp trajectory.h configuration.o
	$(COMP) trajectory.cpp
configuration.o: configuration.cpp configuration.h species.h utilities.h
	$(COMP) configuration.cpp

clean:
	@rm -fv *.o *~ *.pyc $(BINS) $(OBJLIBS)

