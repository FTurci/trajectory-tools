include Makefile.inc

BINS	= neighcorr config-test
OBJLIBS	= 
OBJS	= trajectory.o configuration.o

all: $(BINS) $(OBJLIBS)

neighcorr: neighcorr.o
	$(LINK) neighcorr neighcorr.o $(OBJS)
neighcorr.o: neighcorr.cc optionparser.h trajectory.o
	$(COMP) neighcorr.cc

config-test: config-test.o
	$(LINK) config-test config-test.o $(OBJS)
config-test.o: config-test.cc optionparser.h trajectory.o
	$(COMP) config-test.cc

trajectory.o: trajectory.cc trajectory.h configuration.o
	$(COMP) trajectory.cc
configuration.o: configuration.cc configuration.h species.h utilities.h
	$(COMP) configuration.cc

clean:
	@rm -fv *.o *~ *.pyc $(BINS) $(OBJLIBS)

