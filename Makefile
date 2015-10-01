include Makefile.inc

BINS	= neighcorr
OBJLIBS	= 
OBJS	= neighcorr.o trajectory.o configuration.o

all: $(BINS) $(OBJLIBS)

neighcorr: $(OBJS)
	$(LINK) neighcorr $(OBJS)

neighcorr.o: neighcorr.cpp trajectory.o
	$(COMP) neighcorr.cpp

trajectory.o: trajectory.cpp trajectory.h configuration.o
	$(COMP) trajectory.cpp

configuration.o: configuration.cpp configuration.h utilities.h
	$(COMP) configuration.cpp

clean:
	@rm -fv *.o *~ *.pyc $(BINS) $(OBJLIBS)

