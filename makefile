include makefile.inc

#BINS	= neighcorr traj-box3d go_big minimal
BINS	= go_big
OBJLIBS	=
#OBJS	= trajectory.o configuration.o container.o
OBJS	= configuration.o

all: $(BINS) $(OBJLIBS)

neighcorr: neighcorr.o
	$(LINK) neighcorr neighcorr.o $(OBJS)
neighcorr.o: neighcorr.cc optionparser.h trajectory.o
	$(COMP) neighcorr.cc


traj-box3d: traj-box3d.o
	$(LINK) traj-box3d traj-box3d.o $(OBJS)
traj-box3d.o: traj-box3d.cc trajectory.o
	$(COMP) traj-box3d.cc


go_big: go_big.o
	$(LINK) go_big go_big.o $(OBJS)
go_big.o: go_big.cc configuration.o
	$(COMP) go_big.cc


minimal: minimal.o
	$(LINK) minimal minimal.o $(OBJS)
minimal.o: minimal.cc trajectory.o
	$(COMP) minimal.cc

trajectory.o: trajectory.cc trajectory.h configuration.o
	$(COMP) trajectory.cc
configuration.o: configuration.cc configuration.h configuration_base.h configuration_view.h subsystems.h
	$(COMP) configuration.cc
subsystems.o: subsystems.cc subsystems.h
	$(COMP) subsystems.cc
container.o: container.cc container.h
	$(COMP) container.cc

clean:
	@rm -fv *.o *~ *.pyc $(BINS) $(OBJLIBS)
