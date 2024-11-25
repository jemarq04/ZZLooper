CC=g++
FLAGS=`root-config --cflags --libs` -lGenVector -o 
CORR_FLAGS=`correction config --cflags --ldflags --rpath`
CFLAGS=`root-config --cflags` -c -g
CORR_CFLAGS=`correction config --cflags`
LFLAGS=`root-config --libs` -lGenVector -g -o
CORR_LFLAGS=`correction config --ldflags --rpath`
EXEC=loop

ifndef NOSF
	FLAGS:=$(CORR_FLAGS) $(FLAGS)
	CFLAGS:=$(CORR_CFLAGS) $(CFLAGS)
	LFLAGS:=$(CORR_LFLAGS) $(LFLAGS)
else
	FLAGS:=-DNOSF $(FLAGS)
	CFLAGS:=-DNOSF $(CFLAGS)
endif

all: ZZSlimmer loopers

loopers: ZZLooper CompLooper

ZZSlimmer: ZZSlimmer.cc interface/ZZSlimmerBase.h interface/argparse.h; $(CC) $(FLAGS) $@ $<

ZZLooper: ZZLooper.cc interface/ZZLooperBase.h interface/argparse.h; $(CC) $(FLAGS) $@ $<

CompLooper: CompLooper.o RatioPlotter.o; $(CC) $(LFLAGS) $@ $^

CompLooper.o: CompLooper.cc interface/CompLooperBase.h interface/argparse.h; $(CC) $(CFLAGS) $<

RatioPlotter.o: RatioPlotter.cc RatioPlotter.h; $(CC) $(CFLAGS) $<

cleanall: clean rootclean

clean:
	rm -f *.o
	rm -f ZZSlimmer ZZLooper CompLooper

rootclean:
	rm -f ZZSlimmer_cc*.* ZZLooper_cc*.* CompLooper_cc*.*
