CC=g++
FLAGS=`root-config --cflags --libs` -lGenVector `correction config --cflags --ldflags --rpath` -o
CFLAGS=`root-config --cflags` `correction config --cflags` -c -g
LFLAGS=`root-config --libs` -lGenVector `correction config --ldflags --rpath` -g -o
EXEC=loop

all: ZZSlimmer loopers

loopers: ZZLooper CompLooper

ZZSlimmer: ZZSlimmer.cc interface/ZZSlimmerBase.h interface/argparse.h; $(CC) $(FLAGS) $@ $<

ZZLooper: ZZLooper.cc interface/ZZLooperBase.h interface/argparse.h; $(CC) $(FLAGS) $@ $<

CompLooper: CompLooper.o RatioPlotter.o; $(CC) $(LFLAGS) $@ $^

CompLooper.o: CompLooper.cc interface/CompLooperBase.h interface/argparse.h; $(CC) $(CFLAGS) $<

RatioPlotter.o: RatioPlotter.cc RatioPlotter.h; $(CC) $(CFLAGS) $<

#CompSlimmer: CompSlimmer.cc interface/CompLooperBase.h; $(CC) $(FLAGS) $@ $<

cleanall: clean rootclean

clean:
	rm -f *.o
	rm -f ZZSlimmer ZZLooper CompLooper

rootclean:
	rm -f ZZSlimmer_cc*.* ZZLooper_cc*.* CompLooper_cc*.*
