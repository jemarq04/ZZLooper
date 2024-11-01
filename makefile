CC=g++
FLAGS=`root-config --cflags --libs` -lGenVector -o
CFLAGS=`root-config --cflags` -c -g
LFLAGS=`root-config --libs` -lGenVector -g -o
EXEC=loop

all: ZZSlimmer loopers

loopers: ZZLooper CompLooper

ZZSlimmer: ZZSlimmer.cc interface/ZZSlimmerBase.h; $(CC) $(FLAGS) $@ $<

ZZLooper: ZZLooper.cc interface/ZZLooperBase.h; $(CC) $(FLAGS) $@ $<

CompLooper: CompLooper.o RatioPlotter.o; $(CC) $(LFLAGS) $@ $^

CompLooper.o: CompLooper.cc interface/CompLooperBase.h; $(CC) $(CFLAGS) $<

RatioPlotter.o: RatioPlotter.cc RatioPlotter.h; $(CC) $(CFLAGS) $<

#CompSlimmer: CompSlimmer.cc interface/CompLooperBase.h; $(CC) $(FLAGS) $@ $<

cleanall: clean rootclean

clean:
	rm -f *.o
	rm -f ZZSlimmer ZZLooper CompLooper

rootclean:
	rm -f ZZSlimmer_cc*.* ZZLooper_cc*.* CompLooper_cc*.*
