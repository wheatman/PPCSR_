.PHONY:	all clean

all:
	./setup.sh make -B -f Makefile2 $(MAKECMDGOALS) -j 24

clean:
	./setup.sh make -f Makefile2 clean
