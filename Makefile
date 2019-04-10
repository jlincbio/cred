CC=gcc
CFLAGS=-Wall -Wno-unused-variable -O3 -DNDEBUG -march=native -std=c99

HTSDIR=htslib
INCLUDES=-I$(HTSDIR)
LFLAGS=-L$(HTSDIR)
LIBS=$(HTSDIR)/libhts.a
LDLIBS=-lm -lz -llzma -lbz2 -lpthread

MAIN = cred

all: HTSLIB
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) $(MAIN).c -o $(MAIN) $(LIBS) $(LDLIBS)

HTSLIB:
	$(MAKE) -C $(HTSDIR)/

install:
	cp $(MAIN) $(PREFIX)/
	cp bcred.pl $(PREFIX)/bcred
	chmod +x $(PREFIX)/bcred

clean:
	rm -f $(MAIN) *.o
	rm -rf $(MAIN).dSYM
	rm -f bcred
	
# DO NOT DELETE THIS LINE -- make depend needs it
