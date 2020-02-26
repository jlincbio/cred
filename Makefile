CC=gcc
OSET=-O3

PREFIX=/usr/local/bin
LOCAL_DIR=$(CURDIR)
HTSDIR=$(CURDIR)/htslib
INCLUDES=-I$(HTSDIR)
LFLAGS=-L$(HTSDIR)
LIBS=$(HTSDIR)/libhts.a

CFLAGS=-Wall -Wno-unused-variable $(OSET) -DNDEBUG -march=native -std=c99

# HTSlib versions 1.10.2 onward appears to require external libcurl linkage
# tested configuration: GNU Make 4.2.1, GCC 9.2.1 
LDLIBS=-lm -lz -llzma -lbz2 -lpthread -lcurl

MAIN=cred

all: HTSLIB 
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) $(MAIN).c -o $(MAIN) $(LIBS) $(LDLIBS)

HTSLIB:
	$(MAKE) -C $(HTSDIR)/

darwin: macprep all
	
macprep:
	cd $(HTSDIR) && make CFLAGS=$(OSET)
	cd $(LOCAL_DIR)

install:
	cp $(MAIN) $(PREFIX)/$(MAIN)
	cp $(LOCAL_DIR)/bcred.pl $(PREFIX)/bcred
	chmod +x $(PREFIX)/bcred

clean:
	cd $(HTSDIR) && make clean
	rm -f $(LOCAL_DIR)/$(MAIN) *.o
	rm -rf $(LOCAL_DIR)/$(MAIN).dSYM
	rm -f $(LOCAL_DIR)/bcred
	
# DO NOT DELETE THIS LINE -- make depend needs it