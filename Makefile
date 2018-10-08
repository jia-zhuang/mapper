CC=gcc
CFLAGS=-W -Wall -O2 -std=c99
LDFLAGS=-lz -lrt

mapper: main.o utils.o utils.h index.o index.h align.o align.h shm.o shm.h
	$(CC) -o mapper main.o utils.o index.o align.o shm.o $(LDFLAGS)

test: utils.o
	$(CC) -o test utils.o $(LDFLAGS)

dump: dump_idx.o utils.o utils.h
	$(CC) -o dump utils.o dump_idx.o $(LDFLAGS)

clean:
	rm -f mapper
	rm -f *.o
