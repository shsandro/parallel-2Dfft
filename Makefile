CC=gcc
MPI=mpicc
CFLAGS=-lm
SERIAL_BIN=seq_2dfft.out
PARALLEL_BIN=par_2dfft.out
SRC=serial_2dfft.c complex.c
INCLUDES=./includes/complex.h
OBJECTS=serial_2dfft.o complex.o

all: serial parallel

serial: $(SERIAL_BIN)

$(SERIAL_BIN): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS)

$(OBJECTS): %: $(SRC)
	$(CC) -c -Wall -o serial_2dfft.o serial_2dfft.c
	$(CC) -c -Wall -o complex.o complex.c

$(SRC): $(INCLUDES)

parallel: $(PARALLEL_BIN)

$(PARALLEL_BIN): par_2dfft.c complex.c
	$(MPI) -o $@ $^ $(CFLAGS)

complex.c: $(INCLUDES)

clean:
	rm -rf *.o *.out