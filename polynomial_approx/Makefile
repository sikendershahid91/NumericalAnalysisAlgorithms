C=gcc
CFLAGS=-Wall -c -g -lm -pthread
LFLAGS =-Wall -g -lm -pthread
OBJS=roots.o

all: roots

roots: roots.o
	$(CC) $(LFLAGS) $(OBJS) -o roots -lm

roots.o: roots.c
	$(CC) $(CFLAGS) roots.c

run1:
	./roots -7.007359873672922568e+00 -6.434467526661742376e+00 -100 100

run2:
	./roots -3.482446274320096791e-01 4.668557702739741799e+00 -9.027222188370192413e+00 5 1000

run3: 
	./roots 1.0e+00 -2.0e+00 0.0e+00 1.0e+00 5.0e+00 0 2

run4:
	./roots 1 -3 0 2 -2 0

run5:
	./roots 1 0.66667 -1 -7 -10 20

.PHONY: clean
clean:
	-rm *.o roots