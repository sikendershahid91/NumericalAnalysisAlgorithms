CC = g++
CFLAGS=-Wall -c -g 
LFLAGS =-Wall -g
OBJS=linearsys.o

all: linearsys
	
linearsys: $(OBJS)
	$(CC) -std=c++11 $(LFLAGS) $(OBJS) -o linearsys

$(OBJS) : linearsys.cpp
	$(CC) -std=c++11 $(CFLAGS) linearsys.cpp


clean:
	rm outfile.txt