CC = gcc
CFLAGS = -Wall -Werror
#si erreur de compilation enlever Werror

all: 2048.x

movement.o : movement.c movement.h
	$(CC) $(CFLAGS) -c movement.c -o movement.o

functions.o : functions.c functions.h movement.c movement.h
	$(CC) $(CFLAGS) -c functions.c -o functions.o

2048-template.o :  2048-template.c functions.c functions.h movement.c movement.h
	$(CC) $(CFLAGS) -c 2048-template.c -o 2048-template.o

2048.x : movement.o functions.o 2048-template.o
	$(CC) $(CFLAGS) -o 2048.x movement.o functions.o 2048-template.o -lm

clean:
	rm -f *.o 

realclean: clean
	rm -f *.x

run : 2048.x
	./2048.x