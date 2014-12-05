CC = gcc
CFLAGS = -g -lm -Wall -std=c99
OBJECTS = precession.o
all : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o precession


precession: .c
	$(CC) $(CFLAGS) -c precession.c
