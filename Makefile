TARGET = euler1d

CC = gcc
CFLAGS = -O3 -Wall
DEBUGFLAGS = -DDEBUG=1 -Wall -ggdb

OBJ =  euler1d.o fluxes.o reconstruction.o  thermodynamics.o iniconditions.o boundaries.o ioutils.o varutils.o 



program: $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ)

debug: $(OBJ)
	$(CC) $(DEBUGFLAGS) -o $(TARGET) $(OBJ)	


%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	-rm -f *.o
	-rm -f $(TARGET)

