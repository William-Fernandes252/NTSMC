CFLAGS =-std=c99 -pedantic -Wall -Werror
GDBFLAGS = -g
CC =gcc
CSOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
LDFLAGS = -lm 

slabs: $(CSOURCES:.c=.o)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

debug_slabs: $(CSOURCES:.c=.o)
	$(CC) $(GDBFLAGS) -o $@ $^ $(CFLAGS) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) -c $< $(CFLAGS) -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm $@ -f *.o *.c~ *.o~