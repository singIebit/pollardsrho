CC = gcc
CFLAGS = -g -Wall -pthread
LDFLAGS = -lm -lgmp

SRCS = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJS = $(SRCS:.c=.o)
EXEC = pollardsrho

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
