# Build flags.
DEBUG   =
#DEBUG	= -pg
#DEBUG	= -g
BFLAGS  = -fopenmp
CFLAGS	= -Ofast -std=c++11 -Wall -Wextra $(BFLAGS) $(DEBUG)
LFLAGS	= $(BFLAGS) $(DEBUG)

# Build tools.
CC	= g++
LD	= g++
COMP	= $(CC) $(CFLAGS) -c
LINK	= $(LD) $(LFLAGS) -o
AR      = ar