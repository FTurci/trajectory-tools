# Build flags.
DEBUG	=
#DEBUG	= -g
CFLAGS	= -Ofast -std=c++0x -Wall -Wextra $(DEBUG)
LFLAGS	= $(DEBUG)

# Build tools.
CC	= g++
LD	= g++
COMP	= $(CC) $(CFLAGS) -c
LINK	= $(LD) $(LFLAGS) -o
AR      = ar

