# Build flags.
DEBUG	=
#DEBUG	= -g
CFLAGS	= -Ofast -std=c++11 -Wall -Wextra $(DEBUG)
LFLAGS	=  $(DEBUG)

# Build tools.
CC	= g++
LD	= g++
COMP	= $(CC) $(CFLAGS) -c
LINK	= $(LD) $(LFLAGS) -o
AR      = ar

