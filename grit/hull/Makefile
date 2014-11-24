# to install executable into $(BINDIR),
# and library into $(LIBDIR),
#       type "make".

CC	= gcc
#CFLAGS	= -fullwarn -g
OBJS	= hull.o ch.o io.o rand.o pointops.o fg.o
HDRS	= hull.h points.h pointsites.h stormacs.h
SRC	= hull.c ch.c io.c rand.c pointops.c fg.c
LIB	= lib$(PROG).a
PROG	= hull

all	: $(PROG)

# the following probably needn't be changed

$(OBJS) : $(HDRS)

hullmain.o	: $(HDRS)

$(PROG)	: $(OBJS) hullmain.o
	$(CC) $(CFLAGS) $(OBJS) hullmain.o -o $(PROG) -lm

clean	:
	-rm -f $(OBJS) hullmain.o core a.out $(PROG) hull.ilk rsites.ilk rsites.suo hull.suo rsites$(EXE)

