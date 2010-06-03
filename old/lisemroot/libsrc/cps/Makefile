# makefile for cygnus gcc, and gnu-make

# library specific part
# name of library
LIB_TARGET=cps
# MAKE_PUBLIC: include files to be linked in by
# lisem
MAKE_PUBLIC=cps.h mprolog.h

OBJS = \
   lddutil.o \
   chkdata.o \
   upstream.o \
   cpsparam.o\
   cpswarn.o\
   cpstest.o\
   cpslang.o\
   cpsio.o\
   cpsmem.o\
   initmaps.o\
   mask.o\
   compress.o\
   stdlibx.o \
   calcheig.o\
   cpstimin.o\
   timeoutp.o\
   kinemati.o\
   sprdldd.o

# PCRTREE debris:
#   chandisp.o\ not used
#   disperse.o\ not used
NO_EXAMPLES=TRUE

# standard part

CC=gcc
CFLAGS+= -ansi -pedantic
CFLAGS+= -DDEBUG -I../../include
# versions:
OPTIMIZE_FLAGS+=-O3 -fomit-frame-pointer -m486
DEBUG_FLAGS+= -g -DDEBUG_DEVELOP
# if building final version:
CFLAGS+=$(OPTIMIZE_FLAGS)
# if building debug version:
# CFLAGS+=$(DEBUG_FLAGS)

TARGET=lib$(LIB_TARGET).a

$(TARGET): $(OBJS)
	@ar rcs $@ $^

.c.o:
	$(CC) $(CFLAGS) -c $<

install: $(TARGET) $(MAKE_PUBLIC)
	cp -f $(TARGET) $(MAKE_PUBLIC) ../../include

clean:
	rm $(TARGET) $(OBJS)
