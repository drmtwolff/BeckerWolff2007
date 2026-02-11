#
#  Makefile for XSPEC model software subdirectory.
#
CC = gcc
FC = gfortran
LD = gfortran
FLD = gfortran
CP = cp -fv
RM = rm -fv
#
#  The variables below need to be changed for each local
#  installation of these codes.
#
#  Location of any extra libraries I need.
#
ifeq ($(OSTYPE),darwin)
	USRLIB = -L/opt/local/lib
endif

ifeq ($(OSTYPE),linux)
	USRLIB = 
endif

#
#  Location of any extra include files.
#
ifeq ($(OSTYPE),darwin)
	USRINC = -I/opt/local/include
endif

ifeq ($(OSTYPE),linux)
	USRINC = 
endif
#
#  Set up Make vaiables for default libraries I will need to 
#  link with:
#
LDLIBS = -lX11 -ltermcap -lm 
#
#  Now set the compiler and linker flags:
#
CFLAGS = -g -Wall $(USRINC)

LDFLAGS = $(USRLIB)

#
#  Set up programs to be compiled here:
#
CPROGS = phys2sim sim2phys 

#
#  Get going here:
#

objdir = ./obj
bindir = ./bin
srcdir = .

vpath   %.c     $(srcdir)
vpath   %.o     $(objdir)
VPATH   =       $(bindir)

$(objdir)/%.o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.o : %.c
	$(CC) -c $(CFLAGS) -o $(objdir)/$@ $<

.PHONY: all clean dirs

all: $(CPROGS)

dirs:
	mkdir -p  obj  bin

clean:
	$(RM) core
	- ( cd $(objdir) ; $(RM) *.o )
	- ( cd $(bindir) ; $(RM) $(CPROGS) )

phys2sim: $(objdir)/phys2sim.o
	$(LD) $(LDFLAGS) -o $(bindir)/$@ $^ $(LDLIBS)

sim2phys: $(objdir)/sim2phys.o
	$(LD) $(LDFLAGS) -o $(bindir)/$@ $^ $(LDLIBS)

