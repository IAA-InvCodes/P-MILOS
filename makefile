CC=mpicc

CFLAGS=
CFLAGS+=-O3 

ifdef develop
	ifeq ($(develop),yes)
		CFLAGS+=-g 
		CFLAGS+=-Wall -Wextra		
	endif
endif
ifdef use_double
	ifeq ($(use_double),yes)
		CFLAGS+=-D USE_DOUBLE_PRECISION=double
	endif
endif
CFLAGS+=-fno-omit-frame-pointer
#CFLAGS+=-Wall -Wextra
#CFLAGS+=-Wconversion
#CFLAGS+=-Wno-unused-but-set-variable -Wno-unused-parameter

SRCDIR= src
DEPEN=$(SRCDIR)/calculosCompartidos.o $(SRCDIR)/fgauss.o $(SRCDIR)/fvoigt.o  $(SRCDIR)/milos.o $(SRCDIR)/me_der.o $(SRCDIR)/mil_sinrf.o $(SRCDIR)/lib.o $(SRCDIR)/create_cuantic.o $(SRCDIR)/utilsFits.o $(SRCDIR)/milosUtils.o $(SRCDIR)/convolution.o $(SRCDIR)/readConfig.o
LDLIBS= -lm -lcfitsio -lnsl -lgsl -lgslcblas -lfftw3 -ldl -lpthread # -shared-intel
BIN= milos milosMPI 


all: $(BIN)

milos: $(DEPEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

milosMPI: $(DEPEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

clean:
	rm -f  $(SRCDIR)/*.o $(BIN)
