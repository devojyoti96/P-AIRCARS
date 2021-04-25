DESTDIR = $(PWD)/
CC=gcc
GSL_INCLUDE_DIR=/data1/devojyoti/PhD/ankflag/gsl/include/
GSL_LIBRARIES=-L/data1/devojyoti/PhD/ankflag/gsl/lib -Wl,"-R /data1/devojyoti/PhD/ankflag/gsl/lib"
CLIBS= $(GSL_LIBRARIES) -lgsl -lgslcblas -lm 
CFLAGS=-I. -I./include -I$(GSL_INCLUDE_DIR) -g -fopenmp
c = $(CC) $(CFLAGS) 
o = -o$(DESTDIR)$@
.c.o: ; $c -c  $< 
.f.o: ; $f -c $< 
.c: ; $c $o $< $(CLIBS)
.f: ; $f $o $< 
make: ankflag
ANKOBJ=ankflg.c exc_fns.c flg_fns.c basic_fns.c aux_fns.c init_pars.c
ankflag:$(ANKOBJ)
	$c $(ANKOBJ) $o $(CLIBS)
clean:
	rm ankflag
