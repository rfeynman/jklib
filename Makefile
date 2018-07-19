
OBJSO := spline_class.o  mxequ.o \
dminv.o  rampFit.o  complex_class.o jams.o \
rfn3.o microsleep.o poly.o amoeba_class.o powell_class.o \
rungeKutta_class.o fourier_class.o bessel.o fish_class.o Faddeeva.o \
grace_class.o xmgr_class.o

CC= clang

all:	$(OBJSO) libjkLib.a  


libjkLib.a:	$(OBJSO)
	ar -cr libjkLib.a $(OBJSO)

minim_class.o:		minim_class.cxx minim_class.hxx
	$(CC) minim_class.cxx -c -g -O0 
amoeba_class.o:		amoeba_class.cxx amoeba_class.hxx
	$(CC) amoeba_class.cxx -c -g -O0 
complex_class.o:	complex_class.cxx complex_class.hxx
	$(CC) complex_class.cxx -c -g -O0 
dminv.o:		dminv.cxx
	$(CC) dminv.cxx -c -g -O0 
jams.o:			jams.cxx jams.hxx
	$(CC) jams.cxx -c -g -O0 
knobcontainer.o:	knobcontainer.cxx knobcontainer.hxx
	$(CC) knobcontainer.cxx -c -g -O0 
microsleep.o:		microsleep.cxx
	$(CC) microsleep.cxx -c -g -O0 
mxequ.o:		mxequ.cxx 
	$(CC) mxequ.cxx -c -g -O0 
powell_class.o:		powell_class.cxx powell_class.hxx
	$(CC) powell_class.cxx -c -g -O0 
poly.o:			poly.cxx poly.hxx
	$(CC) poly.cxx -c -g -O0 
rampFit.o:		rampFit.cxx rampFit.hxx
	$(CC) rampFit.cxx -c -g -O0 
rampcontainer.o:	rampcontainer.cxx rampcontainer.hxx
	$(CC) rampcontainer.cxx -c -g -O0 
rfn3.o:			rfn3.cxx rfn3.hxx
	$(CC) rfn3.cxx -c -g -O0 
grace_class.o:		grace_class.cxx grace_class.hxx
	$(CC) grace_class.cxx -c -g -O0 
xmgr_class.o:		xmgr_class.cxx xmgr_class.hxx
	$(CC) xmgr_class.cxx -c -g -O0 
fish_class.o:		fish_class.cxx fish_class.hxx
	$(CC) fish_class.cxx -c -g -O0 
spinor_class.o:		spinor_class.cxx spinor_class.hxx
	$(CC) spinor_class.cxx -c -g -O0 
spline_class.o:		spline_class.cxx spline_class.hxx
	$(CC) spline_class.cxx -c -g -O0 
stonelist.o:		stonelist.cxx stonelist.hxx
	$(CC) stonelist.cxx -c -g -O0 
wireup.o:		wireup.cxx wireup.hxx
	$(CC) wireup.cxx -c -g -O0 
rungeKutta_class.o:		rungeKutta_class.cxx rungeKutta_class.hxx
	$(CC) rungeKutta_class.cxx -c -g -O0 
fourier_class.o:		 fourier_class.cxx  fourier_class.hxx
	$(CC)  fourier_class.cxx -c -g -O0 
bessel.o:		 bessel.cxx  bessel.hxx
	$(CC)  bessel.cxx -c -g -O0 

Faddeeva.o:	Faddeeva.cxx Faddeeva.hxx
	$(CC) Faddeeva.cxx -c -g -O0

clean:	
	rm -f *.o  libjkLib.a
