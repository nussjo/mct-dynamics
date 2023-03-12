#! /usr/bin/gnuplot

set terminal tikz solid standalone header '\usepackage{amsmath,lmodern}' #sclae 0.67 0.67
set output "kerplot.tex"
#
set logscale x;
set logscale y;
set format '$10^{%L}$';
set yrange [1e-50:1e2];
set style data lines;
set key autotitle columnhead inside top right title '$q$';
set key font ',7';
#
f(x) = a*x**(-3/2.);
a = 2e-4;
set fit quiet;
fit [1e8:1e11] f(x) 'CKernel.dat' u 1:2 via a;
#
pf = 2.53578e-4;
#
set ylabel '$m_{q}[t]$';
set xlabel '$t$';
plot [1e-8:1e24] for [col=2:22:5] 'CKernel.dat'	u 1:col w l, \
				  'CObs.dat' u 1:2 w l ti '$m_{q\rightarrow0}$', \
		 		  f(x) ti sprintf('$%.6e*t^{-3/2}$',a), \
				  pf*x**(-3./2) ti sprintf('$%.6e*t^{-3/2}$',pf);

