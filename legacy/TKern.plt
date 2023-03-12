#! /usr/bin/gnuplot

set terminal tikz solid standalone header '\usepackage{amsmath,lmodern}' #sclae 0.67 0.67
set output "kerplot.tex"

set logscale;
set format '$10^{%L}$';
set yrange [1e-50:1e2];
set style data lines;
set key autotitle columnhead inside top right title '$qd$';
set key font ',7';
#
f(x) = a*x**(-5/2.);
a = 3.86447;
set ylabel '$M^{S}_q[t]$';
set xlabel '$t$';
plot [1e-8:1e24] for [col=2:22] 'TKernel.dat' u 1:col w l, \
		 f(x) ti '$2.56E(-4)*t^{-5/2}$';

