#! /usr/bin/gnuplot

set logscale;
set format '10^{%L}';
set style data lines;
#
set title 'Velocity Autocorrelation function';
set ylabel '$VACF$';
set xlabel 't';
#
plot "MSD.dat" u 1:(abs($6)) w l, x**(-2.5) lc 'blue';
