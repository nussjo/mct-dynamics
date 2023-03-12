#! /usr/bin/gnuplot

set logscale x;
set logscale y;
set format x '10^{%L}';
set style data lines;
#
set title 'Mean-Squared Displacement';
set ylabel '$MSD$';
set xlabel 't';
#
plot "MSD.dat" u 1:2 w l,\
	"../jn/tmp.msd" u 1:2 w l,\
	6./160*x ti '6Dt';
