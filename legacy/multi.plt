#! /usr/bin/gnuplot

set multiplot layout 1,2;
set xtics rotate;
set bmargin 5;
#
set logscale x;
set yrange [0:1];
set format x '10^{%L}';
set style data lines;
set nokey;
#
set title 'Collective';
set ylabel '$phi_q$';
set xlabel 't';
plot [1e-8:1e24] for [col=2:120:10] 'CCorr.dat' u 1:col;
#
set title 'Incoherent';
set ylabel '$phis_q$';
set xlabel 't';
plot [1e-8:1e24] for [col=2:120:10] 'TCorr.dat' u 1:col;
#
unset multiplot;
