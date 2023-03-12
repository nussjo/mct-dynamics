#! /usr/bin/gnuplot

set logscale x;
set yrange [0:1];
set format x '10^{%L}';
set style data lines;
set nokey;
#
set ylabel '$phi_q$';
set xlabel 't';
plot [1e-8:1e24] for [col=2:200] 'TCorr.dat' u 1:col;
