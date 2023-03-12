#! /usr/bin/gnuplot

set logscale x;
set logscale y;
set yrange [1e-40:1e3];
set format x '10^{%L}';
set style data lines;
set key autotitle columnhead inside top right title '$q$';
#
f(x)=a*x**(-3/2.);
a = 2e-4;
set fit quiet;
fit [2e7:1e13] f(x) 'CKernel.dat' u 1:2 via a;
#
set title 'Collective';
set ylabel '$m_q$';
set xlabel 't';
plot [1e-8:1e24] for [col=2:800:80] 'CKernel.dat' u 1:col w l notitle, \
				    f(x) ti sprintf("%.6e*t^(-3/2)",a);
#
