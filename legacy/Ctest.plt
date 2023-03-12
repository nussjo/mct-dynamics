#! /usr/bin/gnuplot

set logscale x;
set logscale y;
set yrange [1e-40:1e2];
set format x '10^{%L}';
set style data lines;
set key autotitle columnhead inside top right title '$q$';
#
f(x)=a*x**(-3/2.);
a = 2e-4;
set fit quiet;
fit [1e5:1e10] f(x) 'CKernel.dat' u 1:2 via a;
set style data lines;
D = 1./160.;
#
set title 'Collective';
set ylabel '$m_q$';
set xlabel 't';
ares = a*(D**(3./2.));
pf = 2.48591e-4*(D**(3./2.));
plot [1e-8:1e24] for [col=2:22:5] 'CKernel.dat' u ($1*D):col w l notitle, \
				    f(x)*(D**(3./2.)) ti sprintf("%.6e*t^(-3/2)",ares), \
				    pf*x**(-3./2) ti sprintf("%.6e*t^(-3/2)",pf);
#
