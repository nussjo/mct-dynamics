#! /usr/bin/gnuplot

set logscale x;
set logscale y;
set yrange [1e-40:1e2];
set format x '10^{%L}';
set style data lines;
set key autotitle columnhead inside top right title '$q$';
#
f(x)=a*x**(-5/2.);
a = 1e-5;
set fit quiet;
fit [1e5:1e10] f(x) 'TKernel.dat' u 1:2 via a;
set style data lines;
#
set title 'Incoherent';
set ylabel '$q^{2} m_q$';
set xlabel 't';
D = 1./160.;
ares = a*(D**(5./2.));
pf = 20.78283*(D**(5./2.));
plot [1e-8:1e24] for [col=2:22:5] 'TKernel.dat' u (D*$1):col w lines notitle, \
				    f(x)*(D**(5./2.)) ti sprintf("%.6e*t^(-5/2)",ares), \
				    pf*x**(-5./2) ti sprintf("%.6e*t^(-5/2)",pf);
#
