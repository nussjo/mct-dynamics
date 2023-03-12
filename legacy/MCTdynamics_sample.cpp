// Correlator.cpp : 
//
// --------------------------------------------------------------------------------------
//
// AUTHORS :
// NAME 		- EMAIL				- ABBREVIATION
// Jonas Nußdorfer	- jonas.nussdorfer@gmail.com	- jn
//
// --------------------------------------------------------------------------------------
//
// CHANGELOG :
// 2016-10-17 : implemented solution of MCT problem with decimation technique; it seems
// 		that there is something wrong with the calculation of the memory kernel
// 		since most of the correlators evaluate to garbage NaN (jn)
// 2016-11-10 : program works for F12 model, but still needs a fix for the real memory
// 		kernel (jn)
// 2016-12-15 : program works for real memory kernel but only for equidistant grids and
//		adjusted critical packing fraction (jn)
// 2017-04-04 : program works for collective dynamics and a backup is found in Version 3.3
// 		the next step contains creating a solution for the tagged particle
// 		correlator (jn)
//
// --------------------------------------------------------------------------------------
//
// Usage : ./MCTdynamics <packing_fraction> <tau> <delta>
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - implement Filon-Fourier transform
// - implement MSD and VACF/memory kernel from the paper "Long-Time Tails in the VACF
//   for Brownian Dynamics"
//
// ======================================================================================


#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Correlator.h"
#include "Mesh.h"
#include "mdefs.h"

int main(int argc, const char **argv)
{
	const int EQUIDISTANT_POINTS = 50;	// number of grid points for equidistant part
	const int TIME_STEPS = 256;		// number of time steps per decimation block
	const int DECIMATION_STEPS = 100;	// number of decimation steps
	const double Q_MAX = 40;		// highest wavenumber in grid
	const double H_INITIAL = 1e-9;		// initial spacing on temporal domain
	const double PRECISION = 1e-15;		// calculation precision
	const int MAX_ITERATIONS = 10000;	// maximum number of iterations for solution
	const int SKIP_VALUES = 16;		// number of values to skip when outputting data to file
	const int GRID_ITERATIONS = 20;		// number of divisions for wavevector grid, determines smallest wavenumber

	const double V_THERM = 2.5;
	const double NUVAL = 1000;

	// start bisection
	double deltaphi = .7;
	double phimin = 0.2, phimax = 0.6;
	double phiold = phimin;
	double packing_fraction = 1.;
	int runs = 0;

	while (deltaphi>1e-11)
	{

	//double crit = 0.5114802; // pretty close

	packing_fraction = .5*(phimin+phimax);
	deltaphi = fabs(packing_fraction-phiold);
	std::cout << "\n\n# Number of runs = " << ++runs << "\n";
	std::cout << "# Current phi = " << packing_fraction << "\told phi = " << phiold << "\tdelta phi = " << deltaphi << "\n";

	// uncomment after bisection
	//double packing_fraction = (argc==2) ? atof(argv[1]) : crit;
	double tau = (argc==3) ? atof(argv[2]) : 0;
	double Delta = (argc==4) ? atof(argv[3]) : 0.1;

	Mesh Q(Q_MAX, EQUIDISTANT_POINTS, GRID_ITERATIONS, packing_fraction, tau, Delta);

	Collective C(TIME_STEPS, Q, H_INITIAL, V_THERM, NUVAL);
	C.initializeCorrelator();

	Tagged T(TIME_STEPS, Q, H_INITIAL, V_THERM, NUVAL, C);
	T.initializeCorrelator();

	Correlator *Corr[] = { &C, &T };

	MSD MSD(Q, C, T);
	MSD.initialize();

	std::vector<int> qout(Q.getSize(),0);
	for (int i=0; i<Q.getSize(); ++i)
	{
		qout[i] = i;
	}

	Corr[0]->openFiles(DECIMATION_STEPS, V_THERM, "CCorr.dat", "CKernel.dat");
	Corr[1]->openFiles(DECIMATION_STEPS, V_THERM, "TCorr.dat", "TKernel.dat");

	Corr[0]->writeColumnNames(qout);
	Corr[1]->writeColumnNames(qout);

	MSD.openFile("MSD.dat");

	// temporary values
	const array2D &tmp = C.getArrayValues();
	double lastvalue = 1;

	for (int d=1; d<=DECIMATION_STEPS; ++d)
	{	
		std::cout << "# Decimation step D = " << d << "/" << DECIMATION_STEPS << "  ";

		Corr[0]->calcSolution(PRECISION, MAX_ITERATIONS);
		Corr[1]->calcSolution(PRECISION, MAX_ITERATIONS);
		
		MSD.calculateMSD();
		MSD.calculateVACF();

		Corr[0]->writeDatatoFile(qout, SKIP_VALUES);
		Corr[1]->writeDatatoFile(qout, SKIP_VALUES);

		MSD.writeDatatoFile(SKIP_VALUES);

		lastvalue = tmp.back().back();

		Corr[0]->decimate();
		Corr[1]->decimate();
		
		MSD.decimate();

		std::cout << "\r";
		std::cout.flush();
	}

	
	
	std::cout << "\n";

	// end bisection
	if (lastvalue!=0) { phimax = packing_fraction; }
	else if (lastvalue==0) { phimin = packing_fraction; }
	phiold = packing_fraction;
			
	}

	// result

	std::cout << "\n\nfinal phi = " << packing_fraction << "!\n";

	return 0;
}
