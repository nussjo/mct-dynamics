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
	const int EQUIDISTANT_POINTS = 100;	// number of grid points for equidistant part
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

	//double crit = 0.515805; // pretty close for 500 equidistant gridpoints
	double crit = 0.5155251; // pretty close for 100 equidistant gridpoints

	//double packing_fraction = (argc==2) ? crit+atof(argv[1]) : crit;
	double tau = (argc==3) ? atof(argv[2]) : 0;
	double Delta = (argc==4) ? atof(argv[3]) : 0.1;

for (int j=0; j<30; ++j)
{
double a=0;
std::string close = "";

long long diff = j-14;

if (j<14 && j>=0)
{
	a = -pow(10,-(14-j)/2.);
	close = "L"+std::to_string(-diff);
}
else if (j<29 && j>=15 )
{
	a = pow(10,-(j-14)/2.);
	close = "G"+std::to_string(diff);
}
else
{
	close = "crit";
}

double packing_fraction = crit*(1+a);

std::cout << "# Current phi: " << packing_fraction << "\n";

	Mesh Q(Q_MAX, EQUIDISTANT_POINTS, GRID_ITERATIONS, packing_fraction, tau, Delta);
	Q.openFileandWrite("Mesh_"+close+".dat");

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

	Corr[0]->openFiles(DECIMATION_STEPS, V_THERM, "CCorr_"+close+".dat", "CKernel_"+close+".dat");
	Corr[1]->openFiles(DECIMATION_STEPS, V_THERM, "TCorr_"+close+".dat", "TKernel_"+close+".dat");

	Corr[0]->writeColumnNames(qout);
	Corr[1]->writeColumnNames(qout);

	MSD.openFile("MSD_"+close+".dat");

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

		Corr[0]->decimate();
		Corr[1]->decimate();
		
		MSD.decimate();

		std::cout << "\r";
		std::cout.flush();
	}
	
	std::cout << "\n";
}
	return 0;
}
