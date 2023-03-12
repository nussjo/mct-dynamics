// constants.h : this header file includes all necessary calculation parameters which do
//		 NOT change during calculation ---> all variables here will be of type
// 		 "const <variable type>"
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
// 2016-08-23 : created namespace grid_sizes with size of wavevector grid, starting-
// 		and endpoint of the wavevector grid, array sizes for calculations and
// 		the number of decimations (jn)
// 2016-08-24 : created namespace system_parameters; contains physical parameters &
// 		constants for the calculation (jn)
//
// ======================================================================================

#ifndef CONSTANTS_H
#define CONSTANTS_H


namespace computation
{
	const int EQUIDISTANT_POINTS = 100;	// number of grid points for equidistant part
	const int TIME_STEPS = 256;		// number of time steps per decimation block
	const int NUMBER_OF_DECIMATIONS = 100;	// number of decimation steps
	const double Q_LOWER = 4;		// lowest wavenumber in grid --- obsolete!!!
	const double Q_UPPER = 40;		// highest wavenumber in grid
	const double H_INITIAL = 1e-9;		// initial spacing on temporal domain
	const double PRECISION = 1e-15;		// calculation precision
	const int MAX_ITERATIONS = 10000;	// maximum number of iterations for solution
	const int SKIP_VALUES = 16;		// number of values to skip when outputting data to file
	const int GRID_ITERATIONS = 20;
}


namespace physics
{
	const double V_THERM = 2.5;
	const double NUVAL = 1000;
	const double DIFFUSION_CONSTANT = V_THERM*V_THERM/NUVAL;
}


#endif
