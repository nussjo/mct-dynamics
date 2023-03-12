// Utils.h :
// Class declaration which contains the wavevector grid, as well as the static structure
// factor and the direct correlation function. Also contains data on the vertices.
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
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
//
// ======================================================================================

#include "Utils.h"

double getInterpolationValue(double x, double leftIP, double leftIV, double rightIP, double rightIV)
{
	//return leftIV + (x-leftIP)*(rightIV-leftIV)/(rightIP-leftIP);
	double f = (x-leftIP)/(rightIP-leftIP);

	if (leftIV<=0 || rightIV<=0)
	{
		return 0;
	}
	else
	{
		return pow(fabs(leftIV),1-f)*pow(fabs(rightIV),f);
	}
}
