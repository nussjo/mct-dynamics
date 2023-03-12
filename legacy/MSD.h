// MSD.h :
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

#ifndef MSD_H
#define MSD_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Mesh.h"
#include "Correlator.h"

// Typedef definition for multidimensional arrays with STL vectors
typedef std::vector<std::vector<double>> array2D;

class MSD
{
protected:
	// Member variables :
	Mesh &m_qmesh;				// reference to a Mesh class object (contains all wavevector relevant data)
	int m_tsteps;				// size of one solution block
	double m_spacing;			// initial spacing in temporal domain
	double m_D0;				// short-time diffusion constant
	std::vector<double> m_time;		// array for time grid
	std::vector<double> m_MSD, m_mMSD;	// Mean-Squared Displacement and its memory kernel
	std::vector<double> m_dMSD, m_dmMSD;	// differences of array values above
	const array2D &m_phiA, &m_phiB;		// references to arrays for values of collective and tagged particle correlators
	std::string m_filename;			// filename for the MSD data
	std::ofstream m_file;			// files for the correlators and kernels
		
public:
	// Constructor & Destructor :
	MSD(Mesh&, const Correlator&, const Correlator&, std::string);
	MSD() = delete;
	virtual ~MSD();

	// Getter-functions for member variables :
	int getTimesteps();	// returns the size of time block
	double getSpacing();	// return current time spacing

	// Data->File functionality :
	void writeFileHeader();		// writes a header to the file which contains all relevant data used
	void writeColumnNames();	// writes column names to make distinctions between different data easier
	void writeDatatoFile(int);	// actual data->file function

	// Reference functions to return references to member arrays :
	// returned as const reference (read-only)!
	const std::vector<double>& getMSDValues();			// returns the 2D-array with correlator values
	const std::vector<double>& getKernelValues();			// returns the 2D-array with the kernel values

	// Main calculation functionality :
	void decimate();				// perform a decimation procedure on the current arrays

	// Calculation of the MSD :
	double calculateMSDKernel(int);	// calculate memory kernel of MSD at given time
	void calculateMSD();		// calculate time evolution of MSD
};


#endif
