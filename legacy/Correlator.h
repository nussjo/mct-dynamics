// Correlator.h :
// Declaration of correlator class which contains arrays for a time block of a system of
// density correlators which are wavevector-dependent, as well as corresponding arrays
// for their respective memory kernels and physical quantities. The evolution of the
// solution is carried out with a decimation procedure.
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
// 2016-11-30 : current implementation works for equidistant wavevector grids and yields
// 		the same results as the old program by M. Sperl (jn)
// 2017-04-06 : split up Correlator class into a base part (general functionality) and
// 		derived classes (collective/incoherent dynamics) (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - Fourier transformation (jn)
//
// ======================================================================================

#ifndef CORRELATOR_H
#define CORRELATOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"

// Typedef definition for multidimensional arrays with STL vectors
typedef std::vector<std::vector<std::vector<double>>> array3D;
typedef std::vector<std::vector<double>> array2D;

class Decimation
{
protected:
	double m_spacing;		// initial spacing in temporal domain
	int m_tsteps;			// size of time block per decimation step T
	int m_qsteps;			// number of rows for 2D array M
	std::vector<double> m_time;	// array for time grid
	array2D m_values, m_kernel;	// MxT-sized 2D-array for the array and kernel values
	array2D m_dP, m_dM;		// MxT-sized 2D-array to store the differences of neighbouring values

public:
	// Constructor & Destructor :
	Decimation(int, int, double);
	Decimation() = delete;
	virtual ~Decimation();

	// Getter-functions for member variables :
	int getTimesteps() const;	// returns the size of time block
	int getQgridsteps() const;	// returns the number of correlators
	double getSpacing() const;	// return current time spacing

	// Reference functions to return references to member arrays :
	// returned as const reference (read-only)!
	const array2D& getArrayValues() const;	// returns the 2D-array with correlator values
	const array2D& getKernelValues() const;	// returns the 2D-array with the kernel values

	// Decimation functionality :
	virtual void decimate();	// perform a decimation procedure on the current arrays
};


class Correlator : public Decimation
{
protected:
	// Member variables :
	Mesh &m_qmesh;				// reference to a Mesh class object (contains all wavevector relevant data)
	double m_D0;				// short-time diffusion constant
	double m_nu;				// damping
	bool m_isNewtonian;			// boolean values to switch on/off newtonian/brownian dynamics
	std::vector<double> m_Omega;		// wavevector-dependent frequency
	std::vector<double> m_vertInt;		// array to temporarily store values of the inner memory kernel integral
	std::ofstream m_phifile, m_kernelfile;	// files for the correlators and kernels

	// Constructor & Destructor :
	Correlator(int, Mesh&, double, double, double);
	Correlator() = delete;
	virtual ~Correlator();
		
public:
	// Getter-functions for member variables :
	double getDiffusionConstant() const;	// return short-time diffusion constant
	double getFrictionCoefficient() const;	// return friction coefficient for given index
	bool isNewtonian() const;		// return bool whether dynamics are brownian or newtonian

	// Pure virtual functions :
	virtual double calculateLowQKernel(int) = 0;	// calculate the q->0 memory kernel at t
	virtual double calculateLowQPrefactor() = 0;	// calculate the long-time prefactor of the q->0 kernel
	virtual void initializeCorrelator() = 0;	// initialize correlators: use outside of constructor after object creation

	// Data->File functionality :
	void openFiles(int, double, std::string, std::string);	// opens files and calls function to write column names
	void writeFileHeader(int, double);			// writes a header to the file which contains all relevant data used
	virtual void writeColumnNames(std::vector<int>&);	// writes column names to make distinctions between different data easier
	virtual void writeDatatoFile(std::vector<int>&, int);	// actual data->file function

	// Reference functions to return references to member arrays :
	// returned as const reference (read-only)!
	const std::vector<double>& getLowQKernelValues() const;	// returns the 1D-array with the q->0 kernel values
	const std::vector<double>& getFrequencies() const;	// returns the 1D-array with the q->0 kernel values

	// Main calculation functionality :
	double calculateKernel(int, int, const array3D&, const array2D&, const array2D&);
	virtual void calcSolution(const double, const int);	// calculate solution
};


class Collective : public Correlator
{
public:
	// Constructor & Destructor :
	Collective(int, Mesh&, double, double, double);
	Collective() = delete;
	virtual ~Collective();

	// Declarations for pure virtual functions from base class :
	// See above! Implementation for collective correlator dynamics.
	virtual double calculateLowQKernel(int);
	virtual double calculateLowQPrefactor();
	virtual void initializeCorrelator();
};


class Tagged : public Correlator
{
private:
	// Member variables :
	double m_DL;
	const Correlator &m_collective;	// const reference to the values of the collective correlators

public:
	// Constructor & Destructor :
	Tagged(int, Mesh&, double, double, double, const Correlator&);
	Tagged() = delete;
	virtual ~Tagged();

	// Data->File functionality : 
	// Overrides for base class methods due to slight differences for tagged particles
	virtual void writeColumnNames(std::vector<int>&);	// writes column names to make distinctions between different data easier
	virtual void writeDatatoFile(std::vector<int>&, int);	// actual data->file function

	// Declarations for pure virtual functions from base class :
	// See above! Implementation for tagged particle correlator dynamics.
	virtual double calculateLowQKernel(int);
	virtual double calculateLowQPrefactor();
	virtual void initializeCorrelator();

	// virtual function overload :
	virtual void calcSolution(const double, const int);
};


class MSD : public Decimation
{
protected:
	// Member variables :
	Mesh &m_qmesh;				// reference to a Mesh class object (contains all wavevector relevant data)
	std::vector<double> m_diffusionConstant;// short-time diffusion constant
	double m_nu;				// friction coefficient
	bool m_isNewtonian;			// boolean variable indicating newtonian or brownian dynamics
	const array2D &m_phiA, &m_phiB;		// references to arrays for values of collective and tagged particle correlators
	const std::vector<double> &m_kernelB;	// q->0 memory kernel of self motion correlator
	std::ofstream m_file;			// files for the correlators and kernels
		
public:
	// Constructor & Destructor :
	MSD(Mesh&, const Collective&, const Tagged&);
	MSD() = delete;
	virtual ~MSD();

	// Initialization :
	void initialize();	// initialize MSD separately

	// Data->File functionality :
	void openFile(std::string);	// open file
	void writeFileHeader();		// writes a header to the file which contains all relevant data used
	void writeColumnNames();	// writes column names to make distinctions between different data easier
	void writeDatatoFile(int);	// actual data->file function

	// Calculation of the MSD :
	double calculateMSDKernel(int);	// calculate memory kernel of MSD at given time
	void calculateMSD();		// calculate time evolution of MSD
	void calculateVACF();		// calculate time evolution of VACF
	void calculateDiffusionConstant();

	// Getter-functions for diffusion constant :
	const std::vector<double>& getDiffusionConstant();
	double getDiffusionConstant(int);
};

#endif
