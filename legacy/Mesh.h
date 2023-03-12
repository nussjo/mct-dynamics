// Mesh.h :
// Class declaration which contains the wavevector grid, as well as the static structure
// factor and the direct correlation function. Also contains data on the vertices.
//
// --------------------------------------------------------------------------------------
//
// AUTHORS :
// NAME 		- EMAIL				- ABBREVIATION
// Matthias Sperl 	- matthias.sperl@dlr.de		- ms
// Jonas Nußdorfer	- jonas.nussdorfer@gmail.com	- jn
//
// --------------------------------------------------------------------------------------
//
// CHANGELOG :
// 2016-11-30 : current implementation works for equidistant wavevector grids; vertex
// 		integration/calculation still needs to be figured out for non-equidistant
// 		grids in terms of integration spacings (jn)
// 2017-04-07 : implementation of non-equidistant grid and both vertices functions done (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - generalize grid creation to a function call; if possible with the use of a function
//   pointer to guarantee usability with different types of grids (jn)
//
// ======================================================================================

#ifndef MESH_H
#define MESH_H

#include <fstream>
#include <vector>
#include <string>

// Typedef definition for multidimensional arrays with STL vectors
typedef std::vector<std::vector<std::vector<double>>> array3D;
typedef std::vector<std::vector<double>> array2D;

class Mesh
{
private:
	// Member variables :
	double m_eta;			// packing fraction of the system
	double m_tau;			// Tau value for calculation of direct correlation function
	double m_delta;			// Delta value for calculation of direct correlation function
	double m_particleDensity;	// particle density
	int m_size;			// size of the grid
	std::vector<double> m_qgrid;	// contains all wavenumbers of the grid stored in ascending order
	std::vector<double> m_Sq;	// corresponding values of the static structure factor
	std::vector<double> m_cq;	// corresponding values of the direct correlation function
	array3D m_CVertex;		// contains the values of the collective vertex
	array3D m_SVertex;		// contains the values of the incoherent vertex
	array2D m_BorderIndices;	// index lookup-table for floor/ceil indices of |q-k| / q+k respectively
	std::ofstream m_file;		// filestream for the file which is used to store data

	// Grid generation functions :
	void sizeUpdate();
	void generateGrid(double, int, int);

protected:
	// Internal initialization functions :
	void determineIntegralBounds();		// evaluate |q-k| / q+k for all possible combinations and store floor/ceil index in array
	void initializeVertex();		// initializes both vertex functions

	// Internal calculation functions for member arrays :
	double calculateDirectCorrelation(double);	// calculates the direct correlation function value-wise from an analytic expression
	double calculateStaticStructure(double);	// calculates the static structure faction value-wise from the direct correlation function using the Ornstein-Zernicke relation
	double calcCollectiveVertex(int, int, int);	// calculate element on 3D-array of the collective vertex with all 3 wavevectors on the grid (index as input)
	double calcCollectiveVertex(int, int, double);	// calculate element on 3D-array of the collective vertex with the last wavevector not on the grid
	double calcTaggedVertex(int, int, int);		// calculate element on 3D-array of the incoherent vertex with all 3 wavevectors on the grid (index as input)
	double calcTaggedVertex(int, int, double);	// calculate element on 3D-array of the incoherent vertex with the last wavevector not on the grid

public:
	// Constructor and Destructor
	Mesh(double, int, int, double, double, double);
	Mesh() = delete;
	virtual ~Mesh();

	// Getter-functions for member variables :
	double qback() const;			// last q-grid value
	double qfront() const;			// first q-grid value
	double getq(const int) const;		// q-grid value at given index
	double getSq(const int) const;		// value of static structure function at given index
	double getcq(const int) const;		// value of direct correlation function at given index
	double getPackingFraction() const;	// packing fraction of this system
	double getTau() const;			// Tau value in direct correlation function calculation
	double getDelta() const;		// Delta value in direct correlation function calculation
	double getParticleDensity() const;	// returns particle density of the system
	int getSize() const;			// size of the q-grid

	// Data->File :
	void writeDatatoFile();							// writes all data to file after calculation (except vertex data)
	void writeDatatoFile(std::vector<double>&, std::vector<double>&);	// writes data to file including frequencies from correlator
	void openFile(std::string);						// open file to write data

	// Reference functions to return references to member arrays :
	// returned as const reference (read-only)!
	const array3D& getCollectiveVertex();	// vertex for collective dynamics
	const array3D& getTaggedVertex();	// vertex for incoherent dynamics
	const array3D& getVertex(bool);		// choose which vertex to return
	const array2D& getBorderIndices();	// returns array with indices marking the floor/ceil index of |q-k| and q+k respectively
	const std::vector<double>& getQgrid();	// returns the wavevector grid
	const std::vector<double>& getcq();	// returns grid containing the values of the direct correlation function
	const std::vector<double>& getSq();	// returns grid containing the values of the static structure function
		
};

// Overloaded << operator to easily write q/Cq/Sq-data to a stream (console/file)
// in one step
std::ostream& operator << (std::ostream&, const Mesh&);

#endif
