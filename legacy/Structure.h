// Mesh.h : class declaration which contains the wavevector grid, as well as the
// 	    static structure factor and the direct correlation function.
// 	    Implementation of the vertex matrices as described in the diploma thesis
//	    by M. Sperl ("Glass Transition in Colloids with Attractive Interaction",
// 	    p.84f) ---> only suitable for equidistant grids (jn)
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
// 2016-11-30 : current implementation works for equidistant wavevector grids; vertex
// 		integration/calculation still needs to be figured out for non-equidistant
// 		grids in terms of integration spacings (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - implementation for non-equidistant grids still missing
//
// ======================================================================================

#ifndef STRUCTURE_H
#define STRUCTURE_H


#include <fstream>
#include <vector>
#include <string>

typedef std::vector<std::vector<std::vector<double>>> array3D;
typedef std::vector<std::vector<double>> array2D;

class Structure
{
private:
	double m_eta, m_tau, m_delta;
	int m_size;
	std::vector<double> m_qgrid, m_Sq, m_cq;
	array3D m_CVertex, m_SVertex;
	array2D m_BorderIndices;
	std::string m_filename;
	std::ofstream m_file;	

public:
	Structure(double, int, int, double, double, double);
	virtual ~Structure();

	double qback() const;
	double qfront() const;
	double getq(const int) const;
	double getSq(const int) const;
	double getcq(const int) const;
	double getPackingFraction() const;
	double getTau() const;
	double getDelta() const;
	double getInterpolationValue(double, double, double, double, double);

	const array3D& getVertex();
	const array2D& getBorderIndices();
	const std::vector<double>& getQgrid();
	const std::vector<double>& getcq();
	const std::vector<double>& getSq();

	void writeDatatoFile();
	void determineIntegralBounds();
	void initializeVertex();

	double calculateDirectCorrelation(double);
	double calculateStaticStructure(double);

	int getSize() const;
	double calcVertex(int, int, int);
	double calcVertex(int, int, double);
};

std::ostream& operator << (std::ostream&, const Structure&);

#endif
