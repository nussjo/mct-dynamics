// Mesh.cpp : 
// Definition of the class member functions used in "Mesh.h"; main functionality for
// structure/wavevector related quantities. Implementation of the vertex matrices as
// described in the diploma thesis by M. Sperl ("Glass Transition in Colloids with
// Attractive Interaction", p.84f) are only suitable for equidistant grids. Hence, a more
// direct approach was chosen here.
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
// 2016-10-17 : moved c_swf_first & s_swf_first (by ms) functions into the Mesh
//		class to calculate cq/Sq; slightly changed function implementation (jn)
// 2016-11-16 : changed vertices matrices from copied form to the form presented in the
// 		diploma thesis of M. Sperl (jn)
// 2017-03-23 : divide et impera-styled grid implemented and currently in use (jn)
// 2017-04-07 : implementation of non-equidistant grid and both vertices functions done (jn)
// 2017-05-09 : generalized grid generation to be used with a function pointer. This 
// 		allows to easier switch between different grids (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - adapt all functionality to be used with any grid (especially integration bound
//   indices) and try to out-source the grid creation to outside the class (jn)
//
// ======================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include "mdefs.h"
#include "Mesh.h"

// Constructor for Mesh class objects. This constructor initializes all member variables
// with the initializer list and opens a file "Mesh.dat". In this file, all data of the
// wavevector grid, direct correlation function and static structure factor will be
// stored. The file will be closed in the destructor.
// After member initialization, the wavevector grid has to be generated. Note that we
// still need to resize all member arrays to the grid size since in the current implementation
// the final grid size is not given.
// The final step consists of calculating the direct correlation function and static structure
// factor and write the data to the file, as well as an initialization of integration border
// indices and vertex matrices.
Mesh::Mesh(double max, int eqpoints, int geompoints, double eta, double tau, double delta) :
	m_eta(eta), m_tau(tau), m_delta(delta), m_particleDensity(6*eta/Pi),
	m_size(eqpoints+geompoints+1), m_qgrid(m_size), m_Sq(m_size),
	m_cq(m_size), m_CVertex(m_size), m_SVertex(m_size), m_BorderIndices(m_size),
	m_file()
{
	std::cout << "# Generating iterative wavevector grid with:\n";
	std::cout << "# M = " << eqpoints << " equidistant points\n";
	std::cout << "# G = " << geompoints << " geometric points...\n";
	std::cout << "# Total number of grid points S = " << m_size << "\n";

	generateGrid(max, eqpoints, geompoints);
	sizeUpdate();

	for (int q=0; q<m_size; ++q)
	{
		m_cq[q] = calculateDirectCorrelation(m_qgrid[q]);
		m_Sq[q] = calculateStaticStructure(m_cq[q]);
	}

	determineIntegralBounds();
	initializeVertex();
}

// Only responsible for closing the file again.
Mesh::~Mesh()
{
	m_file.close();
	std::cout << "# Closing Mesh data file...  \n";
}

double Mesh::qback() const
{
	return m_qgrid.back();
}

// We actually want the smallest non-zero wavevector to be returned as q[0]=0 (always)!
double Mesh::qfront() const
{
	return m_qgrid[1];
}

double Mesh::getq(const int index) const
{
	return m_qgrid[index];
}

double Mesh::getSq(const int index) const
{
	return m_Sq[index];
}

double Mesh::getcq(const int index) const
{
	return m_cq[index];
}

double Mesh::getPackingFraction() const
{
	return m_eta;
}

double Mesh::getTau() const
{
	return m_tau;
}

double Mesh::getDelta() const
{
	return m_delta;
}

double Mesh::getParticleDensity() const
{
	return m_particleDensity;
}

int Mesh::getSize() const
{
	return m_size;
}

const array3D& Mesh::getCollectiveVertex()
{
	return m_CVertex;
}

const array3D& Mesh::getTaggedVertex()
{
	return m_SVertex;
}

const array3D& Mesh::getVertex(bool returnCollective)
{
	if(returnCollective)
	{
		return m_CVertex;
	}
	else
	{
		return m_SVertex;
	}
}

const array2D& Mesh::getBorderIndices()
{
	return m_BorderIndices;
}

const std::vector<double>& Mesh::getQgrid()
{
	return m_qgrid;
}

const std::vector<double>& Mesh::getcq()
{
	return m_cq;
}

const std::vector<double>& Mesh::getSq()
{
	return m_Sq;
}

// This function is used to determine the indices of the integration bounds |q-k| and q+k.
// These integration bounds are later needed to retrieve the index range for any given
// combination of q,k when calculating the memory kernel.
void Mesh::determineIntegralBounds()
{
	// In order to guarantee that all of those wavevectors will be on the grid (especially for
	// any non-equidistant grid) all wavevectors are transformed to scaled indices. The spacing
	// of those indices is with respect to the smallest discretization step.
	std::vector<long long> scaled_indices(m_size);	
	for (int q=0; q<m_size; ++q)
	{
		m_BorderIndices[q].resize(m_size);
		scaled_indices[q] = m_qgrid[q]/(m_qgrid[1]-m_qgrid[0])+.5;
	}	
	
	// Loop through q and k for all possible combinations of |q-k| and q+k.
	for (int q=0; q<m_size; ++q)
	{		
		for (int k=0; k<m_size; ++k)
		{
			// Determine |q-k| and q+k in terms of scaled indices. If q+k
			// is greater than the highest grid value, set it to the highest
			// grid value.
			long long lower = fabs(scaled_indices[q]-scaled_indices[k]);
			long long upper = (scaled_indices[q]+scaled_indices[k] < scaled_indices.back()) ? scaled_indices[q]+scaled_indices[k] : scaled_indices.back();

			// Find nearest true scaled index on grid using std::lower_bound functions. The index found here will
			// always be either the actual index (perfect match) or the next scaled index on the grid.
			lower = std::lower_bound(scaled_indices.begin(), scaled_indices.end(), lower) - scaled_indices.begin();
			upper = std::lower_bound(scaled_indices.begin(), scaled_indices.end(), upper) - scaled_indices.begin();

			// For the lower bounds, instead of the next scaled index we take
			// the preceding index
			if (abs(scaled_indices[q]-scaled_indices[k])!=scaled_indices[lower])
			{
				--lower;
			}

			// Store both upper and lower indices for the integration
			// simultaneously. The lower strict triangular half contains
			// the lower integration bound indices, the upper triangular half
			// (including the diagonal) the upper integration bound indices.
			// Note that the lower integration bound indices on the diagonal
			// are always 0!
			if (k==q)
			{
				m_BorderIndices[q][q] = upper;
			}
			else if (k>q)
			{
				m_BorderIndices[q][k] = upper;
				m_BorderIndices[k][q] = lower;
			}
			else
			{
				m_BorderIndices[q][k] = lower;
				m_BorderIndices[k][q] = upper;
			}
		}
	}
}

// This function initializes (and terminally calculates) the vertices for both types
// of correlators. Those vertex arrays are 3D matrices which store the value for
// V(q,k,p) for any given q,k and p=|q-k|...q+k.
void Mesh::initializeVertex()
{
	for (int q=0; q<m_size; ++q)
	{
		// Resize the second dimension of the vertex matrix to be of the same
		// size as the wavevector grid
		m_CVertex[q].resize(m_size);
		m_SVertex[q].resize(m_size);

		for (int k=0; k<m_size; ++k)
		{
			// Using the previously determined integration bound indices,
			// retrieve the integration range indices.
			int lower = 0, upper = 0;
			if (k==q)
			{
				upper = m_BorderIndices[k][k];
			}
			else if (k>q)
			{
				upper = m_BorderIndices[q][k];
				lower = m_BorderIndices[k][q];
			}
			else
			{
				upper = m_BorderIndices[k][q];
				lower = m_BorderIndices[q][k];
			}

			// Implementation of q->0 limit for the vertices.
			if (q==0)
			{
				m_CVertex[q][k].resize(1);
				m_SVertex[q][k].resize(1);

				double dc = 0;
				if (k>0)
				{
					dc = (m_cq[k]-m_cq[k-1])/(m_qgrid[k]-m_qgrid[k-1]);
				}

				m_CVertex[q][k][0] = pow(m_cq[k],2)+2*m_qgrid[k]*m_cq[k]*dc/3;
				m_CVertex[q][k][0] += pow(m_qgrid[k]*dc,2)/5;
				m_CVertex[q][k][0] *= m_Sq[0]*pow(m_Sq[k]*m_qgrid[k],2);
				m_CVertex[q][k][0] *= m_particleDensity/(4*pow(Pi,2));

				m_SVertex[q][k][0] = pow(m_cq[k]*m_qgrid[k]*m_qgrid[k],2);
				m_SVertex[q][k][0] *= m_Sq[k]*m_particleDensity/(2*pow(Pi,2));
			}
			
			// For all other q,k-combinations, calculate the vertex according
			// to the formula. For the first and last entry, use the function
			// which uses the actual wavenumber for p since it might not exist
			// on the grid. For the rest, you may use the function which passes
			// q,k,p as index.
			// TODO: Reference formula in thesis.
			else
			{
				// Resize the last dimension to be of the size of
				// the integration range.
				m_CVertex[q][k].resize(upper-lower+1);
				m_SVertex[q][k].resize(upper-lower+1);

				double plower = fabs(m_qgrid[q]-m_qgrid[k]);
				double pupper = (m_qgrid[q]+m_qgrid[k] < m_qgrid[m_size-1]) ? m_qgrid[q]+m_qgrid[k] : m_qgrid.back();

				m_CVertex[q][k][0] = calcCollectiveVertex(q, k, plower);
				m_CVertex[q][k].back() = calcCollectiveVertex(q, k, pupper);

				m_SVertex[q][k][0] = calcTaggedVertex(q, k, plower);
				m_SVertex[q][k].back() = calcTaggedVertex(q, k, pupper);

				for (int p=lower+1; p<upper; ++p)
				{
					m_CVertex[q][k][p-lower] = calcCollectiveVertex(q, k, p);
					m_SVertex[q][k][p-lower] = calcTaggedVertex(q, k, p);
				}
			}
		}
	}

	std::cout << "# Vertex initialization complete...\n";
}

// Calculate the collective vertex for a given wavevector combination q,k,p which all lie
// on the grid. Hence, these are passed as index. For the calculation the formula is used.
// TODO: Reference formula in thesis.
double Mesh::calcCollectiveVertex(int q, int k, int p)
{
	double a = (m_qgrid[q]+m_qgrid[p])*(m_qgrid[q]-m_qgrid[p])+pow(m_qgrid[k],2); // q^2-p^2+k^2;
	double b = (m_qgrid[q]+m_qgrid[k])*(m_qgrid[q]-m_qgrid[k])+pow(m_qgrid[p],2); // q^2+p^2-k^2;

	double V = pow(a*m_cq[k]+b*m_cq[p],2);
	V *= m_Sq[q]*m_Sq[k]*m_Sq[p]*m_particleDensity;
	V *= m_qgrid[p]*m_qgrid[k];

	return V/(32*Pi*Pi*pow(m_qgrid[q],5));
}

// This is the same function as above, but with p passed as the actual wavenumber. Thus,
// the direct correlation function and static structure factor at wavenumber p have
// to be evaluated directly. The calculation procedure stays the same, though.
double Mesh::calcCollectiveVertex(int q, int k, double p)
{
	double a = (m_qgrid[q]+p)*(m_qgrid[q]-p)+pow(m_qgrid[k],2); // q^2-p^2+k^2;
	double b = (m_qgrid[q]+m_qgrid[k])*(m_qgrid[q]-m_qgrid[k])+pow(p,2); // q^2+p^2-k^2;

	double CP = calculateDirectCorrelation(p);
	double SP = calculateStaticStructure(CP);

	double V = pow(a*m_cq[k]+b*CP,2);
	V *= m_Sq[q]*m_Sq[k]*SP*m_particleDensity;
	V *= p*m_qgrid[k];

	return V/(32*Pi*Pi*pow(m_qgrid[q],5));
}

// Calculate the incoherent vertex for a given wavevector combination q,k,p which all lie
// on the grid. Hence, these are passed as index. For the Calculation the formula is used.
// TODO: Reference formula in thesis.
double Mesh::calcTaggedVertex(int q, int k, int p)
{
	double a = (m_qgrid[q]+m_qgrid[p])*(m_qgrid[q]-m_qgrid[p])+pow(m_qgrid[k],2); // q^2-p^2+k^2;

	double V = pow(a*m_cq[k],2);
	V *= m_Sq[k]*m_particleDensity;
	V *= m_qgrid[p]*m_qgrid[k];

	return V/(16*Pi*Pi*pow(m_qgrid[q],5));
}

// This is the same function as above, but with p passed as the actual wavenumber. Thus,
// the direct correlation function and static structure factor at wavenumber p have
// to be evaluated directly. The calculation procedure stays the same, though.
double Mesh::calcTaggedVertex(int q, int k, double p)
{
	double a = (m_qgrid[q]+p)*(m_qgrid[q]-p)+pow(m_qgrid[k],2); // q^2-p^2+k^2;

	double V = pow(a*m_cq[k],2);
	V *= m_Sq[k]*m_particleDensity;
	V *= p*m_qgrid[k];

	return V/(16*Pi*Pi*pow(m_qgrid[q],5));
}

// Implementation of the output stream operator overloading. The output consists of all
// values of the direct correlation function and static structure factor at a given q value,
// separated by a tab.
std::ostream& operator << (std::ostream &out, const Mesh &array)
{
	for (int q=0; q<array.getSize(); ++q)
	{
		out << array.getq(q) << "\t" << array.getcq(q);
		out << "\t" << array.getSq(q) << "\n";
	}	

	return out; 
}

void Mesh::openFile(std::string filename)
{
	m_file.open(filename);
	std::cout << "# Opening Mesh data file " << filename << "...\n";
}

// This function writes the data to the file using the previously defined << operator.
void Mesh::writeDatatoFile()
{
	std::cout << "# Writing Mesh data to file...\n";
	m_file << "# eta = " << m_eta;
	m_file << "\n# tau = " << m_tau;
	m_file << "\n# delta = " << m_delta;
	m_file << "\n# rho = " << m_particleDensity;
	m_file << "\n# q\tc[q]\tS[q]\n";

	m_file << *this;
}

void Mesh::writeDatatoFile(std::vector<double>& omegaC, std::vector<double>& omegaT)
{
	std::cout << "# Writing Mesh data to file...\n";
	m_file << "# eta = " << m_eta;
	m_file << "\n# tau = " << m_tau;
	m_file << "\n# delta = " << m_delta;
	m_file << "\n# rho = " << m_particleDensity;
	m_file << "\n# q\tc[q]\tS[q]\n";
	
	for (int q=0; q<m_size; ++q)
	{
		m_file << m_qgrid[q] << "\t" << m_cq[q] << "\t";
		m_file << m_Sq[q] << "\t";
		m_file << omegaC[q] << "\t" << omegaT[q] <<"\n";
	}
}

// Calculation of the direct correlation function at any given wavevector (actual wavenumber,
// not the index!) using the formula generated by a Mathematica notebook (by M. Sperl).
double Mesh::calculateDirectCorrelation(double wavevector)
{
	double q = !(wavevector<0.02) ? wavevector : 0.02;
	double a0, b0, c0, a1, b1, c1, a, b, Sqinv;
	a0 = (1+2*m_eta)/pow(1-m_eta,2)-12*m_tau*m_eta/(1-m_eta);
	b0 = -3*m_eta/(2*pow(1-m_eta,2))+6*m_tau*m_eta/(1-m_eta);;
	c0 = m_tau-a0/2-b0;
	a1 = m_tau*6*m_eta/pow(1-m_eta,2)*(-2+m_eta*(5-12*c0)+12*pow(m_eta,2)*c0);
	b1 = m_tau*9*m_eta/pow(1-m_eta,2)*(1-2*m_eta+4*m_eta*c0-4*pow(m_eta,2)*c0);
	c1 = -a1/2-b1+m_tau/2+6*m_tau*m_eta*c0;
	a = a0+m_delta*a1;
	b = b0+m_delta*b1;

	Sqinv = 1 + (144*Power(a,2)*Power(m_eta,2))/Power(q,6)
		+ (144*Power(a,2)*Power(m_eta,2))/Power(q,4)
		+ (288*a*b*Power(m_eta,2))/Power(q,4)
		+ (144*Power(b,2)*Power(m_eta,2))/Power(q,4)+(24*b*m_eta)/Power(q,2)
		+ (36*Power(a,2)*Power(m_eta,2))/Power(q,2)
		+ (144*a*b*Power(m_eta,2))/Power(q,2)
		+ (144*Power(b,2)*Power(m_eta,2))/Power(q,2)
		- (288*a*Power(m_eta,2)*m_tau)/Power(q,4)
		- (144*a*m_delta*Power(m_eta,2)*m_tau)/Power(q,4)
		- (1728*a*c0*m_delta*Power(m_eta,3)*m_tau)/Power(q,4)
		- (144*a*Power(m_eta,2)*m_tau)/Power(q,2)
		- (288*b*Power(m_eta,2)*m_tau)/Power(q,2)
		- (72*a*m_delta*Power(m_eta,2)*m_tau)/Power(q,2)
		- (144*b*m_delta*Power(m_eta,2)*m_tau)/Power(q,2)
		- (864*a*c0*m_delta*Power(m_eta,3)*m_tau)/Power(q,2)
		- (1728*b*c0*m_delta*Power(m_eta,3)*m_tau)/Power(q,2)
		+ (3456*b*Power(m_eta,3)*Power(m_tau,2))/(Power(m_delta,2)*Power(q,6))
		+ (3456*a*Power(m_eta,3)*Power(m_tau,2))/(m_delta*Power(q,6))
		+ (288*Power(m_eta,2)*Power(m_tau,2))/(Power(m_delta,2)*Power(q,4))
		- (1728*b*Power(m_eta,3)*Power(m_tau,2))/Power(q,4)
		+ (1728*a*Power(m_eta,3)*Power(m_tau,2))/(m_delta*Power(q,4))
		+ (3456*b*Power(m_eta,3)*Power(m_tau,2))/(m_delta*Power(q,4))
		- (576*a*m_delta*Power(m_eta,3)*Power(m_tau,2))/Power(q,4)
		+ (144*m_delta*Power(m_eta,2)*Power(m_tau,2))/Power(q,2)
		- (288*a*m_delta*Power(m_eta,3)*Power(m_tau,2))/Power(q,2)
		- (576*b*m_delta*Power(m_eta,3)*Power(m_tau,2))/Power(q,2)
		+ (1728*c0*m_delta*Power(m_eta,3)*Power(m_tau,2))/Power(q,2)
		- (1728*Power(m_eta,3)*Power(m_tau,3))/Power(q,4)
		- (3456*Power(m_eta,3)*Power(m_tau,3))/(m_delta*Power(q,4))
		- (20736*c0*Power(m_eta,4)*Power(m_tau,3))/Power(q,4)
		+ (576*m_delta*Power(m_eta,3)*Power(m_tau,3))/Power(q,2)
		+ (20736*Power(m_eta,4)*Power(m_tau,4))/(Power(m_delta,4)*Power(q,8))
		- (1728*Power(m_eta,4)*Power(m_tau,4))/Power(q,4)
		- (288*Power(a,2)*Power(m_eta,2)*Cos(q))/Power(q,6)
		- (144*Power(a,2)*Power(m_eta,2)*Cos(q))/Power(q,4)
		- (576*a*b*Power(m_eta,2)*Cos(q))/Power(q,4)
		- (288*Power(b,2)*Power(m_eta,2)*Cos(q))/Power(q,4)
		- (24*a*m_eta*Cos(q))/Power(q,2)
		- (24*b*m_eta*Cos(q))/Power(q,2)
		- (288*a*Power(m_eta,2)*m_tau*Cos(q))/(m_delta*Power(q,6))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Cos(q))/(m_delta*Power(q,6))
		+ (288*a*Power(m_eta,2)*m_tau*Cos(q))/Power(q,4)
		- (144*a*Power(m_eta,2)*m_tau*Cos(q))/(m_delta*Power(q,4))
		- (576*b*Power(m_eta,2)*m_tau*Cos(q))/(m_delta*Power(q,4))
		+ (144*a*m_delta*Power(m_eta,2)*m_tau*Cos(q))/Power(q,4)
		- (3456*b*c0*Power(m_eta,3)*m_tau*Cos(q))/Power(q,4)
		+ (1728*a*c0*Power(m_eta,3)*m_tau*Cos(q))/(m_delta*Power(q,4))
		+ (3456*b*c0*Power(m_eta,3)*m_tau*Cos(q))/(m_delta*Power(q,4))
		+ (1728*a*c0*m_delta*Power(m_eta,3)*m_tau*Cos(q))/Power(q,4)
		- (24*m_eta*m_tau*Cos(q))/(m_delta*Power(q,2))
		- (288*c0*Power(m_eta,2)*m_tau*Cos(q))/Power(q,2)
		- (3456*a*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/(Power(m_delta,2)*Power(q,6))
		- (3456*b*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/(Power(m_delta,2)*Power(q,6))
		- (3456*a*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/(m_delta*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Cos(q))/Power(q,4)
		+ (288*Power(m_eta,2)*Power(m_tau,2)*Cos(q))/(m_delta*Power(q,4))
		+ (1728*a*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/Power(q,4)
		+ (1728*b*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/Power(q,4)
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/(m_delta*Power(q,4))
		+ (576*a*m_delta*Power(m_eta,3)*Power(m_tau,2)*Cos(q))/Power(q,4)
		- (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Cos(q))/Power(q,4)
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(q))/(Power(m_delta,3)*Power(q,6))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(q))/(Power(m_delta,2)*Power(q,6))
		+ (576*Power(m_eta,3)*Power(m_tau,3)*Cos(q))/Power(q,4)
		+ (1728*Power(m_eta,3)*Power(m_tau,3)*Cos(q))/(m_delta*Power(q,4))
		+ (13824*c0*Power(m_eta,4)*Power(m_tau,3)*Cos(q))/Power(q,4)
		+ (144*Power(a,2)*Power(m_eta,2)*Power(Cos(q),2))/Power(q,6)
		+ (144*Power(a,2)*Power(m_eta,2)*Power(Cos(q),2))/Power(q,4)
		+ (288*a*b*Power(m_eta,2)*Power(Cos(q),2))/Power(q,4)
		+ (144*Power(b,2)*Power(m_eta,2)*Power(Cos(q),2))/Power(q,4)
		+ (288*a*Power(m_eta,2)*m_tau*Power(Cos(q),2))/(m_delta*Power(q,6))
		- (3456*a*c0*Power(m_eta,3)*m_tau*Power(Cos(q),2))/(m_delta*Power(q,6))
		+ (288*a*Power(m_eta,2)*m_tau*Power(Cos(q),2))/(m_delta*Power(q,4))
		+ (288*b*Power(m_eta,2)*m_tau*Power(Cos(q),2))/(m_delta*Power(q,4))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Power(Cos(q),2))/Power(q,4)
		+ (3456*b*c0*Power(m_eta,3)*m_tau*Power(Cos(q),2))/Power(q,4)
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Cos(q),2))/(Power(m_delta,2)*Power(q,6))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Cos(q),2))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Cos(q),2))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Cos(q),2))/(Power(m_delta,2)*Power(q,4))
		+ (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Cos(q),2))/(m_delta*Power(q,4))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Cos(q),2))/Power(q,4)
		- (3456*b*Power(m_eta,3)*Power(m_tau,2)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,4))
		- (41472*Power(m_eta,4)*Power(m_tau,4)*Cos(m_delta*q))/(Power(m_delta,4)*Power(q,8))
		+ (20736*Power(m_eta,4)*Power(m_tau,4)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (3456*a*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (3456*b*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(q)*Cos(m_delta*q))/(Power(m_delta,3)*Power(q,6))
		+ (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos(q)*Cos(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(m_eta,4)*Power(m_tau,4)*Power(Cos(m_delta*q),2))/(Power(m_delta,4)*Power(q,8))
		+ (288*a*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,6))
		- (3456*a*c0*Power(m_eta,3)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,6))
		+ (288*b*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q))/Power(q,4)
		+ (144*a*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		+ (576*b*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (1728*a*c0*Power(m_eta,3)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (3456*b*c0*Power(m_eta,3)*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		+ (24*m_eta*m_tau*Cos((1 + m_delta)*q))/Power(q,2)
		+ (24*m_eta*m_tau*Cos((1 + m_delta)*q))/(m_delta*Power(q,2))
		- (144*Power(m_eta,2)*Power(m_tau,2)*Cos((1 + m_delta)*q))/Power(q,4)
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		+ (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Cos((1 + m_delta)*q))/Power(q,4)
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Cos((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,6))
		+ (6912*Power(m_eta,3)*Power(m_tau,3)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (2304*Power(m_eta,3)*Power(m_tau,3)*Cos((1 + m_delta)*q))/Power(q,4)
		- (1728*Power(m_eta,3)*Power(m_tau,3)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		+ (6912*c0*Power(m_eta,4)*Power(m_tau,3)*Cos((1 + m_delta)*q))/Power(q,4)
		- (288*a*Power(m_eta,2)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,6))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,6))
		- (288*a*Power(m_eta,2)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/Power(q,4)
		- (288*b*Power(m_eta,2)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/Power(q,4)
		- (288*a*Power(m_eta,2)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (288*b*Power(m_eta,2)*m_tau*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		+ (6912*c0*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (41472*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,4))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/Power(q,4)
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Cos((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(m_delta*q)*Cos((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,6))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(m_delta*q)*Cos((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/Power(q,4)
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,4))
		+ (288*Power(m_eta,2)*Power(m_tau,2)*Power(Cos((1 + m_delta)*q),2))/(m_delta*Power(q,4))
		- (288*Power(a,2)*Power(m_eta,2)*Sin(q))/Power(q,5)
		+ (24*a*m_eta*Sin(q))/Power(q,3)
		- (144*Power(a,2)*Power(m_eta,2)*Sin(q))/Power(q,3)
		- (432*a*b*Power(m_eta,2)*Sin(q))/Power(q,3)
		- (288*Power(b,2)*Power(m_eta,2)*Sin(q))/Power(q,3)
		- (288*a*Power(m_eta,2)*m_tau*Sin(q))/(m_delta*Power(q,5))
		+ (288*b*Power(m_eta,2)*m_tau*Sin(q))/(m_delta*Power(q,5))
		- (3456*a*c0*Power(m_eta,3)*m_tau*Sin(q))/Power(q,5)
		- (3456*b*c0*Power(m_eta,3)*m_tau*Sin(q))/(m_delta*Power(q,5))
		+ (24*m_eta*m_tau*Sin(q))/(m_delta*Power(q,3))
		+ (288*a*Power(m_eta,2)*m_tau*Sin(q))/Power(q,3)
		+ (288*b*Power(m_eta,2)*m_tau*Sin(q))/Power(q,3)
		- (144*a*Power(m_eta,2)*m_tau*Sin(q))/(m_delta*Power(q,3))
		- (288*b*Power(m_eta,2)*m_tau*Sin(q))/(m_delta*Power(q,3))
		- (288*c0*Power(m_eta,2)*m_tau*Sin(q))/(m_delta*Power(q,3))
		+ (144*a*m_delta*Power(m_eta,2)*m_tau*Sin(q))/Power(q,3)
		+ (144*b*m_delta*Power(m_eta,2)*m_tau*Sin(q))/Power(q,3)
		- (1728*a*c0*Power(m_eta,3)*m_tau*Sin(q))/Power(q,3)
		- (3456*b*c0*Power(m_eta,3)*m_tau*Sin(q))/Power(q,3)
		+ (1728*a*c0*m_delta*Power(m_eta,3)*m_tau*Sin(q))/Power(q,3)
		+ (1728*b*c0*m_delta*Power(m_eta,3)*m_tau*Sin(q))/Power(q,3)
		+ (3456*a*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/(Power(m_delta,2)*Power(q,7))
		- (1728*a*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/Power(q,5)
		- (3456*a*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/(m_delta*Power(q,5))
		- (3456*b*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/(m_delta*Power(q,5))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (288*Power(m_eta,2)*Power(m_tau,2)*Sin(q))/(m_delta*Power(q,3))
		+ (5184*c0*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (576*a*m_delta*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (576*b*m_delta*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (1728*c0*m_delta*Power(m_eta,3)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (20736*Power(c0,2)*m_delta*Power(m_eta,4)*Power(m_tau,2)*Sin(q))/Power(q,3)
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(q))/(Power(m_delta,3)*Power(q,7))
		- (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Sin(q))/(Power(m_delta,3)*Power(q,7))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(q))/(Power(m_delta,2)*Power(q,5))
		- (1728*Power(m_eta,3)*Power(m_tau,3)*Sin(q))/(m_delta*Power(q,5))
		- (20736*c0*Power(m_eta,4)*Power(m_tau,3)*Sin(q))/(m_delta*Power(q,5))
		+ (576*Power(m_eta,3)*Power(m_tau,3)*Sin(q))/Power(q,3)
		+ (6912*c0*m_delta*Power(m_eta,4)*Power(m_tau,3)*Sin(q))/Power(q,3)
		- (3456*a*Power(m_eta,3)*Power(m_tau,2)*Cos(m_delta*q)*Sin(q))/(Power(m_delta,2)*Power(q,7))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(m_delta*q)*Sin(q))/(Power(m_delta,3)*Power(q,7))
		+ (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos(m_delta*q)*Sin(q))/(Power(m_delta,3)*Power(q,7))
		+ (288*a*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q)*Sin(q))/Power(q,5)
		- (288*b*Power(m_eta,2)*m_tau*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		+ (3456*b*c0*Power(m_eta,3)*m_tau*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		+ (288*Power(m_eta,2)*Power(m_tau,2)*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		- (6912*c0*Power(m_eta,3)*Power(m_tau,2)*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		+ (41472*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Cos((1 + m_delta)*q)*Sin(q))/(m_delta*Power(q,5))
		+ (144*Power(a,2)*Power(m_eta,2)*Power(Sin(q),2))/Power(q,6)
		+ (144*Power(a,2)*Power(m_eta,2)*Power(Sin(q),2))/Power(q,4)
		+ (288*a*b*Power(m_eta,2)*Power(Sin(q),2))/Power(q,4)
		+ (144*Power(b,2)*Power(m_eta,2)*Power(Sin(q),2))/Power(q,4)
		+ (288*a*Power(m_eta,2)*m_tau*Power(Sin(q),2))/(m_delta*Power(q,6))
		- (3456*a*c0*Power(m_eta,3)*m_tau*Power(Sin(q),2))/(m_delta*Power(q,6))
		+ (288*a*Power(m_eta,2)*m_tau*Power(Sin(q),2))/(m_delta*Power(q,4))
		+ (288*b*Power(m_eta,2)*m_tau*Power(Sin(q),2))/(m_delta*Power(q,4))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Power(Sin(q),2))/Power(q,4)
		+ (3456*b*c0*Power(m_eta,3)*m_tau*Power(Sin(q),2))/Power(q,4)
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Sin(q),2))/(Power(m_delta,2)*Power(q,6))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Sin(q),2))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Sin(q),2))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Sin(q),2))/(Power(m_delta,2)*Power(q,4))
		+ (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Sin(q),2))/(m_delta*Power(q,4))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Sin(q),2))/Power(q,4)
		- (3456*a*Power(m_eta,3)*Power(m_tau,2)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,7))
		- (1728*a*Power(m_eta,3)*Power(m_tau,2)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,5))
		- (3456*b*Power(m_eta,3)*Power(m_tau,2)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,5))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,5))
		+ (1728*Power(m_eta,3)*Power(m_tau,3)*Sin(m_delta*q))/(m_delta*Power(q,5))
		+ (20736*c0*Power(m_eta,4)*Power(m_tau,3)*Sin(m_delta*q))/(m_delta*Power(q,5))
		- (41472*Power(m_eta,4)*Power(m_tau,4)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,7))
		+ (6912*Power(m_eta,4)*Power(m_tau,4)*Sin(m_delta*q))/(m_delta*Power(q,5))
		+ (3456*a*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,7))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(q)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,7))
		- (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos(q)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,7))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Cos((1 + m_delta)*q)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,7))
		+ (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos((1 + m_delta)*q)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,7))
		+ (3456*a*Power(m_eta,3)*Power(m_tau,2)*Sin(q)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (3456*b*Power(m_eta,3)*Power(m_tau,2)*Sin(q)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(q)*Sin(m_delta*q))/(Power(m_delta,3)*Power(q,6))
		+ (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Sin(q)*Sin(m_delta*q))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(m_eta,4)*Power(m_tau,4)*Power(Sin(m_delta*q),2))/(Power(m_delta,4)*Power(q,8))
		+ (288*a*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/Power(q,5)
		+ (288*a*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (288*b*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		+ (3456*b*c0*Power(m_eta,3)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (24*m_eta*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,3))
		+ (144*a*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/Power(q,3)
		+ (288*b*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/Power(q,3)
		+ (144*a*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,3))
		+ (288*b*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,3))
		+ (288*c0*Power(m_eta,2)*m_tau*Sin((1 + m_delta)*q))/(m_delta*Power(q,3))
		- (432*Power(m_eta,2)*Power(m_tau,2)*Sin((1 + m_delta)*q))/Power(q,3)
		- (288*Power(m_eta,2)*Power(m_tau,2)*Sin((1 + m_delta)*q))/(m_delta*Power(q,3))
		- (144*m_delta*Power(m_eta,2)*Power(m_tau,2)*Sin((1 + m_delta)*q))/Power(q,3)
		- (1728*c0*Power(m_eta,3)*Power(m_tau,2)*Sin((1 + m_delta)*q))/Power(q,3)
		- (1728*c0*m_delta*Power(m_eta,3)*Power(m_tau,2)*Sin((1 + m_delta)*q))/Power(q,3)
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Sin((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,7))
		+ (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Sin((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,7))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,5))
		+ (5184*Power(m_eta,3)*Power(m_tau,3)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (20736*c0*Power(m_eta,4)*Power(m_tau,3)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (576*Power(m_eta,3)*Power(m_tau,3)*Sin((1 + m_delta)*q))/Power(q,3)
		- (576*m_delta*Power(m_eta,3)*Power(m_tau,3)*Sin((1 + m_delta)*q))/Power(q,3)
		- (288*a*Power(m_eta,2)*m_tau*Cos(q)*Sin((1 + m_delta)*q))/Power(q,5)
		+ (288*b*Power(m_eta,2)*m_tau*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (3456*a*c0*Power(m_eta,3)*m_tau*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (3456*b*c0*Power(m_eta,3)*m_tau*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		+ (6912*c0*Power(m_eta,3)*Power(m_tau,2)*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		- (41472*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Cos(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,5))
		+ (3456*Power(m_eta,3)*Power(m_tau,3)*Cos(m_delta*q)*Sin((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,7))
		- (41472*c0*Power(m_eta,4)*Power(m_tau,3)*Cos(m_delta*q)*Sin((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,7))
		- (288*a*Power(m_eta,2)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,6))
		+ (3456*a*c0*Power(m_eta,3)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,6))
		- (288*a*Power(m_eta,2)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/Power(q,4)
		- (288*b*Power(m_eta,2)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/Power(q,4)
		- (288*a*Power(m_eta,2)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (288*b*Power(m_eta,2)*m_tau*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		+ (6912*c0*Power(m_eta,3)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (41472*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,4))
		- (288*Power(m_eta,2)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/Power(q,4)
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Sin(q)*Sin((1 + m_delta)*q))/(m_delta*Power(q,4))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(m_delta*q)*Sin((1 + m_delta)*q))/(Power(m_delta,3)*Power(q,6))
		- (3456*Power(m_eta,3)*Power(m_tau,3)*Sin(m_delta*q)*Sin((1 + m_delta)*q))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		- (3456*c0*Power(m_eta,3)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		+ (20736*Power(c0,2)*Power(m_eta,4)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,6))
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/Power(q,4)
		+ (144*Power(m_eta,2)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/(Power(m_delta,2)*Power(q,4))
		+ (288*Power(m_eta,2)*Power(m_tau,2)*Power(Sin((1 + m_delta)*q),2))/(m_delta*Power(q,4));

	return (1-Sqinv)/m_particleDensity;
}

// Calculation of the static structure factor by passing the value of the direct correlation
// function at a given wavenumber using the Ornstein-Zernicke relation.
double Mesh::calculateStaticStructure(double directcorrfct)
{
	return 1/(1-directcorrfct*m_particleDensity);
}

// Update all grid sizes to be equal to the wavevector grid size. This needs to be done
// after the grid generation was called.
void Mesh::sizeUpdate()
{
	m_size = m_qgrid.size();
	m_Sq.resize(m_size);
	m_cq.resize(m_size);
	m_BorderIndices.resize(m_size);
	m_CVertex.resize(m_size);
	m_SVertex.resize(m_size);
}

// This function works very similar to the one above, but instead there will only be one
// large equidistant grid generated. Then, the smallest wavenumber will be halved and
// inserted in the front etc...
// This produces grids with smaller sizes (m_size = size + iterations + 1).
void Mesh::generateGrid(double max, int size, int iterations)
{
	double spacing = max/size;
	m_qgrid.resize(m_size);

	m_qgrid[0] = 0;

	for (int s=1; s<iterations+1; ++s)
	{
		m_qgrid[s] = spacing/pow(2,iterations+1-s);
	}

	for (int q=iterations+1; q<m_size; ++q)
	{
		m_qgrid[q] = spacing;
		spacing += max/size;
	}
}
