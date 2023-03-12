// Structure.cpp : definition of the class member functions used in "Structure.h"; main
//		   functionality for structure/wavevector related quantities (jn)
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
// 2016-10-17 : moved c_swf_first & s_swf_first (by ms) functions into the Structure
//		class to calculate cq/Sq; slightly changed function implementation (jn)
// 2016-11-16 : changed vertices matrices from copied form to the form presented in the
// 		diploma thesis of M. Sperl (jn)
// 2016-11-30 : current implementation works for equidistant wavevector grids; vertex
// 		integration/calculation still needs to be figured out for non-equidistant
// 		grids in terms of integration spacings (jn)
// 2017-01-10 : implemented the logarithmic grid for the wavevectors due to problems
// 		with the equidistant grid iterationn method (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
// - figure out vertex integration/adapt vertex calculation for non-equidistant grids
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
#include "constants.h"
#include "Structure.h"

Structure::Structure(double max, int size, int iterations, double eta, double tau, double delta) :
	m_eta(eta), m_tau(tau), m_delta(delta),
	m_size(size), m_qgrid(m_size), m_Sq(m_size),
	m_cq(m_size), m_Vertex(m_size), m_BorderIndices(m_size),
	m_filename("Structure.dat"), 
	m_file(m_filename.c_str(), std::ios::out | std::ios::trunc)
{
	std::cout << "# Opening file " << m_filename << "...\n";
	std::cout << "# Generating iterative wavevector grid with " << size;
	std::cout << " wavevectors and " << iterations << " iterations...\n";

	double spacing = max/size;

	for (int q=0; q<size; ++q)
	{
		m_qgrid[q] = spacing;
		spacing += max/size;
	}

	for (int s=0; s<iterations; ++s)
	{
		max /= 2;
		m_qgrid.erase(m_qgrid.begin(), m_qgrid.begin()+size/2);
		std::vector<double> tmp(size);
		spacing = max/size;
		for (int q=0; q<size; ++q)
		{
			tmp[q] = spacing;
			spacing += max/size;
		}
		
		m_qgrid.insert(m_qgrid.begin(), tmp.begin(), tmp.end());
	}

	m_qgrid.insert(m_qgrid.begin(), 0);

	m_size = m_qgrid.size();
	m_Sq.resize(m_size);
	m_cq.resize(m_size);
	m_BorderIndices.resize(m_size);
	m_Vertex.resize(m_size);

	for (int q=0; q<m_size; ++q)
	{
		m_cq[q] = calculateDirectCorrelation(m_qgrid[q]);
		m_Sq[q] = calculateStaticStructure(m_cq[q]);
	}

	writeDatatoFile();
	determineIntegralBounds();
	initializeVertex();
}

Structure::~Structure()
{
	m_file.close(); std::cout << "# Closing " << m_filename << "...  \n";
}

double Structure::qback() const
{
	return m_qgrid.back();
}

double Structure::qfront() const
{
	return m_qgrid[1];
}

double Structure::getq(const int index) const
{
	return m_qgrid[index];
}

double Structure::getSq(const int index) const
{
	return m_Sq[index];
}

double Structure::getcq(const int index) const
{
	return m_cq[index];
}

double Structure::getPackingFraction() const
{
	return m_eta;
}

double Structure::getTau() const
{
	return m_tau;
}

double Structure::getDelta() const
{
	return m_delta;
}

int Structure::getSize() const
{
	return m_size;
}

const array3D& Structure::getVertex()
{
	return m_Vertex;
}

const array2D& Structure::getBorderIndices()
{
	return m_BorderIndices;
}

const std::vector<double>& Structure::getQgrid()
{
	return m_qgrid;
}

const std::vector<double>& Structure::getcq()
{
	return m_cq;
}

const std::vector<double>& Structure::getSq()
{
	return m_Sq;
}

void Structure::determineIntegralBounds()
{
	std::vector<int> scaled_indices(m_size);	
	for (int q=0; q<m_size; ++q)
	{
		m_BorderIndices[q].resize(m_size);
		scaled_indices[q] = m_qgrid[q]/(m_qgrid[1]-m_qgrid[0])+.5;
	}	

	for (int q=0; q<m_size; ++q)
	{
		for (int k=0; k<m_size; ++k)
		{
			int lower = fabs(scaled_indices[q]-scaled_indices[k]);
			int upper = (scaled_indices[q]+scaled_indices[k] < scaled_indices.back()) ? scaled_indices[q]+scaled_indices[k] : scaled_indices.back();

			lower = std::lower_bound(scaled_indices.begin(), scaled_indices.end(), lower) - scaled_indices.begin();
			upper = std::lower_bound(scaled_indices.begin(), scaled_indices.end(), upper) - scaled_indices.begin();

			if (fabs(scaled_indices[q]-scaled_indices[k])!=scaled_indices[lower])
			{
				--lower;
			}

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

void Structure::initializeVertex()
{
	std::ofstream file("Vertex.dat");
	for (int q=0; q<m_size; ++q)
	{
		m_Vertex[q].resize(m_size);

		for (int k=0; k<m_size; ++k)
		{
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

			m_Vertex[q][k].resize(upper-lower+1);

			if (q==0)
			{
				for (int p=0; p<upper-lower+1; ++p)
				{
					m_Vertex[q][k][p] = 0;
				}
			}
			else
			{
				double plower = fabs(m_qgrid[q]-m_qgrid[k]);
				double pupper = (m_qgrid[q]+m_qgrid[k] < m_qgrid[m_size-1]) ? m_qgrid[q]+m_qgrid[k] : m_qgrid.back();

				m_Vertex[q][k][0] = calcVertex(q, k, plower);
				m_Vertex[q][k].back() = calcVertex(q, k, pupper);

				for (int p=lower+1; p<upper; ++p)
				{
					m_Vertex[q][k][p-lower] = calcVertex(q, k, p);
				}
			}
		}
	}

	file.close();

	std::cout << "# Vertex initialization complete...\n";
}

double Structure::getInterpolationValue(double x, double leftIP, double leftIV, double rightIP, double rightIV)
{
	return leftIV + (x-leftIP)*(rightIV-leftIV)/(rightIP-leftIP);
}

double Structure::calcVertex(int q, int k, int p)
{
	double a = (m_qgrid[q]+m_qgrid[p])*(m_qgrid[q]-m_qgrid[p])+pow(m_qgrid[k],2); // q^2-p^2+k^2;
	double b = (m_qgrid[q]+m_qgrid[k])*(m_qgrid[q]-m_qgrid[k])+pow(m_qgrid[p],2); // q^2+p^2-k^2;
	double particleDensity = 6*m_eta/Pi;

	double V = pow(a*m_cq[k]+b*m_cq[p],2);
	V *= m_Sq[q]*m_Sq[k]*m_Sq[p]*particleDensity;
	V *= m_qgrid[p]*m_qgrid[k];

	return V/(32*Pi*Pi*pow(m_qgrid[q],5));
}

double Structure::calcVertex(int q, int k, double p)
{
	double a = (m_qgrid[q]+p)*(m_qgrid[q]-p)+pow(m_qgrid[k],2); // q^2-p^2+k^2;
	double b = (m_qgrid[q]+m_qgrid[k])*(m_qgrid[q]-m_qgrid[k])+pow(p,2); // q^2+p^2-k^2;
	double particleDensity = 6*m_eta/Pi;

	double CP = calculateDirectCorrelation(p);
	double SP = calculateStaticStructure(CP);

	double V = pow(a*m_cq[k]+b*CP,2);
	V *= m_Sq[q]*m_Sq[k]*SP*particleDensity;
	V *= p*m_qgrid[k];

	return V/(32*Pi*Pi*pow(m_qgrid[q],5));
}

std::ostream& operator << (std::ostream &out, const Structure &array)
{
	out << "# q\tc[q]\tS[q]\n";
	for (int q=0; q<array.getSize(); ++q)
	{
		out << array.getq(q) << "\t" << array.getcq(q) << "\t" << array.getSq(q) << "\n";
	}
	return out; 
}

void Structure::writeDatatoFile()
{
	std::cout << "# Writing data to file " << m_filename << "...\n";
	m_file << *this;
}

double Structure::calculateDirectCorrelation(double wavevector)
{
	double q = !(wavevector<0.01) ? wavevector : 0.01;
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

	return (1-Sqinv)*(Pi/(6*m_eta));
}

double Structure::calculateStaticStructure(double directcorrfct)
{
	return 1/(1-directcorrfct*(6*m_eta)/Pi);
}
