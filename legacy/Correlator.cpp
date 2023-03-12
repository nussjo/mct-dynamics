// Correlator.cpp :
// Implementation of functionalities for Correlator class and subclasses.
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
// 2016-11-18 : the problem with the calculation of the memory kernel is due to scaling
// 		inconsistencies ---> working with indices of wavenumbers like in the M.
// 		Sperl code works but ONLY for equidistant grids (iterative grid?) (jn)
// 2017-04-05 : calculation of memory kernel integral now done with a trapezoidal rule and
// 		is working properly (jn)
//
// --------------------------------------------------------------------------------------
//
// TODO-LIST :
//
// ======================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "mdefs.h"
#include "Correlator.h"
#include "Utils.h"

Decimation::Decimation(int tsteps, int qsteps, double h) :
	m_spacing(h), m_tsteps(tsteps), m_qsteps(qsteps),
	m_time(m_tsteps,0),
	m_values(m_qsteps, std::vector<double>(m_tsteps,1)),
	m_kernel(m_qsteps, std::vector<double>(m_tsteps,1)),
	m_dP(m_qsteps, std::vector<double>(m_tsteps,0)),
	m_dM(m_qsteps, std::vector<double>(m_tsteps,0))
{
	for (int i=1; i<m_tsteps; ++i)	
	{
		m_time[i] = i*m_spacing;
	}
}

void Decimation::decimate()
{
	for (int i=0; i<m_tsteps/2; ++i)
	{
		m_time[i] = m_time[2*i];

		for (int q=0; q<m_qsteps; ++q)
		{
			m_values[q][i] = m_values[q][2*i];
			m_kernel[q][i] = m_kernel[q][2*i];
			m_dP[q][i] = .5*(m_dP[q][2*i]+m_dP[q][2*i-1]);
			m_dM[q][i] = .5*(m_dM[q][2*i]+m_dM[q][2*i-1]);
		}
	}

	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = 0;

		for (int q=0; q<m_qsteps; ++q)	
		{
			m_values[q][i] = 0; m_kernel[q][i] = 0;
			m_dP[q][i] = 0; m_dM[q][i] = 0;
		}
	}

	m_spacing *= 2;
}

int Decimation::getTimesteps() const
{
	return m_tsteps;
}

int Decimation::getQgridsteps() const
{
	return m_qsteps;
}

double Decimation::getSpacing() const
{
	return m_spacing;
}

const array2D& Decimation::getArrayValues() const
{
	return m_values;
}

const array2D& Decimation::getKernelValues() const
{
	return m_kernel;
}

Decimation::~Decimation()
{
}

Correlator::Correlator(int tsteps, Mesh &Q, double h, double vtherm, double nu) :
	Decimation(tsteps, Q.getSize(), h), m_qmesh(Q),
	m_D0(vtherm*vtherm/nu),	m_nu(nu), m_isNewtonian(false),
	m_Omega(m_qsteps,vtherm), m_vertInt(m_qsteps,0),
	m_phifile(), m_kernelfile()
{
}

Correlator::~Correlator()
{
}

void Correlator::openFiles(int decimations, double thermalvelocity, std::string phifilename, std::string kernelfilename)
{
	m_phifile.open(phifilename);
	m_kernelfile.open(kernelfilename);

	std::cout << "# Opening file " << phifilename << "...\n";
	std::cout << "# Opening file " << kernelfilename << "...\n";

	writeFileHeader(decimations, thermalvelocity);
}

void Correlator::writeFileHeader(int decimations, double thermalvelocity)
{
	m_phifile << "# eta = " << m_qmesh.getPackingFraction();
	m_phifile << ", tau = " << m_qmesh.getTau();
	m_phifile << ", delta = " << m_qmesh.getDelta() << "\n";
	m_phifile << "# N = " << m_tsteps << ", D = " << decimations;
	m_phifile << ", S = " << m_qsteps << "\n";
	m_phifile << "# q_min = " << m_qmesh.qfront() << ", q_max = ";
	m_phifile << m_qmesh.qback() << "\n";
	m_phifile << "# h_initial = " << m_spacing << ", v_therm = ";
	m_phifile << thermalvelocity << "\tD0 = " << m_D0 << "\n";

	m_kernelfile << "# eta = " << m_qmesh.getPackingFraction();
	m_kernelfile << ", tau = " << m_qmesh.getTau();
	m_kernelfile << ", delta = " << m_qmesh.getDelta() << "\n";
	m_kernelfile << "# N = " << m_tsteps << ", D = " << decimations;
	m_kernelfile << ", S = " << m_qsteps << "\n";
	m_kernelfile << "# q_min = " << m_qmesh.qfront() << ", q_max = ";
	m_kernelfile << m_qmesh.qback() << "\n";
	m_kernelfile << "# h_initial = " << m_spacing << ", v_therm = ";
	m_kernelfile << thermalvelocity << "\tD0 = " << m_D0 << "\n";
	m_kernelfile << "# Prefactor m(q->0) c = " << calculateLowQPrefactor() << "\n";
	m_kernelfile << "# m(q->0) results stored in memory kernel array m_kernel[0][]\n";
}

void Correlator::writeColumnNames(std::vector<int> &indices)
{
	m_phifile << "#t";
	m_kernelfile << "#t";

	for (unsigned int i=0; i<indices.size(); ++i)
	{
		m_phifile << "\tPhi[" << indices[i] << "]";
		m_kernelfile << "\tm[" << indices[i] << "]";
	}

	m_phifile << "\nt";
	m_kernelfile << "\nt";

	for (unsigned int i=0; i<indices.size(); ++i)
	{
		m_phifile << "\t" << m_qmesh.getq(indices[i]);
		m_kernelfile << "\t" << m_qmesh.getq(indices[i]);
	}

	m_phifile << "\n";
	m_kernelfile << "\n";
}

void Correlator::writeDatatoFile(std::vector<int> &indices, int skip)
{
	for (int i=m_tsteps/2; i<m_tsteps; i=i+skip)
	{
		m_phifile << m_time[i];
		m_kernelfile << m_time[i];

		for (unsigned int index=0; index<indices.size(); ++index)
		{
			m_phifile << "\t" << m_values[indices[index]][i];
			m_kernelfile << "\t" << m_kernel[indices[index]][i];
		}

		m_phifile << "\n";
		m_kernelfile << "\n";
	}
}

double Correlator::getDiffusionConstant() const
{
	return m_D0;
}

double Correlator::getFrictionCoefficient() const
{
	return m_nu;
}

bool Correlator::isNewtonian() const
{
	return m_isNewtonian;
}

const std::vector<double>& Correlator::getLowQKernelValues() const
{
	return m_kernel[0];
}

const std::vector<double>& Correlator::getFrequencies() const
{
	return m_Omega;
}

double Correlator::calculateKernel(int t, int q, const array3D &Vertex, const array2D &fk, const array2D &fp)
{
	double F = 0;
	const std::vector<double> &Grid = m_qmesh.getQgrid();
	const array2D &Bounds = m_qmesh.getBorderIndices();

	m_vertInt[0] = 0;

	for (int k=0; k<m_qsteps; ++k)
	{
		int lower = 0, upper = 0;

		if (k==q)
		{
			upper = Bounds[k][k];
		}

		else if (k>q)
		{
			upper = Bounds[q][k];
			lower = Bounds[k][q];
		}

		else
		{
			upper = Bounds[k][q];
			lower = Bounds[q][k];
		}

		double plower = fabs(Grid[q]-Grid[k]);
		double pupper = (Grid[q]+Grid[k] < Grid[m_qsteps-1]) ? Grid[q]+Grid[k] : Grid.back();

		m_vertInt[k] = 0;

		for (int p=lower+1; p<upper+1; ++p)
		{
			double left = 0, right = 0;
			double dq = 0;

			if (p>lower+1 && p<upper)
			{
				left = fp[p-1][t]*Vertex[q][k][p-lower-1];
				right = fp[p][t]*Vertex[q][k][p-lower];
				dq = Grid[p]-Grid[p-1];
			}

			else if (p==lower+1)
			{
				double phiIP = getInterpolationValue(plower, Grid[lower], fp[lower][t], Grid[lower+1], fp[lower+1][t]);
				left = phiIP*Vertex[q][k][0];
				right = fp[p][t]*Vertex[q][k][1];
				dq = Grid[p]-plower;
			}

			else if (p==upper)
			{
				double phiIP = getInterpolationValue(pupper, Grid[upper-1], fp[upper-1][t], Grid[upper], fp[upper][t]);
				left = fp[p][t]*Vertex[q][k][p-lower-1];
				right = phiIP*Vertex[q][k][p-lower];
				dq = pupper-Grid[p-1];

			}

			m_vertInt[k] += .5*dq*(left+right);
		}

		m_vertInt[k] *= fk[k][t];

		F += .5*(Grid[k]-Grid[k-1])*(m_vertInt[k]+m_vertInt[k-1]);
	}

	return F;
}

void Correlator::calcSolution(const double eps, const int itermax)
{
	const array3D &Vertex = m_qmesh.getCollectiveVertex();

	for (int i=m_tsteps/2; i<m_tsteps; ++i) 
	{
		m_time[i] = i*m_spacing;
		int lbar = i/2;

		for (int q=0; q<m_qsteps; ++q)
		{
			m_values[q][i] = m_values[q][i-1];
		}

		for (int q=1; q<m_qsteps; ++q)
		{
			double A = m_nu/(2*m_spacing*m_Omega[q]*m_Omega[q]);
			double B = m_isNewtonian/pow(m_spacing*m_Omega[q],2);
			double C = m_kernel[q][i-lbar]*m_values[q][lbar];

			C -= m_dP[q][1]*m_kernel[q][i-1];
			C -= m_dM[q][1]*m_values[q][i-1];
				
			for (int l=2; l<=lbar; ++l)
			{
				C += m_dP[q][l]*(m_kernel[q][i-l+1]-m_kernel[q][i-l]);
			}
			
			for (int l=2; l<=i-lbar; ++l)
			{
				C += m_dM[q][l]*(m_values[q][i-l+1]-m_values[q][i-l]);
			}

			C -= A*(4*m_values[q][i-1]-m_values[q][i-2]);
			C -= B*(5*m_values[q][i-1]-4*m_values[q][i-2]+m_values[q][i-3]);
	
			m_kernel[q][i] = calculateKernel(i,q,Vertex,m_values,m_values);
			int iter = 0;

			while (++iter<=itermax)
			{
				double phi_old = m_values[q][i];
				m_values[q][i] = (1-m_dP[q][1])*m_kernel[q][i]-C;
				m_values[q][i] /= 1+m_dM[q][1]+3*A+2*B;
				m_kernel[q][i] = calculateKernel(i,q,Vertex,m_values,m_values);

				if (fabs(phi_old-m_values[q][i])<eps)
				{
					break;
				}
			}

			m_dP[q][i] = .5*(m_values[q][i]+m_values[q][i-1]);
			m_dM[q][i] = .5*(m_kernel[q][i]+m_kernel[q][i-1]);	
		}

		m_values[0][i] = m_values[1][i];
		m_kernel[0][i] = calculateLowQKernel(i);
	}
}

Collective::Collective(int tsteps, Mesh &Q, double h, double vtherm, double nu) :
	Correlator(tsteps, Q, h, vtherm, nu)		
{
}

Collective::~Collective()
{
	m_phifile.close();
	m_kernelfile.close();

	std::cout << "# Closing all coherent correlator files...\n";
}

void Collective::initializeCorrelator()
{
	const array3D &Vertex = m_qmesh.getCollectiveVertex();

	for (int i=0; i<m_tsteps/2; ++i)
	{		
		for (int q=1; q<m_qsteps; ++q)
		{
			m_kernel[q][i] = calculateKernel(i,q,Vertex,m_values,m_values);
		}		
	}

	for (int q=0; q<m_qsteps; ++q)
	{
		m_Omega[q] *= m_qmesh.getq(q)/sqrt(m_qmesh.getSq(q));
	}

	m_kernel[0][0] = calculateLowQKernel(0);

	for (int i=1; i<m_tsteps; ++i)
	{
		m_kernel[0][i] = calculateLowQKernel(i);

		for (int q=1; q<m_qsteps; ++q)
		{
			m_dP[q][i] = .5*(m_values[q][i]+m_values[q][i-1]);
			m_dM[q][i] = .5*(m_kernel[q][i]+m_kernel[q][i-1]);
		}
	}

	std::cout << "# Initialization of coherent correlators and kernels complete...\n";
}

double Collective::calculateLowQKernel(int t)
{
	double kernel = 0;
	const std::vector<double> &Grid = m_qmesh.getQgrid();
	const array3D &Vertex = m_qmesh.getCollectiveVertex();

	double left = 0;

	for (int k=1; k<m_qsteps; ++k)
	{
		double right = Vertex[0][k][0]*pow(m_values[k][t],2);
		kernel += .5*(Grid[k]-Grid[k-1])*(left+right);
		left = right;
	}

	return kernel;
}

double Collective::calculateLowQPrefactor()
{
	double particleDensity = 6*m_qmesh.getPackingFraction()/Pi;
	const std::vector<double> &S = m_qmesh.getSq();
	const std::vector<double> &c = m_qmesh.getcq();

	return particleDensity*pow(S[0],3)*pow(c[0],2)*pow(2*Pi*m_D0/S[0],-3./2.)/16;
}

Tagged::Tagged(int tsteps, Mesh &Q, double h, double vtherm, double nu, const Correlator &C) :
	Correlator(tsteps, Q, h, vtherm, nu),
	m_DL(m_D0), m_collective(C)
{
}

Tagged::~Tagged()
{
	m_phifile.close();
	m_kernelfile.close();

	std::cout << "# Closing all incoherent correlator files...\n";
}

void Tagged::initializeCorrelator()
{
	const array3D &Vertex = m_qmesh.getTaggedVertex();
	const array2D &collectivePhi = m_collective.getArrayValues();

	for (int i=0; i<m_tsteps/2; ++i)
	{	
		for (int q=1; q<m_qsteps; ++q)
		{
			m_kernel[q][i] = calculateKernel(i,q,Vertex,collectivePhi,m_values);
		}	
	}

	for (int q=0; q<m_qsteps; ++q)
	{
		m_Omega[q] *= m_qmesh.getq(q);
	}

	m_kernel[0][0] = calculateLowQKernel(0);

	for (int i=1; i<m_tsteps; ++i)
	{
		m_kernel[0][i] = calculateLowQKernel(i);

		for (int q=1; q<m_qsteps; ++q)
		{
			m_dP[q][i] = .5*(m_values[q][i]+m_values[q][i-1]);
			m_dM[q][i] = .5*(m_kernel[q][i]+m_kernel[q][i-1]);
		}
	}

	std::cout << "# Initialization of incoherent correlators and kernels complete...\n";
}

void Tagged::calcSolution(const double eps, const int itermax)
{
	const array3D &Vertex = m_qmesh.getTaggedVertex();
	const array2D &collectivePhi = m_collective.getArrayValues();

	for (int i=m_tsteps/2; i<m_tsteps; ++i) 
	{
		m_time[i] = i*m_spacing;
		int lbar = i/2;

		for (int q=0; q<m_qsteps; ++q)
		{
			m_values[q][i] = m_values[q][i-1];
		}

		for (int q=1; q<m_qsteps; ++q)
		{
			double A = m_nu/(2*m_spacing*m_Omega[q]*m_Omega[q]);
			double B = m_isNewtonian/pow(m_spacing*m_Omega[q],2);
			double C = m_kernel[q][i-lbar]*m_values[q][lbar];

			C -= m_dP[q][1]*m_kernel[q][i-1];
			C -= m_dM[q][1]*m_values[q][i-1];
				
			for (int l=2; l<=lbar; ++l)
			{
				C += m_dP[q][l]*(m_kernel[q][i-l+1]-m_kernel[q][i-l]);
			}
			
			for (int l=2; l<=i-lbar; ++l)
			{
				C += m_dM[q][l]*(m_values[q][i-l+1]-m_values[q][i-l]);
			}

			C -= A*(4*m_values[q][i-1]-m_values[q][i-2]);
			C -= (5*m_values[q][i-1]-4*m_values[q][i-2]+m_values[q][i-3])*B;
	
			m_kernel[q][i] = calculateKernel(i,q,Vertex,collectivePhi,m_values);
			int iter = 0;

			while (++iter<=itermax)
			{
				double phi_old = m_values[q][i];
				m_values[q][i] = (1-m_dP[q][1])*m_kernel[q][i]-C;
				m_values[q][i] /= 1+m_dM[q][1]+3*A+2*B;
				m_kernel[q][i] = calculateKernel(i,q,Vertex,collectivePhi,m_values);

				if (fabs(phi_old-m_values[q][i])<eps)
				{
					break;
				}
			}

			m_dP[q][i] = .5*(m_values[q][i]+m_values[q][i-1]);
			m_dM[q][i] = .5*(m_kernel[q][i]+m_kernel[q][i-1]);	
		}

		m_values[0][i] = m_values[1][i];
		m_kernel[0][i] = calculateLowQKernel(i);
	}
}

double Tagged::calculateLowQKernel(int t)
{
	double kernel = 0;
	const std::vector<double> &Grid = m_qmesh.getQgrid();
	const array3D &SVertex = m_qmesh.getTaggedVertex();
	const array2D &collectivePhi = m_collective.getArrayValues();

	double left = 0;

	for (int k=1; k<m_qsteps; ++k)
	{
		double right = SVertex[0][k][0]*m_values[k][t]*collectivePhi[k][t];
		kernel += .5*(Grid[k]-Grid[k-1])*(left+right);
		left = right;
	}

	return kernel;
}

double Tagged::calculateLowQPrefactor()
{
	double particleDensity = 6*m_qmesh.getPackingFraction()/Pi;
	const std::vector<double> &S = m_qmesh.getSq();
	const std::vector<double> &c = m_qmesh.getcq();
	double collectiveD0 = m_collective.getDiffusionConstant(); 

	return 3*particleDensity*S[0]*pow(c[0],2)*pow(Pi,-3./2.)*pow(collectiveD0/S[0]+m_DL,-5./2.)/16;
}

void Tagged::writeColumnNames(std::vector<int> &indices)
{
	m_phifile << "#t";
	m_kernelfile << "#t";

	for (unsigned int i=0; i<indices.size(); ++i)
	{
		m_phifile << "\tPhi[" << indices[i] << "]";
		m_kernelfile << "\tq^2 m[" << indices[i] << "]";
	}

	m_phifile << "\nt";
	m_kernelfile << "\nt";

	for (unsigned int i=0; i<indices.size(); ++i)
	{
		m_phifile << "\t" << m_qmesh.getq(indices[i]);
		m_kernelfile << "\t" << m_qmesh.getq(indices[i]);
	}

	m_phifile << "\n";
	m_kernelfile << "\n";
}

void Tagged::writeDatatoFile(std::vector<int> &indices, int skip)
{
	const std::vector<double> &Grid = m_qmesh.getQgrid();

	for (int i=m_tsteps/2; i<m_tsteps; i=i+skip)
	{
		m_phifile << m_time[i];
		m_kernelfile << m_time[i];

		for (unsigned int index=0; index<indices.size(); ++index)
		{
			m_phifile << "\t" << m_values[indices[index]][i];

			if (indices[index]!=0)
			{	
				m_kernelfile << "\t" << m_kernel[indices[index]][i]*pow(Grid[indices[index]],2);
			}
			else
			{
				m_kernelfile << "\t" << m_kernel[indices[index]][i];
			}
		}

		m_phifile << "\n";
		m_kernelfile << "\n";
	}
}

MSD::MSD(Mesh &Q, const Collective &A, const Tagged &B) :
	Decimation(A.getTimesteps(), 2, A.getSpacing()),
	m_qmesh(Q), m_diffusionConstant(m_tsteps, B.getDiffusionConstant()),
	m_nu(B.getFrictionCoefficient()), m_isNewtonian(B.isNewtonian()),
	m_phiA(A.getArrayValues()), m_phiB(B.getArrayValues()),
	m_kernelB(B.getLowQKernelValues()), m_file()
	
{
}

MSD::~MSD()
{
	m_file.close();
	std::cout << "# Closing MSD data file...\n";
}

void MSD::initialize()
{
	double D0 = m_diffusionConstant[0];
	calculateDiffusionConstant();

	for (int i=0; i<m_tsteps/2; ++i)
	{		
		m_values[0][i] = 6*D0*m_time[i];
		m_kernel[0][i] = m_kernelB[i]/3;
		m_kernel[1][i] = m_kernelB[i]/3;
		m_values[1][i] = m_nu*D0/3;
	}
	
	for (int i=1; i<m_tsteps; ++i)
	{
		for (int q=0; q<m_qsteps; ++q)
		{
			m_dP[q][i] = .5*(m_values[q][i]+m_values[q][i-1]);
			m_dM[q][i] = .5*(m_kernel[q][i]+m_kernel[q][i-1]);
		}
	}
}

void MSD::openFile(std::string filename)
{
	m_file.open(filename);

	std::cout << "# Opening file " << filename << "...\n";

	writeFileHeader();
	writeColumnNames();
}

void MSD::writeFileHeader()
{
	m_file << "# eta = " << m_qmesh.getPackingFraction();
	m_file << ", tau = " << m_qmesh.getTau();
	m_file << ", delta = " << m_qmesh.getDelta() << "\n";
	m_file << "# N = " << m_tsteps << "\n";
	m_file << "# h_initial = " << m_spacing;
	m_file << ", D0 = " << m_diffusionConstant[0] << "\n";
}

void MSD::writeColumnNames()
{
	m_file << "#t\tMSD\tMSD'\tD\tmMSD\tVACF\n";
}

void MSD::writeDatatoFile(int skip)
{
	for (int i=m_tsteps/2; i<m_tsteps; i+=skip)
	{
		m_file << m_time[i] << "\t" << m_values[0][i] << "\t";
		m_file << 6*m_diffusionConstant[i] << "\t";
		m_file << m_diffusionConstant[i] << "\t";
		m_file << m_kernel[0][i] << "\t";
		m_file << m_values[1][i] << "\n";
	}
}

double MSD::calculateMSDKernel(int t)
{
	return m_kernelB[t]/3;
}

void MSD::calculateMSD()
{
	double D0 = m_diffusionConstant[0];

	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = i*m_spacing;
		int lbar = i/2;

		double C = m_kernel[0][i-lbar]*m_values[0][lbar];
		C -= m_dP[0][1]*m_kernel[0][i-1]+m_dM[0][1]*m_values[0][i-1];

		for (int l=2; l<=lbar; ++l)
		{
			C += m_dP[0][l]*(m_kernel[0][i-l+1]-m_kernel[0][i-l]);
		}
		
		for (int l=2; l<=i-lbar; ++l)
		{
			C += m_dM[0][l]*(m_values[0][i-l+1]-m_values[0][i-l]);
		}

		m_kernel[0][i] = calculateMSDKernel(i);
		C += m_kernel[0][i]*m_dP[0][1];
		C *= D0;

		C += (m_values[0][i-2]-4*m_values[0][i-1])/(2*m_spacing);

		m_values[0][i] = (6*D0-C)/(3./(2*m_spacing)+D0*m_dM[0][1]);
		//m_values[0][i] = (6*D0-C)/(3./(2*m_spacing)+D0*(m_dM[0][1]+m_kernel[0][1]));

		m_dP[0][i] = .5*(m_values[0][i-1]+m_values[0][i]);
		m_dM[0][i] = .5*(m_kernel[0][i-1]+m_kernel[0][i]);
	}

	calculateDiffusionConstant();
}

void MSD::calculateVACF()
{
	/*for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = i*m_spacing;

		double A = m_isNewtonian*(4*m_values[1][i-1]-m_values[1][i-2])/(2*m_spacing);

		double B = 0;

		for (int l=2; l<i; ++l)
		{
			B += m_dM[1][l]*m_dP[1][i-l+1];
		}

		B *= m_spacing;

		double C = m_isNewtonian*3/(2*m_spacing)+m_nu+m_dM[1][1]*m_spacing;

		m_values[1][i] = (A-B)/C;

		m_dP[1][i] = .5*(m_values[1][i-1]+m_values[1][i]);
		m_dM[1][i] = .5*(m_kernel[1][i-1]+m_kernel[1][i]);

	}*/

	// via equation of motion
	/*double D0 = m_diffusionConstant[0];

	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = i*m_spacing;
		int lbar = i/2;

		m_kernel[1][i] = calculateMSDKernel(i);
		double C = m_kernel[1][i]*m_dP[1][1];
		C += m_kernel[1][i-lbar]*m_values[1][lbar];

		for (int l=2; l<=lbar; ++l)
		{
			C += m_dP[1][l]*(m_kernel[1][i-l+1]-m_kernel[1][i-l]);
		}
		
		for (int l=2; l<=i-lbar; ++l)
		{
			C += m_dM[1][l]*(m_values[1][i-l+1]-m_values[1][i-l]);
		}

		C -= m_dP[1][1]*m_kernel[1][i-1];
		C -= m_dM[1][1]*m_values[1][i-1];
		C *= D0;

		double dm = (3*m_kernel[1][i]-4*m_kernel[1][i-1]+m_kernel[1][i-2])/(2*m_spacing);

		m_values[1][i] = 6*D0*D0*dm-C;
		m_values[1][i] /= (1+D0*(m_dM[1][1]+m_kernel[1][1]));

		m_dP[1][i] = .5*(m_values[1][i-1]+m_values[1][i]);
		m_dM[1][i] = .5*(m_kernel[1][i-1]+m_kernel[1][i]);
	}*/

	/*for (int i=m_tsteps/2; i<m_tsteps-1; ++i)
	{
		double diff1 = 6*m_diffusionConstant[i+1]-6*m_diffusionConstant[i];
		double dt1 = m_time[i+1]-m_time[i];
		double diff2 = 6*m_diffusionConstant[i]-6*m_diffusionConstant[i-1];
		double dt2 = m_time[i]-m_time[i-1];

		m_values[1][i] = (diff1/dt1+diff2/dt2)/12;
	}

	m_values[1][m_tsteps-1] = m_values[1][m_tsteps-2];*/

	// via 2nd derivative of MSD
	for (int i=m_tsteps/2; i<m_tsteps-2; ++i)
	{
		double dt1 = m_time[i+1]-m_time[i];
		double dt2 = m_time[i]-m_time[i-1];

		m_values[1][i] = -1/12*m_values[0][i-2]+4/3*m_values[0][i-1]-5/2*m_values[0][i]+4/3*m_values[0][i+1]-1/12*m_values[0][i+2];
		m_values[1][i] /= dt1*dt2;
	}

	m_values[1][m_tsteps] = m_values[1][m_tsteps-1];
	m_values[1][m_tsteps] = m_values[1][m_tsteps-2];

	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_dP[0][i] = .5*(m_values[0][i-1]+m_values[0][i]);
		m_dM[0][i] = .5*(m_kernel[0][i-1]+m_kernel[0][i]);
	}
}

const std::vector<double>& MSD::getDiffusionConstant()
{
	return m_diffusionConstant;
}

double MSD::getDiffusionConstant(int index)
{
	return m_diffusionConstant[index];
}

void MSD::calculateDiffusionConstant()
{
	for (int i=1; i<m_tsteps-1; ++i)
	{
		double diff1 = m_values[0][i+1]-m_values[0][i];
		double dt1 = m_time[i+1]-m_time[i];
		double diff2 = m_values[0][i]-m_values[0][i-1];
		double dt2 = m_time[i]-m_time[i-1];
		m_diffusionConstant[i] = (diff1/dt1+diff2/dt2)/12;
	}

	m_diffusionConstant.back() = m_diffusionConstant[m_tsteps-2];
}
