// MSD.cpp :
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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include "constants.h"
#include "mdefs.h"
#include "MSD.h"
#include "Correlator.h"

MSD::MSD(Mesh &Q, const Correlator &A, const Correlator &B, std::string filename) :
	m_qmesh(Q), m_tsteps(A.getTimesteps()), m_spacing(A.getSpacing()),
	m_D0(A.getDiffusionConstant()), m_time(m_tsteps),
	m_MSD(m_tsteps), m_mMSD(m_tsteps), m_dMSD(m_tsteps), m_dmMSD(m_tsteps),
	m_phiA(A.getCorrelatorValues()), m_phiB(B.getCorrelatorValues()),
	m_filename(filename), m_file(m_filename.c_str(), std::ios::out | std::ios::trunc)
	
{
	std::cout << "# Opening file " << m_filename << "...\n";

	writeFileHeader();
	writeColumnNames();

	for (int i=0; i<m_tsteps; ++i)
	{
		m_time[i] = i*m_spacing;
		m_MSD[i] = 6*m_D0*m_time[i];
		m_mMSD[i] = calculateMSDKernel(i);
	}
	
	for (int i=1; i<m_tsteps; ++i)
	{
		m_dMSD[i] = .5*(m_MSD[i]+m_MSD[i-1]);
		m_dmMSD[i] = .5*(m_mMSD[i]+m_mMSD[i-1]);
	}
}

MSD::~MSD()
{
	m_file.close(); std::cout << "# Closing " << m_filename << "...\n";
}

void MSD::writeFileHeader()
{
	m_file << "# eta = " << m_qmesh.getPackingFraction();
	m_file << ", tau = " << m_qmesh.getTau();
	m_file << ", delta = " << m_qmesh.getDelta() << "\n";
	m_file << "# N = " << m_tsteps << "\n";
	m_file << "# h_initial = " << m_spacing << ", D0 = " << m_D0 << "\n";
}

void MSD::writeColumnNames()
{
	m_file << "#t\tMSD\tmMSD\tVACF\n";
}

void MSD::writeDatatoFile(int skip)
{
	for (int i=m_tsteps/2; i<m_tsteps; i=i+skip)
	{
		m_file << m_time[i] << "\t" << m_MSD[i] << "\t";
		m_file << m_mMSD[i] << "\t" << "\n";
	}
}

void MSD::decimate()
{
	for (int i=0; i<m_tsteps/2; ++i)
	{
		m_time[i] = m_time[2*i];

		m_MSD[i] = m_MSD[2*i];
		m_mMSD[i] = m_mMSD[2*i];

		m_dMSD[i] = .5*(m_dMSD[2*i]+m_dMSD[2*i-1]);
		m_dmMSD[i] = .5*(m_dmMSD[2*i]+m_dmMSD[2*i-1]);
	}

	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = 0;
		
		m_MSD[i] = 0;
		m_mMSD[i] = 0;

		m_dMSD[i] = 0;
		m_dmMSD[i] = 0;
	}

	m_spacing *= 2;
}

int MSD::getTimesteps()
{
	return m_tsteps;
}

double MSD::getSpacing()
{
	return m_spacing;
}

const std::vector<double>& MSD::getMSDValues()
{
	return m_MSD;
}

const std::vector<double>& MSD::getKernelValues()
{
	return m_mMSD;
}	

double MSD::calculateMSDKernel(int t)
{
	double FMSD = 0;
	const std::vector<double> &Grid = m_qmesh.getQgrid();
	const std::vector<double> &S = m_qmesh.getSq();
	const std::vector<double> &c = m_qmesh.getcq();
	double particleDensity = m_qmesh.getPackingFraction()*6/Pi;

	double left = S[1]*pow(c[1],2)*pow(Grid[1],4)*m_phiA[1][t]*m_phiB[1][t];
	FMSD += .5*Grid[1]*left;

	for (int k=2; k<m_qmesh.getSize()-1; k++)
	{
		double right = S[k]*pow(c[k],2)*pow(Grid[k],4)*m_phiA[k][t]*m_phiB[k][t];
		FMSD += .5*(Grid[k]-Grid[k-1])*(left+right);
		left = right;
	}

	FMSD *= particleDensity/(6*Pi*Pi);

	return FMSD;
}

void MSD::calculateMSD()
{
	for (int i=m_tsteps/2; i<m_tsteps; ++i)
	{
		m_time[i] = i*m_spacing;

		double C = m_mMSD[i-m_tsteps/2]*m_MSD[m_tsteps/2]-m_dMSD[1]*m_mMSD[i-1]-m_dmMSD[1]*m_MSD[i-1];

		for (int l=2; l<=m_tsteps/2; ++l)
		{
			C += m_dMSD[l]*(m_mMSD[i-l+1]-m_mMSD[i-l]);
		}

		for (int l=2; l<=i-m_tsteps/2; ++l)
		{
			C += m_dmMSD[l]*(m_MSD[i-l+1]-m_MSD[i-l]);
		}

		m_mMSD[i] = calculateMSDKernel(i);
		C += m_mMSD[i]*m_dMSD[1];
		C *= m_D0;

		C += (m_MSD[i-2]-4*m_MSD[i-1])/(2*m_spacing);

		m_MSD[i] = (6*m_D0-C)/(3/(2*m_spacing)+m_D0*m_dmMSD[1]);

		m_dMSD[i] = .5*(m_MSD[i-1]+m_MSD[i]);
		m_dmMSD[i] = .5*(m_mMSD[i-1]+m_mMSD[i]);
	}
}
