//================================================================
// $Id: PublicID.hpp,v 1.1 2008-07-08 13:58:41 dneilsen Exp $
//================================================================
///
/// \file
/// Defines the interface to interpolate BBH initial data.

#ifndef PublicID_hpp
#define PublicID_hpp

#include <vector>

/// Read the initial-data set from disk and stores it in memory.
/// Reads the spectral representation of the spatial metric, extrinsic
/// curvature, lapse and shift into memory and stores it.  The data is
/// read from the \em current directory.  This function must be called
/// precisely once before InterpolateData() is called.  The number
/// Omega is the orbital corotational frequency, which is given in the
/// file 'Omega' accompanying each data-set. If the wrong 'Omega' is
/// given, then the gauge will not correspond to the corotating
/// coordinate system.  Omega=0 corresponds to the inertial coordinate
/// system.
void ReadData(const double Omega);

/// Interpolate the data to a set of Cartesian coordinates (x[i],y[i],z[i]).
/// The Cartesian components of the metric g_{ij} and extrinsic curvature
/// K_{ij} (LOWER indices), as well as shift beta^i (UPPER index) and 
/// the Lapse N are returned in the corresponding arrays.
///
/// NOTE: This function can be called arbitrarily often. 
/// There is a modest per-call computational overhead, and a fairly
/// significant memory requirement per each interpolated point.
/// Therefore, as a guide, one should interpolate between 100 and 100000 
/// coordinate points in each call to 'InterpolateData'. 
void InterpolateData(const std::vector<double>& x,
		     const std::vector<double>& y,
		     const std::vector<double>& z,

		     std::vector<double>& gxx,
		     std::vector<double>& gxy,
		     std::vector<double>& gxz,
		     std::vector<double>& gyy,
		     std::vector<double>& gyz,
		     std::vector<double>& gzz,

		     std::vector<double>& Kxx,
		     std::vector<double>& Kxy,
		     std::vector<double>& Kxz,
		     std::vector<double>& Kyy,
		     std::vector<double>& Kyz,
		     std::vector<double>& Kzz,

		     std::vector<double>& Betax,
		     std::vector<double>& Betay,
		     std::vector<double>& Betaz,

		     std::vector<double>& N);

/// This releases the data being stored in memory and frees the
/// memory.  It's polite to tidy up when one is done, isn't it?
void ReleaseData(); 



/// PARALLELIZATION NOTE
/// 
/// The functions ReadData(), InterpolateData() and ReleaseData()
/// behave serially, and do not perform any MPI-communication.  If
/// desired, interpolation can be done in parallel for speed-up gains.
/// To do so, every processor must call ReadData() precisely once.
/// Then every processor will store its local copy of all data.  (a
/// full dataset is usually less than 200MB, so memory is not an
/// issue).  Subsequently, InterpolateData() can be called on each
/// processor independently with different coordinates (x[i],y[i],z[i]).


#endif

