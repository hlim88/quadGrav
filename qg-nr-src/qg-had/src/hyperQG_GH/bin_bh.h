/*
 *  Definition of class Bin_BH (binary black hole exportation)
 *
 */

/*
 *   Copyright (c) 2001  Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef __BIN_BH_H_
#define __BIN_BH_H_ 

/*
 * $Id: bin_bh.h,v 1.12 2014/10/13 08:54:05 j_novak Exp $
 * $Log: bin_bh.h,v $
 * Revision 1.12  2014/10/13 08:54:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2014/10/06 15:13:25  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.10  2010/07/14 16:47:30  e_gourgoulhon
 * Corrected error in the documentation for K_xx, K_xy, etc...:
 * the components are the covariant ones, not the contravariant ones.
 *
 * Revision 1.9  2009/09/22 09:23:23  p_grandclement
 * forgot to commit
 *
 * Revision 1.8  2009/09/10 10:05:35  p_grandclement
 * slight change to read different mass BH
 *
 * Revision 1.7  2003/10/24 15:49:03  e_gourgoulhon
 * Updated documentation.
 *
 * Revision 1.6  2003/01/09 11:08:00  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.5  2002/03/20 08:27:45  e_gourgoulhon
 * Added the derivatives of Psi.
 *
 * Revision 1.4  2002/02/06 14:54:44  e_gourgoulhon
 * Update of bibliographical references
 *
 * Revision 1.3  2002/01/10 14:06:58  e_gourgoulhon
 * Modif commentaries.
 *
 * Revision 1.2  2001/12/19 10:14:31  e_gourgoulhon
 * Updated documentation
 *
 * Revision 1.1  2001/12/19 10:08:31  e_gourgoulhon
 * Exporting Lorene structures
 *
 *
 *
 * $Header: /cvsroot/Lorene/Export/C++/Include/bin_bh.h,v 1.12 2014/10/13 08:54:05 j_novak Exp $
 *
 */

// Headers C
#include <cstdio>

#include <iostream>
#include <fstream>
using namespace std ;

namespace Lorene {
/**
 * Binary black hole configuration on a Cartesian grid.
 *
 * A binary black hole system is constructed on a Cartesian grid from
 * data stored in a file resulting from a computation by Grandclement,
 * Gourgoulhon and Bonazzola, Phys. Rev. D 65, 044021 (2002). 
 *
 * All the quantities are in units derived from the length
 * scale defined by the coordinate radius $a$ of black hole 1 apparent horizon
 * (throat).
 *
 * Importation of Lorene data is performed by means of the constructor
 * {\tt Bin\_BH::Bin\_BH(int, const double*, const double*, const double*, const char*)}.
 * This constructor takes general arrays for the location of the Cartesian coordinates
 * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one. Note also
 * that these arrays are 1-D, as well as all the metric fields,
 * in order to be use with any ordering of the 3-D storage.
 *
 *  This class is very simple, with all data members being public.
 *  A typical example of use is the following one
 *
 *  \begin{verbatim}
 *	    // Define the Cartesian grid by means of the arrays xg, yg, zg:
 *	    for (int i=0; i<nb_points; i++) {
 *           xg[i] = ...
 *           yg[i] = ...
 *           zg[i] = ...
 *	    }
 *
 *	    // Read the file containing the spectral data and evaluate
 *	    //  all the fields on the Cartesian grid :
 *
 *	    Bin_BH binary_system(nb_points, xg, yg, zg, fill, datafile) ;
 *
 *	    // Extract what you need :
 *
 *	    double* gamma_xx = binary_system.g_xx ; // metric coefficient g_xx
 *
 *	    double* shift_x = binary_system.beta_x ; // x comp. of shift vector
 *
 *	    ...
 *
 *	    // Save everything in an ASCII file :
 *
 *	    ofstream file_ini("ini.d") ;
 *	    binary_system.save_form(file_ini) ;
 *	    file_ini.close() ;
 *
 *  \end{verbatim}
 *
 * @version #$Id: bin_bh.h,v 1.12 2014/10/13 08:54:05 j_novak Exp $#
 */

class Bin_BH {

    // Data :
    // -----
    public:
	/// Orbital angular velocity [unit: $a^{-1}$]
	double omega ;

	/** Distance between the coordinate centers of two black
	 *  holes [unit: $a$]
	 */
	double dist ;

	/** Coordinate radius of the apparent horizon (throat) of
	 *   black hole 2 [unit: $a$].
	 *   NB: The coordinate radius of black hole 1 is 1 by definition
	 *       of the length unit.
	 */
	double radius2 ;

	/// Total number of grid points
	int np ;
    
	/// 1-D array storing the values of coordinate x of the {\tt np} grid points [unit: $a$]
	double* xx ;
    
	/// 1-D array storing the values of coordinate y of the {\tt np} grid points [unit: $a$]
	double* yy ; 

	/// 1-D array storing the values of coordinate z of the {\tt np} grid points [unit: $a$]
	double* zz ; 
	
	/// Lapse function $N$ at the {\tt np} grid points (1-D array)
	double* nnn ; 
	
	/// Component $\beta^x$ of the shift vector of corotating coordinates [unit: $c$]
	double* beta_x ; 
	
	/// Component $\beta^y$ of the shift vector of corotating coordinates [unit: $c$]
	double* beta_y ; 
	
	/// Component $\beta^z$ of the shift vector of corotating coordinates [unit: $c$]
	double* beta_z ; 
	
	/// Metric coefficient $\gamma_{xx}$ at the grid points (1-D array)
	double* g_xx ; 

	/// Metric coefficient $\gamma_{xy}$ at the grid points (1-D array)
	double* g_xy ; 

	/// Metric coefficient $\gamma_{xz}$ at the grid points (1-D array)
	double* g_xz ; 

	/// Metric coefficient $\gamma_{yy}$ at the grid points (1-D array)
	double* g_yy ; 

	/// Metric coefficient $\gamma_{yz}$ at the grid points (1-D array)
	double* g_yz ; 

	/// Metric coefficient $\gamma_{zz}$ at the grid points (1-D array)
	double* g_zz ; 

	/// Component $K_{xx}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_xx ; 

	/// Component $K_{xy}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_xy ;

	/// Component $K_{xz}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_xz ;

	/// Component $K_{yy}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_yy ;

	/// Component $K_{yz}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_yz ;

	/// Component $K_{zz}$ of the extrinsic curvature at the grid points (1-D array) [unit: $c/a$]
	double* k_zz ;

        /// First derivative $\partial/\partial x$ of the conformal factor $\Psi$  [unit: $a^{-1}$]
        double* dpsi_x ;

        /// First derivative $\partial/\partial y$ of the conformal factor $\Psi$  [unit: $a^{-1}$]
        double* dpsi_y ;

        /// First derivative $\partial/\partial z$ of the conformal factor $\Psi$  [unit: $a^{-1}$]
        double* dpsi_z ;

        /// Second derivative $\partial^2/\partial x^2$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_xx ;

        /// Second derivative $\partial^2/\partial x\partial y$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_xy ;

        /// Second derivative $\partial^2/\partial x\partial z$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_xz ;

        /// Second derivative $\partial^2/\partial y^2$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_yy ;

        /// Second derivative $\partial^2/\partial y\partial z$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_yz ;

        /// Second derivative $\partial^2/\partial z^2$ of the conformal factor $\Psi$  [unit: $a^{-2}$]
        double* d2psi_zz ;

    // Constructors - Destructor
    // -------------------------
    public:
	/** Constructor from Lorene spectral data.
	 *
	 * This constructor takes general arrays {\tt xi, yi, zi}
	 * for the location of the Cartesian coordinates
	 * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one.
	 * These arrays are 1-D to deal with any ordering of a 3-D storage.
	 *
	 *  @param nbpoints [input] Total number of grid points
	 *  @param xi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate x of the grid points [unit: $a$]
	 *  @param yi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate y of the grid points [unit: $a$]
	 *  @param zi [input] 1-D array (size {\tt nbpoints}) storing the
	 *		values of coordinate z of the grid points [unit: $a$]
	 *  @param fill [input] sets how the hole "interiors" must be
	 *		filled: \\
	 *		fill = 0 : all the fields are set to zero \\
	 *		fill = 1 : the fields are extrapolated from their
	 *		values "outside" the holes, by means of
	 *		parabolas along radial directions	
	 *  @param filename [input] Name of the (binary) file containing the result
	 *		of a computation by means of the multi-domain 
	 *		spectral method. 
	 */
	Bin_BH(int nbpoints, const double* xi, const double* yi, 
	       const double* zi, int fill, const char* filename, bool mdiff=false) ;	
	

	/** Constructor from a binary file 
	 *   (previously created by {\tt save\_bin})
	 */
	Bin_BH(FILE* ) ; 

	/** Constructor from a formatted file 
	 *   (previously created by {\tt save\_form})
	 */
	Bin_BH(ifstream& ) ; 

	/// Destructor    		
	~Bin_BH() ;			
 

    // Memory management
    // -----------------
    private:
	
	/// Allocate the memory for the arrays g_ij, k_ij, etc...
	void alloc_memory() ; 
    
    // Outputs
    // -------
    public:
	/** Save in a binary file.
	 *  This file can be subsenquently read by the evolution code, 
	 *  or by the constructor {\tt Bin\_BH::Bin\_BH(FILE* )}.
	 */
	void save_bin(FILE* ) const ;	    
    
	/** Save in a formatted file.
	 *  This file can be subsenquently read by the evolution code, 
	 *  or by the constructor {\tt Bin\_BH::Bin\_BH(ifstream\& )}.
	 */
	void save_form(ofstream& ) const ;	    

	/// Display
	friend ostream& operator<<(ostream& , const Bin_BH& ) ;	

};

}
#endif
