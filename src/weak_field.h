#ifndef __WEAK_FIELD_H_
#define __WEAK_FIELD_H_


#include <stdio.h>

#ifdef OBSOLETE_HEADERS

#include <iostream.h>
#include <fstream.h>

#else

#include <iostream>
#include <fstream>
using namespace std;
#endif

class Weak_Field {
    public:
        // mass of neutron star [unit: solar mass]
        double mass;

        // radius of NS [unit: km]
        double radius;

        // total number of grid points
        int np;

        /// 1-D array storing the values of coordinate x of the {\tt np} grid points [unit: km]
        double* xx;

        /// 1-D array storing the values of coordinate y of the {\tt np} grid points [unit: km]
        double* yy;

        /// 1-D array storing the values of coordinate z of the {\tt np} grid points [unit: km]
        double* zz;

        /// Lapse function $N$ at the {\tt np} grid points (1-D array)
        double* nnn;

        /// Metric coefficient $\gamma_{xx}$ at the grid points (1-D array)
        double* g_xx;

        /// Metric coefficient $\gamma_{yy}$ at the grid points (1-D array)
        double* g_yy;

    	/// Metric coefficient $\gamma_{zz}$ at the grid points (1-D array)
    	double* g_zz;

    	// Hydro components
        //----------------------

        /** Baryon density in the fluid frame at the {\tt np} grid points (1-D array)
        * [unit: ${\rm kg \, m}^{-3}$]
        */
        double* nbar;

    	/// Specific internal energy at the  {\tt np} grid points (1-D array) [unit: $c^2$]
    	double* ener_spec;

    	/** Component $U^x$ of the fluid 3-velocity with respect to the Eulerian
        * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
        */
        double* u_euler_x;

    	/** Component $U^y$ of the fluid 3-velocity with respect to the Eulerian
        * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
        */
        double* u_euler_y;

        /** Component $U^z$ of the fluid 3-velocity with respect to the Eulerian
        * observer, at the {\tt np} grid points (1-D array) [unit: $c$]
        */
        double* u_euler_z;

        // Constructors - Destructor
        // -------------------------

        /** Constructor from Lorene spectral data.
        *
        * This constructor takes general arrays {\tt xi, yi, zi}
        * for the location of the Cartesian coordinates
        * $(x, y, z)$, i.e. it does not assume that the grid is a uniform one.
        * These arrays are 1-D to deal with any ordering of a 3-D storage.
        *
        *  @param nbpoints [input] Total number of grid points
        *  @param xi [input] 1-D array (size {\tt nbpoints}) storing the
        *		values of coordinate x of the grid points [unit: km]
        *  @param yi [input] 1-D array (size {\tt nbpoints}) storing the
        *		values of coordinate y of the grid points [unit: km]
        *  @param zi [input] 1-D array (size {\tt nbpoints}) storing the
        *		values of coordinate z of the grid points [unit: km]
        *  @param filename [input] Name of the (binary) file containing the result
        *		of a computation by means of the multi-domain
        *		spectral method.
        */

        Weak_Field(int nbpoints, const double* xi, const double* yi, const double* zi, const char* filename);

        /** Constructor from a binary file
        *   (previously created by {\tt save\_bin})
        */
        Weak_Field(FILE* );

        /** Constructor from a formatted file
        *   (previously created by {\tt save\_form})
        */
        Weak_Field(ifstream& );

        /// Destructor
        ~Weak_Field();

        // Outputs
        //--------

        /** Save in a binary file.
        *  This file can be subsenquently read by the evolution code,
        *  or by the constructor {\tt Bin\_NS::Bin\_NS(FILE* )}.
        */
        void save_bin(FILE* ) const;

        /** Save in a formatted file.
        *  This file can be subsenquently read by the evolution code,
        *  or by the constructor {\tt Bin\_NS::Bin\_NS(ifstream\& )}.
        */
        void save_form(ofstream& ) const;

        /// Display
        friend ostream& operator<<(ostream& , const Weak_Field& );

    private:

        // Memory management
        //------------------
        /// Allocate the memory for the arrays g\_ij, k\_ij, etc...
        void alloc_memory();

};

#endif
