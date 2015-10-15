#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

// if changed the config:
// make sim-config options=simfactory/mdb/optionlists/ubuntu.cfg
// otherwise:
// make sim
// mpirun -np 4 exe/cactus_sim arrangements/EinsteinInitialData/IDWeakField/par/weak_field.par

#include "weak_field.h"

using namespace std;

// Stolen from EinsteinInitialData/Meudon_Bin_NS/src/Bin_NS.cc

/** Constructor
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
Weak_Field::Weak_Field(int nbpoints, const double *xi, const double *yi, const double *zi, const char *filename)
{
    np = nbpoints;
    xx = new double[np];
    yy = new double[np];
    zz = new double[np];

    for (int i = 0; i < np; i++) {
        xx[i] = xi[i];
        yy[i] = yi[i];
        zz[i] = zi[i];
    }
    // extract data from file and initialise the other stuff
    //open file
    //ifstream inFile(filename);

    //if (!inFile)
    //{
    //    CCTK_WARN (CCTK_WARN_ABORT, "Weak field inputs file could not be opened");
    //}

    mass = 1.0;
    radius = 10.0;

    // metric stuff
    nnn = new double[np];
    g_xx = new double[np];
    g_yy = new double[np];
    g_zz = new double[np];

    // hydro stuff
    nbar = new double[np];
    ener_spec = new double[np];
    u_euler_x = new double[np];
    u_euler_y = new double[np];
    u_euler_z = new double[np];

    //inFile.close();
}

// Constructor from a binary file
Weak_Field::Weak_Field(FILE *file)
{
    return;
}

// Constructor from a formatted file
Weak_Field::Weak_Field(ifstream &ifs)
{
    return;
}

// Destructor
Weak_Field::~Weak_Field()
{
    delete [] xx;
    delete [] yy;
    delete [] zz;
    delete [] nnn;
    delete [] g_xx;
    delete [] g_yy;
    delete [] g_zz;
    delete [] nbar;
    delete [] ener_spec;
    delete [] u_euler_x;
    delete [] u_euler_y;
    delete [] u_euler_z;
}

void Weak_Field::save_bin(FILE *file) const
{
    return;
}

void Weak_Field::save_form(ofstream &ofs) const
{
    return;
}

ostream &operator<<(ostream &output, const Weak_Field &weak_field)
{
    return output;
}


extern "C"
void IDWeakField_initialise (CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO ("Setting up weak field initial data");

    // From Meudon_Mag_NS

    // Defined constants
    CCTK_REAL const c_light = 299792458.0; // speed of light [m/s]

    // Constants of nature (IAU, CODATA):
    CCTK_REAL const G_grav = 6.67428e-11; // gravitational constant [m^3/kg/s^2]
    CCTK_REAL const M_sun  = 1.98892e+30; // solar mass [kg]

    // Cactus units in terms of SI units:
    // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
    // for the magnetic field)
    CCTK_REAL const cactusM = M_sun;
    CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
    CCTK_REAL const cactusT = cactusL / c_light;

    // Other quantities in terms of Cactus units
    CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km
    CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3
    CCTK_REAL const ener_unit  = pow(cactusL,2);           // from c^2
    CCTK_REAL const vel_unit   = cactusL / cactusT / c_light; // from c

    CCTK_INFO ("Setting up coordinates");

    int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    vector<double> xx(npoints), yy(npoints), zz(npoints);

#pragma omp parallel for
    for (int i=0; i<npoints; ++i) {
        xx[i] = x[i] * coord_unit;
        yy[i] = y[i] * coord_unit;
        zz[i] = z[i] * coord_unit;
    }

    CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);

    try {
    Weak_Field weak_field (npoints, &xx[0], &yy[0], &zz[0], filename);

    CCTK_VInfo (CCTK_THORNSTRING, "mass [M_sun]:       %g", weak_field.mass);
    CCTK_VInfo (CCTK_THORNSTRING, "radius [km]:        %g", weak_field.radius);

    double mass = weak_field.mass / cactusM;
    double radius = weak_field.radius / coord_unit;

    assert (weak_field.np == npoints);

    CCTK_INFO ("Filling in Cactus grid points");

#pragma omp parallel for
    for (int i=0; i<npoints; ++i) {

      double rr = radius + sqrt(xx[i]*xx[i] + yy[i]*yy[i] + zz[i]*zz[i]);

      alp[i] = sqrt(1.0 - 2.0 * mass / rr);

      gxx[i] = sqrt(1.0 + 2.0 * mass / rr);
      gyy[i] = sqrt(1.0 + 2.0 * mass / rr);
      gzz[i] = sqrt(1.0 + 2.0 * mass / rr);

      gxy[i] = 0.0;
      gxz[i] = 0.0;
      gyz[i] = 0.0;

      kxx[i] = 0.0;
      kxy[i] = 0.0;
      kxz[i] = 0.0;
      kyy[i] = 0.0;
      kyz[i] = 0.0;
      kzz[i] = 0.0;


      if (CCTK_EQUALS (initial_hydro, "bubble")) {

          rho[i] = weak_field.nbar[i] / rho_unit;

          eps[i] = rho[i] * weak_field.ener_spec[i] / ener_unit;

          vel[i          ] = 0.0;
          vel[i+  npoints] = 0.0;
          vel[i+2*npoints] = 0.0;

      } else if (CCTK_EQUALS (initial_hydro, "kh")) {

          rho[i] = weak_field.nbar[i] / rho_unit;

          eps[i] = rho[i] * weak_field.ener_spec[i] / ener_unit;
          // check to see if above middle of domain
          vel[i          ] = kh_u1 / vel_unit;
          vel[i+  npoints] = 0.0;
          vel[i+2*npoints] = 0.0;

      } else if (CCTK_EQUALS (initial_hydro, "rt")) {

          rho[i] = weak_field.nbar[i] / rho_unit;

          eps[i] = rho[i] * weak_field.ener_spec[i] / ener_unit;

          vel[i          ] = 0.0;
          vel[i+  npoints] = 0.0;
          vel[i+2*npoints] = 0.0;

      } else {
          CCTK_WARN (CCTK_WARN_ABORT, "incorrect initial_hydro");
      }

      // Especially the velocity is set to strange values outside of the
      // matter region, so take care of this in the following way
      if (rho[i] < 1.e-20) {
        rho[i          ] = 1.e-20;
        vel[i          ] = 0.0;
        vel[i+  npoints] = 0.0;
        vel[i+2*npoints] = 0.0;
        eps[i          ] = 0.0;
      }

    } // for i


    CCTK_INFO ("Done.");
    } catch (ios::failure e) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not read initial data from file '%s': %s", filename, e.what());
    }
}
