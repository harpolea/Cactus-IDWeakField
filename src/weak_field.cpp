#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

#include "weak_field.h"

using namespace std;

// Stolen from EinsteinInitialData/Meudon_Bin_NS/src/Bin_NS.cc

extern "C"
void IDWeakField_initialise (CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO ("Setting up weak field initial data");

    // From Meudon_Mag_NS

    // Defined constants
    CCTK_REAL const c_light = 299792458.0; // speed of light [m/s]
    CCTK_REAL const eps0    = 1 / (mu0 * pow(c_light,2));

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
    Weak_Field weak_ns (npoints, &xx[0], &yy[0], &zz[0], filename);

    CCTK_VInfo (CCTK_THORNSTRING, "mass [M_sun]:       %g", weak_ns.mass);
    CCTK_VInfo (CCTK_THORNSTRING, "radius [km]:        %g",

    assert (bin_ns.np == npoints);

    CCTK_INFO ("Filling in Cactus grid points");

#pragma omp parallel for
    for (int i=0; i<npoints; ++i) {

      alp[i] = mag_ns.nnn[i];

      gxx[i] = weak_ns.g_xx[i];
      gyy[i] = weak_ns.g_yy[i];
      gzz[i] = weak_ns.g_zz[i];

      rho[i] = weak_ns.nbar[i] / rho_unit;

      eps[i] = rho[i] * weak_ns.ener_spec[i] / ener_unit;

      vel[i          ] = weak_ns.u_euler_x[i] / vel_unit;
      vel[i+  npoints] = weak_ns.u_euler_y[i] / vel_unit;
      vel[i+2*npoints] = weak_ns.u_euler_z[i] / vel_unit;

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

    CCTK_INFO ("Setting time derivatives of lapse and shift");
    {
      // These initial data assume stationarity

      if (CCTK_EQUALS (initial_dtlapse, "IDWeak_Field")) {
 #pragma omp parallel for
        for (int i=0; i<npoints; ++i) {
          dtalp[i] = 0.0;
        }
      } else if (CCTK_EQUALS (initial_dtlapse, "none")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }

      if (CCTK_EQUALS (initial_dtshift, "IDWeak_Field")) {
 #pragma omp parallel for
        for (int i=0; i<npoints; ++i) {
          dtbetax[i] = 0.0;
          dtbetay[i] = 0.0;
          dtbetaz[i] = 0.0;
        }
      } else if (CCTK_EQUALS (initial_dtshift, "none")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }
    }


    CCTK_INFO ("Done.");
    } catch (ios::failure e) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not read initial data from file '%s': %s", filename, e.what());
    }
}
