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
    //CCTK_REAL const ener_unit  = pow(cactusL,2);           // from c^2
    CCTK_REAL const vel_unit   = cactusL / cactusT / c_light; // from c

// Problem here: vel_unit is identically 1, as it goes c->c.

    CCTK_INFO ("Setting up coordinates");

    int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    vector<double> xx(npoints), yy(npoints), zz(npoints);

#pragma omp parallel for
    for (int i=0; i<npoints; ++i) {
        xx[i] = x[i];// / cactusL;
        yy[i] = y[i];// / cactusL;
        zz[i] = z[i];// / cactusL;
    }

    // I hope that the coordinates have been ordered logically
    CCTK_INT const size = 3;
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing;
    GetDomainSpecification
      (size,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       & spacing);
    double xmin, xmax, zmin, zmax;

    xmin = physical_min[0];// / cactusL;
    xmax = physical_max[0];// / cactusL;
    zmin = physical_min[2];// / cactusL;
    zmax = physical_max[2];// / cactusL;

    cout << "\n\n\n" << xmax << "\n\n";

    try {
    CCTK_VInfo (CCTK_THORNSTRING, "mass [M_sun]:       %g", mass);
    CCTK_VInfo (CCTK_THORNSTRING, "radius [km]:        %g", radius);

    double RR = radius / coord_unit;
    double zcntr = 0.5 * (zmax + zmin);
    double z_smooth = 0.04 * (zmax - zmin);

    // convert stuff to Cactus units
    double _rho0 = rho0 / rho_unit;
    double _rho1 = rho1 / rho_unit;
    double _rho2 = rho2 / rho_unit;
    double _kh_u1 = kh_u1 / vel_unit;
    double _kh_u2 = kh_u2 / vel_unit;

    CCTK_INFO ("Filling in Cactus grid points");

#pragma omp parallel for
    for (int i=0; i<npoints; ++i) {

      double rr = RR + sqrt(xx[i]*xx[i] + yy[i]*yy[i] + zz[i]*zz[i]);

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

      double g = mass / rr;

      // y velocity zero everywhere
      vel[i+  npoints] = 0.0;

      if (CCTK_EQUALS (initial_hydro, "bubble")) {

          rho[i] = _rho0 * exp(-g * zz[i] / (eos_gamma * RR * alp[i]*alp[i]));

          eps[i] = pow(rho[i], eos_gamma - 1.0) / (eos_gamma - 1.0);

          double r_coord = sqrt(pow(xx[i]-bubble_x_pert*(xmax-xmin),2) + pow(zz[i]-bubble_z_pert*(zmax-zmin),2));

          if (r_coord <= bubble_r_pert)
          {
              eps[i] += eps[i] * (bubble_amp - 1.0) * 0.5 * (1.0 + tanh((2.0 - r_coord/(0.9 * bubble_r_pert))));
              rho[i] = pow(rho[i], eos_gamma) / (eps[i] * (eos_gamma - 1.0));
          }

          // velocity zero everywhere
          vel[i          ] = 0.0;
          vel[i+2*npoints] = 0.0;

      } else if (CCTK_EQUALS (initial_hydro, "kh")) {

          if (zz[i] < zcntr) {
              rho[i] = _rho1 - (_rho2-_rho1) * exp((zz[i]-zcntr)/(0.025*(xmax-xmin)));

              // u
              vel[i] = _kh_u1 - (_kh_u2-_kh_u1) * exp((zz[i]-zcntr)/(0.025*(xmax-xmin)));


          } else {
              rho[i] = _rho2 + (_rho2-_rho1) * exp((-zz[i]+zcntr)/(0.025*(xmax-xmin)));

              // u
              vel[i] = _kh_u2 + (_kh_u2-_kh_u1) * exp((-zz[i]+zcntr)/(0.025*(xmax-xmin)));

          }

          // v
          vel[i+2*npoints] = 0.5 * _kh_u1 * sin(4.0 * M_PI * (xx[i] + 0.5 * (xmax-xmin))/(xmax-xmin));

          eps[i] = pow(rho[i], eos_gamma - 1.0) / (eos_gamma - 1.0);

      } else if (CCTK_EQUALS (initial_hydro, "rt")) {

          rho[i] = _rho1 + (_rho2 - _rho1) * 0.5 * (1.0 + tanh((zz[i]-zcntr) / 0.9 * z_smooth));

          rho[i] *= exp(-g * zz[i] / (eos_gamma * RR * alp[i]*alp[i]));

          vel[i          ] = 0.0;
          vel[i+2*npoints] = rt_amp * cos(2.0 * M_PI * xx[i] / (xmax-xmin)) * exp(-(zz[i]-zcntr)*(zz[i]-zcntr)/(rt_sigma*rt_sigma));

          eps[i] = pow(rho[i], eos_gamma - 1.0) / (eos_gamma - 1.0);

      } else {
          CCTK_WARN (CCTK_WARN_ABORT, "incorrect initial_hydro");
      }

      // Check density not too small
      if (rho[i] < 1.e-20) {
        rho[i          ] = 1.e-20;
        vel[i          ] = 0.0;
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
