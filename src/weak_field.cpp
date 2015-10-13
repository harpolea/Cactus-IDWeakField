#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "weak_field.h"

// How do I c++???

// Stolen from EinsteinInitialData/Meudon_Bin_NS/src/Bin_NS.cc

extern "C"
void IDWeakField_initialise (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up weak field initial data");

  // Then cheat by calling your own C++ routine from here.

  CCTK_INFO ("Done setting up weak field initial data");
}
