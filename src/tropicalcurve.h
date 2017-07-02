#ifndef TROPICALCURVE_H_
#define TROPICALCURVE_H_

/*
 * This is the header file for the tropical curve algorithm by Andrew Chan.
 * The original tropical curve algorithm [Bogart et al] is located bergman.cpp.
 */

#include "vektor.h"
#include "polynomial.h"

/**
 * The input must define a variety of dimension 1 in the torus. The output is a list rays whose union is the tropical
 * variety of I. It is allowed that the ideal is homogeneous. This routine does not compute multiplicities.
 *
 * If earlyExit is set to true, this procedure will only compute a single ray.
 */

IntegerVectorList tropicalCurve(PolynomialSet const &I, bool earlyExit=false);

#endif
