#ifndef SINGULARCONVERSION_H_INCLUDED
#define SINGULARCONVERSION_H_INCLUDED

#include <assert.h>

#define OM_NDEBUG
#define NDEBUG
#include "singular/mod2.h"
#include "singular/structs.h" // Singular structs
#include "singular/ring.h"
#include "singular/numbers.h"
#include "singular/polys.h"
#include "singular/longrat.h"
#include "singular/ideals.h"
#include "singular/kstd1.h"
#include "singular/options.h"

//#include <libsingular.h>
#undef NDEBUG


#include "polynomialring.h"
#include "polynomial.h"
#include "field_rationals.h"
#include "printer.h"
#include "log.h"
#include <iostream>

ring singularRing(PolynomialRing const &r);
void freeSingularRing(ring R);
poly singularPolynomial(Polynomial const &p);
ideal singularPolynomialSet(PolynomialSet const &g);
FieldElement fromSingularCoefficient(PolynomialRing const &r, number c, const ring rSing);
Polynomial fromSingularPolynomial(PolynomialRing const &r, poly &p, const ring rSing);
PolynomialSet fromSingularIdeal(PolynomialRing const &r, ideal i, const ring rSing);

#endif
