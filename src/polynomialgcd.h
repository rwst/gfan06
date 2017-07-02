/*
 * polynomialgcd.h
 *
 *  Created on: Apr 26, 2014
 *      Author: anders
 */

#ifndef POLYNOMIALGCD_H_
#define POLYNOMIALGCD_H_

#include "polynomial.h"

Polynomial polynomialGCD(Polynomial const &a, Polynomial const &b);
Polynomial NonMonomialPolynomialGCDForZModP(PolynomialSet p);
Polynomial NonMonomialPolynomialGCDForZ(PolynomialSet p);
#endif /* POLYNOMIALGCD_H_ */
