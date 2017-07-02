/*
 * traverser_resultantfan.h
 *
 *  Created on: Apr 29, 2011
 *      Author: anders
 */

#ifndef TRAVERSER_RESULTANTFAN_H_INCLUDED
#define TRAVERSER_RESULTANTFAN_H_INCLUDED

#include "symmetrictraversal.h"
#include "triangulation2.h"

IntegerMatrix cayleyConfiguration(PolynomialSet const &g);
IntegerMatrix cayleyConfiguration(list<IntegerVectorList> const &g, int d);
vector<pair<int,int> > polynomialSetToIntervals(PolynomialSet const &g);
set<set<int> > mixedCells(vector<pair<int,int> > const&intervals, Triangulation2 const &t, IntegerMatrix const &A, IntegerVector const &w);
vector<pair<int,int> > subsetToIntervals(PolynomialSet const &g, set<int> const &subset);
IntegerVectorList quickLink(IntegerMatrix const &m, set<int> const &cell, vector<pair<int,int> > const &intervals, IntegerVectorList const &linealitySpace);

class ResultantFanTraverser : public ConeTraverser
{
  Triangulation2 theTriangulation;
  PolyhedralCone theCone;
//  PolynomialSet g;
  IntegerVectorListList tuple;
  IntegerMatrix theConfiguration;
  int n,d;
public:
  ResultantFanTraverser(IntegerVectorListList const &tuple_, IntegerMatrix const &theConfiguration_);
//  ResultantFanTraverser(PolynomialSet const &g_, IntegerMatrix const &theConfiguration_);
  virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
  virtual IntegerVectorList link(IntegerVector const &ridgeVector);
  PolyhedralCone &refToPolyhedralCone();
};

#endif
