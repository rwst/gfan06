/*
 * traverser_bsptree.h
 *
 *  Created on: Aug 23, 2011
 *      Author: anders
 */

#ifndef TRAVERSER_BSPTREE_H_INCLUDED
#define TRAVERSER_BSPTREE_H_INCLUDED

#include "symmetrictraversal.h"
#include "bsptree.h"

/**
 * This class finds the closures of the connected components of the complement of the union of the cones stored in the BSPTree.
 * TODO: remove this class and replace it by BSPTreeTraverser.
 */
  class BSPTreeTraverser: public ConeTraverser
  {
    BSPTree const&tree;
    PolyhedralCone theCone;
  public:
          BSPTreeTraverser(int n_, BSPTree const & tree_);
          virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
          virtual IntegerVectorList link(IntegerVector const &ridgeVector);
          PolyhedralCone & refToPolyhedralCone();
          virtual bool hasNoState()const;
  };

  /**
   * This class has the same functionality as BSPTreeTraverser, but with the advantage that the BSP tree need not have been build.
   * TODO: rename this class and header and source file.
   */
  class BSPTreeTraverser2: public ConeTraverser
  {
    BSPTree const&tree;
    PolyhedralCone theCone;
  public:
    /**
     * Set reconstructable to true if the BSPTree stores multiplicities allowing polytope reconstruction.
     */
    BSPTreeTraverser2(int n_, BSPTree const & tree_, bool vertexReconstructing_);
    virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
    virtual IntegerVectorList link(IntegerVector const &ridgeVector);
    PolyhedralCone & refToPolyhedralCone();
    virtual bool hasNoState()const;
	void calibrate(){orthantCalibrate();}
  };
#endif
