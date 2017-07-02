#ifndef SYMMETRICTRAVERSAL_H_INCLUDED
#define SYMMETRICTRAVERSAL_H_INCLUDED

#include "symmetriccomplex.h"
#include "polyhedralfan.h"

/*
 This file contains the generic algorithm for traversing a connected component of a pure fan up to symmetry.
 This will in time be the algorithm to use for all fan traversals which are not reverse search.
 */
class ConeTraverser
{
protected:
	bool vertexConstructing;
    IntegerVector currentVertexCoordinates;//coordinates of dual object
    int n;
public:
    bool isVertexConstructing()const
    {
    	return vertexConstructing;
    }
    ConeTraverser(int n_):vertexConstructing(false),currentVertexCoordinates(n_),n(n_)
    {
    }
	/**
	 * Go to the cone which is connected to the current facet through the ridge in direction ray.
	 * The "ridge" is a relative interior point of the ridge.
	 */
	virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector)=0;
#if 0
	/**
	 * Returns for the offset for the last "changeCone" in the dual polytope. This is valid if offsetCoordinatesKnown returns true.
	 */
	virtual IntegerVector coordinateOffset(){return IntegerVector();}
	/**
	 * Tells whether the coordinateOffset function will be able to return coordinate offset for a dual object.
	 */
	virtual bool coordinateOffsetsKnown(){return false;}
#endif
/**
 * Compute the link of the fan in the ridge given by the vector ridge IS THIS A FACET NORMAL OR AN INTERIOR POINT?
 * This gives a list of symmetry invariant points under the actions keeping the link fixed.
 */
	virtual IntegerVectorList link(IntegerVector const &ridgeVector)=0;
	virtual PolyhedralCone & refToPolyhedralCone()=0;
/**
 * If there is no cone state data for the traverser, half of the changeCone() calls can be avoided.
 * That this is a valid of optimization for the ConeTraverser is indicated returning true in the following function.
 */
	virtual bool hasNoState()const;
	/**
	 * Return coordinates of vertex
	 */
	virtual IntegerVector getCoordinates()const
	{
		assert(vertexConstructing);
		return currentVertexCoordinates;
	}
	/**
	 * This method figures out the correct coordinatisation of the dual vertex of the cone.
	 * It only makes sense for full-dimensional fans.
	 * This method is only to be called on vertex constructing traversers.
	 * The reason for having this method is that it is non-trivial making a recoordinatisation of the computed vertices
	 * after the enumeration is completed, if this enumeration was done up to symmetry.
	 * While the traverser may change state during this procedure, the traverser is brought back to the current cone
	 * after this calibration.
	 */
	virtual void calibrate()
	{
		assert(0);
	}
protected:
	/**
	 * Does a walk in the fan, ending at the target.
	 */
	void walkTo(IntegerVector const &target);
	/**
	 * This method is an implementation of calibrate() which places the polytope in the non-negative orthant, closest
	 * to the origin.
	 */
	void orthantCalibrate();
};

class SymmetricTarget
{
public:
//	virtual bool process(PolyhedralCone const &cone)=0;
	virtual bool process(ConeTraverser &traverser)=0;
};

class SymmetricTargetNull:public SymmetricTarget
{
public:
	bool process(ConeTraverser &traverser){return true;};
};

class SymmetricTargetCounterInterrupted:public SymmetricTarget
{
	SymmetricTarget &target;
	int64 counter;
	int numberOfCalls;
public:
	SymmetricTargetCounterInterrupted(SymmetricTarget &target_, int64 numberOfCalls_=-1):target(target_),counter(0),numberOfCalls(numberOfCalls_)
		{
		}
	virtual bool process(ConeTraverser &traverser)
	{
		if(numberOfCalls<0)return target.process(traverser);
		if(counter<numberOfCalls&&target.process(traverser))
		{
			counter++;
			return numberOfCalls<0 || counter<numberOfCalls;
		}
		return false;
	}
};

class SymmetricTargetFanBuilder : public SymmetricTarget
{
	PolyhedralFan coneCollection;
public:
	PolyhedralFan const &getFanRef(){return coneCollection;}
//	SymmetricComplex toSymmetricComplex()const;
	SymmetricTargetFanBuilder(int n, SymmetryGroup const &sym);
	/* return false to exit */
//	bool process(PolyhedralCone const &cone);
	bool process(ConeTraverser &traverser);
};

class SymmetricTargetVertexSetBuilder : public SymmetricTarget
{
	IntegerVectorList vertexSet;
public:
	IntegerVectorList getVertexSet();
	SymmetricTargetVertexSetBuilder(int n, SymmetryGroup const &sym);
	bool process(ConeTraverser &traverser);
};

void symmetricTraverse(ConeTraverser &traverser, SymmetricTarget &target, SymmetryGroup const *sym=0);

#endif
