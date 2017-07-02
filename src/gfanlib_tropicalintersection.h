#include "gfanlib_matrix.h"
#include "gfanlib_zcone.h"
#include "gfanlib_zfan.h"

namespace gfan{

class HalfOpenZCone
{
	ZCone lifted;
	HalfOpenZCone(ZCone const &_lifted);
public:
	int ambientDimension()const{return lifted.ambientDimension()-1;}
	HalfOpenZCone fromZCone(ZCone const &C, TermOrder const &t);
	HalfOpenCone(int dimension_, ZMatrix const &equations, ZMatrix const &nonstrict, ZMatrix const &strict, bool findFacets=false, bool canonicalize=false);
	HalfOpenCone(int ambientDimension);//full space
	bool isEmpty();
	friend HalfOpenCone intersection(const HalfOpenZCone &a, const HalfOpenZCone &b, bool findFacets=false);
	friend bool haveEmptyIntersection(const HalfOpenZCone &a, const HalfOpenZCone &b);
	ZCone closure();
	bool contains(IntegerVector const &v)const;
#if 0
	void splitIntoRelativelyOpenCones(list<HalfOpenCone> &l);
	void print(class Printer &p)const;
	friend bool operator<(HalfOpenZCone const &a, HalfOpenZCone const &b);
	void canonicalize(){lifted.canonicalize();}
	  /**
	     Remove all coordinates from the space except those listed in chosen.
	   */
	  HalfOpenCone withChosenCoordinates(list<int> chosen)const;
	  HalfOpenCone rewrite(FieldMatrix const &A, list<int> nonPivots)const;
	  HalfOpenCone rewriteExpand(list<int> pivots, IntegerVectorList const &newEquations)const;
#endif
};

class HalfOpenZConeProcessor
{
public:
	virtual bool process(HalfOpenZCone const &c); // Return true to abort.
};

typedef std::vector<HalfOpenZCone> HalfOpenZConeVector;

void intersectHalfOpenFansUsingIntersectionTable(int ambientDimension,
									std::vector<std::vector<HalfOpenZCone> > const &newtonPolytopes,
									HalfOpenZConeProcessor *processor,
									HalfOpenZCone const *restrictingCone=0,
									class IntersectionTable); // if a cone in the intersection

void intersectHalfOpenFans(int ambientDimension,
									std::vector<std::vector<HalfOpenZCone> > const &newtonPolytopes,
									HalfOpenZConeProcessor *processor,
									HalfOpenZCone const *restrictingCone=0); // if a cone in the intersection

void intersectTropicalHypersurfaces(int ambientDimension,
									std::vector<ZMatrix> const &newtonPolytopes,
									HalfOpenZConeProcessor *processor,
									HalfOpenZCone const *restrictingCone=0); // if a cone in the intersection

std::vector<HalfOpenZCone> tropicalHypersurfaceIntersection(int ambientDimension,
									std::vector<ZMatrix> const &newtonpolytopes,
									HalfOpenZCone const *restrictingCone=0);

std::vector<HalfOpenZCone> decomposedTropicalHypersurface(ZMatrix const &newtonPolytope);

std::vector<HalfOpenZCone> intersection(std::vector<HalfOpenZCone> const &A, std::vector<HalfOpenZCone> const &B);

}
