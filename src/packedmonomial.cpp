/*
 * packedmonomial.cpp
 *
 *  Created on: Feb 8, 2014
 *      Author: anders
 */


#include "packedmonomial.h"
#include "parser.h"

#include <vector>

using namespace std;

template <class MonomialType>
vector<MonomialType> minimized(vector<MonomialType> const &generators)
{
	vector<MonomialType> temp(generators);

	assert(0);
//	sort(temp.begin(),temp.end());
/*	g->sort(polynomialOrder(LexicographicTermOrder()));

	  for(PolynomialSet::const_iterator i=g->begin();i!=g->end();i++)
	    {
	      bool someDivides=false;
	      for(PolynomialSet::const_iterator j=ret.begin();j!=ret.end();j++)
		//      for(PolynomialSet::const_iterator j=g->begin();j!=i;j++) //changed Feb 2009
	        {
	          if(j->getMarked().m.exponent.divides(i->getMarked().m.exponent))
	            {
	              someDivides=true;
	            }
	        }
	      if(!someDivides)
	        ret.push_back(*i);
	    }

	  *g=ret;
*/
}



void packedTest()
{
  PolynomialRing R=StringParser("Q[x,y,z]").parsePolynomialRing();
  Polynomial g=StringParser("x*y-z").parsePolynomial(R);
  vector<int> v;v.push_back(30);v.push_back(20);v.push_back(3);
  PacMan man(R,v,12);
  PackedMonomial<2> pm(342,g.terms.begin()->first.exponent,man);
  man.print(







	);
	debug<<"---------------------------------------------------------------\n"
			"";
	pm.print(man);

	{
	  vector<PackedMonomial<2> > v;
	  v.push_back(pm);
	  minimized<PackedMonomial<2> >(v);
	}
}

