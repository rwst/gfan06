/*
 * gebauermoeller.h
 *
 *  Created on: Apr 13, 2014
 *      Author: anders
 */

#ifndef GEBAUERMOELLER_H_
#define GEBAUERMOELLER_H_

#include "packedmonomial.h"

// Mora's version of Gebauer-Moeller. Source: Mora's book

typedef PackedMonomial<2> MyMonomial; // <--- If we want this to be dynamic, the code gets much more complicated
typedef vector<Polynomial> PolynomialVector;
class SPair
{
public:
  int i;
  int j;
  PolynomialVector const &basis;
  IntegerVector v;
  MyMonomial packedLCM;
  TermOrder const *SPairTermOrder; //Because SPairs can be stored in a set, they must have an operator< so we need to store the term order used for this comparison in the object
								   //It is not enough to have a thread local variable, as Buchberger calls can be nested. Of course a thread local stack of these variables could work.
  void computePacked(PacMan const &pacman)const
  {
	  const_cast<MyMonomial&>(packedLCM)=MyMonomial(0,v,pacman);
  }
  SPair(int i_,int j_, PolynomialVector const &basis_, TermOrder const *SPairTermOrder_):
    i(i_),
    j(j_),
    basis(basis_),
    SPairTermOrder(SPairTermOrder_)
  {
    v=max(basis[i].getMarked().m.exponent,basis[j].getMarked().m.exponent);
  }
  SPair(int i_,int j_, PolynomialVector const &basis_,PacMan const &pacman, TermOrder const *SPairTermOrder_):
	    i(i_),
	    j(j_),
	    basis(basis_),
	    SPairTermOrder(SPairTermOrder_)
  {
	    v=max(basis[i].getMarked().m.exponent,basis[j].getMarked().m.exponent);
	  computePacked(pacman);
  }
  Polynomial sPolynomial_()const{return sPolynomial(basis[i],basis[j]);};
  bool relativelyPrime_()const{return relativelyPrime(basis[i].getMarked().m.exponent,basis[j].getMarked().m.exponent);}
  bool operator<(const class SPair &a)const;
  void print(Printer &p)const
  {
    p.printString("S");
    p.printInteger(i+1);//debug
    p.printString(",");
    p.printInteger(j+1);//debug
    p.printString("(");
    p.printPolynomial(basis[i]);
    p.printString(",");
    p.printPolynomial(basis[j]);
    p.printString(")=");
    p.printPolynomial(sPolynomial_());
    p.printString(" ");
    p.printVector(v);
    p.printString("\n");
  };
};

bool SPair::operator<(const class SPair &a)const //partial order
{
	if(v.sum()<a.v.sum())return true;
	if(v.sum()>a.v.sum())return false;
//	if((*SPairTermOrder)(a.v,v))return true;
//	if((*SPairTermOrder)(v,a.v))return false;
	if((*SPairTermOrder)(v,a.v))return true;
	if((*SPairTermOrder)(a.v,v))return false;
//	if(LexicographicTermOrder()(v,a.v))return true;//change
//	if(LexicographicTermOrder()(a.v,v))return false;//chage



	if(j<a.j)return true;
      if(j>a.j)return false;
      if(i<a.i)return true;
      if(i>a.i)return false;
      assert(0);//<--maybe this is not an error
      return false;
}

typedef set<SPair> SPairSet;

void printSPairSet(SPairSet const &l)
{
  AsciiPrinter p(Stderr);
  for(SPairSet::const_iterator i=l.begin();i!=l.end();i++)i->print(p);
}




class SPairContainer
{
public:
	SPairContainer(TermOrder const &to){}
	void updatePairs(SPairSet &sPolynomials, PolynomialVector const &G, int s);
	bool isEmpty()const;
	pair<int, int> popPair();
};

#if 0
class SPairContainerSimple: public SPairContainer
{
	TermOrder const *SPairTermOrder;
	set<SPair> sPolynomials;
public:
	SPairContainerSimple(TermOrder const &to):
		SPairContainer(to)
	{
		SPairTermOrder=&to;
		::SPairTermOrder=&to;//updating global!
	}
	void updatePairs(PolynomialVector const &G, int s)
	{
		//cerr<<"a";
	  for(SPairSet::iterator ij=sPolynomials.begin();ij!=sPolynomials.end();) // TODO: Rewrite using "remove if"
	    {
	      bool skip=false;
	      if(G[s].getMarked().m.exponent.divides(ij->v))
	    	  if(ij->v!=max(G[s].getMarked().m.exponent,G[ij->i].getMarked().m.exponent))
	    		  if(ij->v!=max(G[s].getMarked().m.exponent,G[ij->j].getMarked().m.exponent))
	    		  {
	    			  SPairSet::iterator temp=ij;
	    			  ij++;skip=true;
	    			  sPolynomials.erase(temp);
	    		  }
	      if(!skip)ij++;
	    }

	  SPairSet S;
	  for(int i=0;i<s;i++)
	  {
		  bool found=false;
		  for(int j=1;j<s;j++)
		  {
			  IntegerVector IS=max(G[i].getMarked().m.exponent,G[s].getMarked().m.exponent);
			  if(G[j].getMarked().m.exponent.divides(IS) &&(IS!=
					  (max(G[j].getMarked().m.exponent,G[s].getMarked().m.exponent)))){found=true;break;}
		  }
		  if(!found)
			  S.insert(SPair(i,s,G));
	  }

	  set<IntegerVector> T;//TODO: split S according to V value in some clever way
	  for(SPairSet::const_iterator i=S.begin();i!=S.end();i++)T.insert(i->v);

	  for(set<IntegerVector>::const_iterator tau=T.begin();tau!=T.end();tau++)
	  {
		  int choseni=0;
		  int chosens=0;
		  bool add=true;
		  for(SPairSet::const_iterator is=S.begin();is!=S.end();is++)
			  if(is->v==*tau)
			  {
				  choseni=is->i;
				  chosens=is->j;
	//			  debug<<(G[is->i].getMarked().m.exponent)<<G[is->j].getMarked().m.exponent<<"\n";
				  if(relativelyPrime(G[is->i].getMarked().m.exponent,G[is->j].getMarked().m.exponent))
				  {
					  add=false;
					  break;
				  }
			  }
		  if(add)
		  {
	//		  cerr<<"adding"<<choseni<<chosens<<endl;
			  sPolynomials.insert(SPair(choseni,chosens,G));
		  }
	  }
	//  cerr<<"b\n";
	}

	bool isEmpty()const
	{
		return sPolynomials.empty();
	}
	pair<int, int> popPair()
	{
		int i=sPolynomials.begin()->i;
		int j=sPolynomials.begin()->j;
	      sPolynomials.erase(sPolynomials.begin());
	      return pair<int,int>(i,j);
		}
	int size()const
	{return sPolynomials.size();}
};
#endif

class SPairContainerOptimized: public SPairContainer
{
	TermOrder const *SPairTermOrder;
	set<SPair> sPolynomials;
	vector<IntegerVector> exponentVectors;
public:
	SPairContainerOptimized(PolynomialRing const &r, TermOrder const &to):
		SPairContainer(to)
	{
		SPairTermOrder=&to;
	//	::SPairTermOrder=&to;//updating global!
	}
	void updatePairs(PolynomialVector const &G, int s)
	{
		exponentVectors.push_back(G[s].getMarked().m.exponent);
		//cerr<<"a";
	  for(SPairSet::iterator ij=sPolynomials.begin();ij!=sPolynomials.end();) // TODO: Rewrite using "remove if"
	    {
	      bool skip=false;
	      if(G[s].getMarked().m.exponent.divides(ij->v))
	    	  if(ij->v!=max(exponentVectors[s],exponentVectors[ij->i]))
	    		  if(ij->v!=max(exponentVectors[s],exponentVectors[ij->j]))
	    		  {
	    			  SPairSet::iterator temp=ij;
	    			  ij++;skip=true;
	    			  sPolynomials.erase(temp);
	    		  }
	      if(!skip)ij++;
	    }

	  map<IntegerVector,pair<int,int> > T;
	  for(int i=0;i<s;i++)
	  {
		  IntegerVector IS=max(exponentVectors[i],exponentVectors[s]);
		  bool found=false;
		  for(int j=1;j<s;j++)
		  {
			  if(exponentVectors[j].divides(IS) &&(IS!=
					  (max(exponentVectors[j],exponentVectors[s])))){found=true;break;}
		  }
		  if(!found)
		  {
			  if(T.count(IS)==1)
			  {
				  if(relativelyPrime(exponentVectors[i],exponentVectors[s]))
					  T[IS].first=-1;
			  }
			  else
			  {
				  if(relativelyPrime(exponentVectors[i],exponentVectors[s]))
					  T[IS]=pair<int,int>(-1,-1);
				  else
					  T[IS]=pair<int,int>(i,s);
			  }
		  }
	  }
	  for(map<IntegerVector,pair<int,int> >::const_iterator is=T.begin();is!=T.end();is++)
	  {
		  if(is->second.first!=-1)
		  {
			  int choseni=is->second.first;
			  int chosens=is->second.second;
			  sPolynomials.insert(SPair(choseni,chosens,G,SPairTermOrder));
		  }
	  }
	}

	bool isEmpty()const
	{
		return sPolynomials.empty();
	}
	pair<int, int> popPair()
	{
		int i=sPolynomials.begin()->i;
		int j=sPolynomials.begin()->j;
	      sPolynomials.erase(sPolynomials.begin());
	      return pair<int,int>(i,j);
		}
	int size()const
	{return sPolynomials.size();}
};


class SPairContainerOptimizedPacked: public SPairContainer
{
	TermOrder const *SPairTermOrder;
	set<SPair> sPolynomials;
	vector<IntegerVector> exponentVectors;
	vector<MyMonomial> packedMonomials;
	PacMan pacman;
	int n;
	PolynomialRing r;
public:
	SPairContainerOptimizedPacked(PolynomialRing const &r_, TermOrder const &to, PolynomialSet const *hint=0):
		SPairContainer(to),
		pacman(r_,vector<int>(r_.getNumberOfVariables(),0),0),n(r_.getNumberOfVariables()),
		r(r_)
	{
		SPairTermOrder=&to;
//		::SPairTermOrder=&to;//updating global!

		if(hint)
		{
			IntegerVector maxexp(r_.getNumberOfVariables());
			for(PolynomialSet::const_iterator j=hint->begin();j!=hint->end();j++)
				maxexp=max(maxexp,j->getMarked().m.exponent);
			vector<int> maxexp2;for(int i=0;i<maxexp.size();i++)maxexp2.push_back(maxexp[i]);
			pacman=PacMan(r_,maxexp2,0);
		}
	}
	void grow()
	{
//		cerr<<"GROWING\n";
		IntegerVector maxexp(n);
		for(vector<IntegerVector>::const_iterator i=exponentVectors.begin();i!=exponentVectors.end();i++)
			maxexp=max(maxexp,*i);
		vector<int> maxexp2;for(int i=0;i<maxexp.size();i++)maxexp2.push_back(maxexp[i]);
		pacman=PacMan(r,maxexp2,0);
		packedMonomials=vector<MyMonomial>(0);
		for(vector<IntegerVector>::const_iterator i=exponentVectors.begin();i!=exponentVectors.end();i++)
			packedMonomials.push_back(MyMonomial(0,*i,pacman));
//		pacman.print();


		for(SPairSet::iterator ij=sPolynomials.begin();ij!=sPolynomials.end();ij++)
			ij->computePacked(pacman);
	}
	void updatePairs(PolynomialVector const &G, int s)
	{
		exponentVectors.push_back(G[s].getMarked().m.exponent);
		if(pacman.fits(G[s].getMarked().m.exponent))
			packedMonomials.push_back(MyMonomial(0,G[s].getMarked().m.exponent,pacman));
		else
		{
			grow();
		}
//		cerr<<"a";
	  for(SPairSet::iterator ij=sPolynomials.begin();ij!=sPolynomials.end();) // TODO: Rewrite using "remove if"
	    {
	      bool skip=false;
/*	      if(packedMonomials[s].divides(ij->packedLCM,pacman)!=G[s].getMarked().m.exponent.divides(ij->v))
	      {
	    	  debug<<"PACKED DIVISOR:\n";packedMonomials[s].print(pacman);debug<<"\n";
	    	  debug<<"LCM:\n";ij->packedLCM.print(pacman);debug<<"\n";
	    	  debug<<G[s].getMarked().m.exponent<<"\n";
	    	  debug<<ij->v<<"\n";
	    	  pacman.print();
	    	  assert(0);
	      }*/
	      if(packedMonomials[s].divides(ij->packedLCM,pacman))
	    	  //if(G[s].getMarked().m.exponent.divides(ij->v))
	    	  if(ij->v!=max(exponentVectors[s],exponentVectors[ij->i]))
	    		  if(ij->v!=max(exponentVectors[s],exponentVectors[ij->j]))
	    		  {
	    			  SPairSet::iterator temp=ij;
	    			  ij++;skip=true;
	    			  sPolynomials.erase(temp);
	    		  }
	      if(!skip)ij++;
	    }

	  map<IntegerVector,pair<int,int> > T;
	  for(int i=0;i<s;i++)
	  {
		  bool found=false;
		  for(int j=1;j<s;j++)
		  {
/*			  if(packedMonomials[j].dividesLCM(packedMonomials[i],packedMonomials[s],pacman)!=exponentVectors[j].divides(IS))
		      {
		    	  debug<<"J:\n";packedMonomials[j].print(pacman);debug<<"\n";
		    	  debug<<"I:\n";packedMonomials[i].print(pacman);debug<<"\n";
		    	  debug<<"S:\n";packedMonomials[s].print(pacman);debug<<"\n";
		    	  debug<<G[j].getMarked().m.exponent<<"\n";
		    	  debug<<G[i].getMarked().m.exponent<<"\n";
		    	  debug<<G[s].getMarked().m.exponent<<"\n";
		    	  pacman.print();
		    	  assert(0);
		      }*/
/*			  assert(packedMonomials[j].dividesLCM(packedMonomials[i],packedMonomials[s],pacman)&&packedMonomials[i].strictlyBiggerOnOneCoordinate(packedMonomials[j],packedMonomials[s],pacman)==
					  (exponentVectors[j].divides(IS)&&(IS!=
					  (max(exponentVectors[j],exponentVectors[s])))));
*/
/*			  if(packedMonomials[j].dividesLCM(packedMonomials[i],packedMonomials[s],pacman)&&packedMonomials[i].strictlyBiggerOnOneCoordinate(packedMonomials[j],packedMonomials[s],pacman)!=
					  (exponentVectors[j].divides(IS)&&(IS!=
					  (max(exponentVectors[j],exponentVectors[s])))))
{
	debug<<"Condition1: new"<<(int)packedMonomials[j].dividesLCM(packedMonomials[i],packedMonomials[s],pacman)<<" old"<<(int)exponentVectors[j].divides(IS)<<"\n";
	debug<<"Condition2: new"<<(int)packedMonomials[i].strictlyBiggerOnOneCoordinate(packedMonomials[j],packedMonomials[s],pacman)<<" old"<<(int)(IS!=(max(exponentVectors[j],exponentVectors[s])))<<"\n";
	  debug<<"J:\n";packedMonomials[j].print(pacman);debug<<"\n";
	  debug<<"I:\n";packedMonomials[i].print(pacman);debug<<"\n";
	  debug<<"S:\n";packedMonomials[s].print(pacman);debug<<"\n";
	  debug<<G[j].getMarked().m.exponent<<"\n";
	  debug<<G[i].getMarked().m.exponent<<"\n";
	  debug<<G[s].getMarked().m.exponent<<"\n";
	  pacman.print();
	  assert(0);
}*/

		//	  if(exponentVectors[j].divides(IS) &&(IS!=
		//			  (max(exponentVectors[j],exponentVectors[s]))))
				  if(packedMonomials[j].dividesLCM(packedMonomials[i],packedMonomials[s],pacman)&&packedMonomials[i].strictlyBiggerOnOneCoordinate(packedMonomials[j],packedMonomials[s],pacman))
				  {found=true;break;}
		  }
		  if(!found)
		  {
			  IntegerVector IS=max(exponentVectors[i],exponentVectors[s]);
			  if(T.count(IS)==1)
			  {
				  if(relativelyPrime(exponentVectors[i],exponentVectors[s]))
					  T[IS].first=-1;
			  }
			  else
			  {
				  if(relativelyPrime(exponentVectors[i],exponentVectors[s]))
					  T[IS]=pair<int,int>(-1,-1);
				  else
					  T[IS]=pair<int,int>(i,s);
			  }
		  }
	  }
	  for(map<IntegerVector,pair<int,int> >::const_iterator is=T.begin();is!=T.end();is++)
	  {
		  if(is->second.first!=-1)
		  {
			  int choseni=is->second.first;
			  int chosens=is->second.second;
			  sPolynomials.insert(SPair(choseni,chosens,G,pacman,SPairTermOrder));
		  }
	  }
	}

	bool isEmpty()const
	{
		return sPolynomials.empty();
	}
	pair<int, int> popPair()
	{
		int i=sPolynomials.begin()->i;
		int j=sPolynomials.begin()->j;
	      sPolynomials.erase(sPolynomials.begin());
	      return pair<int,int>(i,j);
		}
	int size()const
	{return sPolynomials.size();}
};

class SPairContainerPair
{
	SPairContainerOptimized A;
	SPairContainerOptimizedPacked B;
public:
	SPairContainerPair(PolynomialRing const &r, TermOrder const &to):
		A(r,to),B(r,to){}
	void updatePairs(PolynomialVector const &G, int s)
	{
		A.updatePairs(G,s);
		B.updatePairs(G,s);
	}
	bool isEmpty()const
	{
		bool a=A.isEmpty();
		bool b=B.isEmpty();
		assert(a==b);
		return a;
	}
	pair<int, int> popPair()
	{
		pair<int,int> a=A.popPair();
		pair<int,int> b=B.popPair();
		assert(a.first==b.first);
		assert(a.second==b.second);
		return a;
	}
	int size()const
	{return A.size();}//debug only
};

#endif /* GEBAUERMOELLER_H_ */
