#include "dimension.h"
#include "buchberger.h"
#include "log.h"
#include "printer.h"

const int bitsPerWord=64;
typedef int64 ExponentType;

ExponentType allOnes(int n)
{
	if(n<64)return (((int64) 1)<<n)-1;
	return -1;
}

int sum(ExponentType v)
{
	int64 t0=v;
	int64 t1=(t0&0x5555555555555555)+((t0>>1)&0x5555555555555555);
	int64 t2=(t1&0x3333333333333333)+((t1>>2)&0x3333333333333333);
	int64 t3=(t2&0x0f0f0f0f0f0f0f0f)+((t2>>4)&0x0f0f0f0f0f0f0f0f);
	int64 t4=(t3&0x00ff00ff00ff00ff)+((t3>>8)&0x00ff00ff00ff00ff);
	int64 t5=(t4&0x0000ffff0000ffff)+((t4>>16)&0x0000ffff0000ffff);
	int64 t6=(t5&0x00000000ffffffff)+((t5>>32)&0x00000000ffffffff);
	return t6;
}

static void rek64(ExponentType &ones, ExponentType &zeros, int nOnes, int nZeros, vector<ExponentType> vectors, int &best, int n)
{
	if(nOnes>best)best=nOnes;
	if(n-nZeros<best)return;
	if(nOnes+nZeros==n)return;

#if 1
	int index=0;
	for(int i=0;i<n;i++,index++)
	//	if((!ones[i])&&(!zeros[i]))break;
		if((!((((int64)1)<<i)&ones))&&(!((((int64)1)<<i)&zeros)))break;
#else
	int index=nOnes+nZeros;
#endif
	assert(index<n);

	ones|=((int64)1)<<index;
	bool good=true;
	for(vector<ExponentType>::const_iterator i=vectors.begin();i!=vectors.end();i++)
		if(((*i)&~ones)==0)
		{
			good=false;
			break;
		}

	if(good)
    {
      rek64(ones,zeros,nOnes+1,nZeros,vectors,best,n);
    }
	ones-=((int64)1)<<index;

	vector<ExponentType> vectorsSubset;
	vectorsSubset.reserve(vectors.size());

	for(vector<ExponentType>::const_iterator i=vectors.begin();i!=vectors.end();i++)
    {
//	      if((*i)[index]==0)vectorsSubset.push_back(*i);
	      if(((*i)&(((int64)1)<<index))==0)vectorsSubset.push_back(*i);
    }

	zeros|=((int64)1)<<index;
	rek64(ones,zeros,nOnes,nZeros+1,vectorsSubset,best,n);
	zeros-=((int64)1)<<index;
}
/*
 * monomialGenerators should have zero-one exponent vectors with at most 64 variables in the ring.
 */
int krullDimensionOfMonomialIdeal64(PolynomialSet const &monomialGenerators)
{
	int n=monomialGenerators.getRing().getNumberOfVariables();
	assert(n<=64);
	vector<ExponentType> generators;
	for(PolynomialSet::const_iterator i=monomialGenerators.begin();i!=monomialGenerators.end();i++)
	{
		ExponentType v=0;
		for(int j=0;j<n;j++)v=2*v+i->getMarked().m.exponent[j];
		generators.push_back(v);
	}

	int best=0;

	ExponentType zeros=0;
	ExponentType ones=0;

	ExponentType possiblyOne=0;
	for(vector<ExponentType>::const_iterator i=generators.begin();i!=generators.end();i++)possiblyOne=possiblyOne|*i;
	ones=allOnes(n)-possiblyOne;

	rek64(ones,zeros,sum(ones),sum(zeros),generators,best,n);

	return best;
}



PolynomialSet radicalOfMonomialIdeal(PolynomialSet const &monomialGenerators)
{
  PolynomialRing theRing=monomialGenerators.getRing();
  PolynomialSet temp=monomialGenerators;

  temp.markAndScale(LexicographicTermOrder()); //just to make sure that some term is marked

  PolynomialSet ret(theRing);
  for(PolynomialSet::const_iterator i=temp.begin();i!=temp.end();i++)
    {
      IntegerVector e=i->getMarked().m.exponent;
      e=e.supportVector();
      ret.push_back(Polynomial(Term(i->getMarked().c,Monomial(theRing,e))));
    }
  return ret;
}

static bool increase(IntegerVector &v, int &numberOfOnes)
{
  int i=0;
  while(i<v.size() && v[i]==1)
    {
      v[i]=0;
      numberOfOnes--;
      i++;
    }
  if(i==v.size())return false;
  v[i]=1;
  numberOfOnes++;
  return true;
}

static int rekCallls;
static void rek(IntegerVector &ones, IntegerVector &zeros, int nOnes, int nZeros, IntegerVectorList const &vectors, int &best)
{
	rekCallls++;
	if(nOnes>best)best=nOnes;
  if(ones.size()-nZeros<best)return;
  if(nOnes+nZeros==ones.size())return;


  log3
    {
      fprintf(Stderr,"Ones:\n");
      AsciiPrinter(Stderr).printVector(ones);
      fprintf(Stderr,"Zeros:\n");
      AsciiPrinter(Stderr).printVector(zeros);

      AsciiPrinter(Stderr).printVectorList(vectors);
    }

  int index=0;
  for(int i=0;i<ones.size();i++,index++)
    if((!ones[i])&&(!zeros[i]))break;

  assert(index<ones.size());

  ones[index]=1;
  bool good=true;
  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)
    if(i->divides(ones))
      {
	good=false;
	break;
      }

  if(good)
    {
      rek(ones,zeros,nOnes+1,nZeros,vectors,best);
    }
  ones[index]=0;

  IntegerVectorList vectorsSubset;

  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)
    {
      if((*i)[index]==0)vectorsSubset.push_back(*i);
    }

  if(0) // It is not clear that this is an improvement, even with an better implementation
  {
	  IntegerVector newOnes(ones.size());
	  IntegerVector possiblyOne(ones.size());
	  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)possiblyOne=max(possiblyOne,*i);

	  if((IntegerVector::allOnes(ones.size())-possiblyOne-ones-zeros).max()>0)
	  {
/*	  	debug<<vectorsSubset;
	  debug<<"ones now\n"<<ones<<"\n";
	  debug<<"zeros now\n"<<zeros<<"\n";
	  debug<<"FORCED ONES:"<<IntegerVector::allOnes(ones.size())-possiblyOne<<"\n";
*/	  newOnes=max(IntegerVector::allOnes(ones.size())-possiblyOne-ones-zeros,IntegerVector(ones.size()));
//	  debug<<"NEW ONES:"<<newOnes<<"\n";
	  }

	  nOnes+=newOnes.sum();
	  ones+=newOnes;

	  zeros[index]=1;
	  rek(ones,zeros,nOnes,nZeros+1,vectorsSubset,best);

	  ones-=newOnes;
	  nOnes-=newOnes.sum();

	  zeros[index]=0;
  }
  else
  {
	  zeros[index]=1;
rek(ones,zeros,nOnes,nZeros+1,vectorsSubset,best);
zeros[index]=0;
  }

  }


int krullDimensionOfMonomialIdeal(PolynomialSet const &monomialGenerators)
{
//	debug<<"Taking radical\n";
	PolynomialSet temp=radicalOfMonomialIdeal(monomialGenerators);
//	debug<<"Minimizing\n";
	minimize(&temp);
//	debug<<"Done\n";
	if(monomialGenerators.getRing().getNumberOfVariables()<=64)return krullDimensionOfMonomialIdeal64(temp);
  IntegerVectorList vectors;
  for(PolynomialSet::const_iterator i=temp.begin();i!=temp.end();i++)
    vectors.push_back(i->getMarked().m.exponent);

  int best=0;
  
  int n=monomialGenerators.getRing().getNumberOfVariables();
  IntegerVector zeros(n);
  IntegerVector ones(n);

  //Preprocessing step. This can be improved.
  IntegerVector possiblyOne(n);
  for(IntegerVectorList::const_iterator i=vectors.begin();i!=vectors.end();i++)possiblyOne=max(possiblyOne,*i);
  ones=IntegerVector::allOnes(n)-possiblyOne;

  rek(ones,zeros,ones.sum(),zeros.sum(),vectors,best);
//debug<<"REKCALLS"<<rekCallls<<"\n";

  return best;
}

int krullDimension(PolynomialSet const &groebnerBasis)
{
  return krullDimensionOfMonomialIdeal(groebnerBasis.markedTermIdeal());
}
