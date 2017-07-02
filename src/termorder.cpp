#include "termorder.h"

#include <assert.h>

#include "printer.h"

// To do: FIX ROWDOT, MUST USE LONGS

void TermOrder::print(Printer &p)const
{
  p.printString("TermOrder");
  p.printNewLine();
}

void TermOrder::printMatrix(Printer &p, int dim)const
{
  IntegerVectorList l;
  for(int i=0;i<dim;i++)
    {
      IntegerVector row(dim);
      for(int j=0;j<dim;j++)
	{
	  row[j]=rowDot(i,IntegerVector::standardVector(dim,j));
	}
      l.push_back(row);
    }
  p.printVectorList(l);
  p.printNewLine();
}

//-----------------------------------------
// LexicographicTermOrder
//-----------------------------------------

LexicographicTermOrder::LexicographicTermOrder(int largest)
{
  this->largest=largest;
}

int LexicographicTermOrder::rowDot(int row, const IntegerVector &v)const
{
  return v[(unsigned int)(row+largest)%(unsigned int)v.size()];
}

bool LexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  if(a.size()!=b.size())
    {
      fprintf(Stderr,"Lexicographic term order compare failed on the following vectors:\n");
      AsciiPrinter(Stderr).printVector(a);
      fprintf(Stderr,"\n");
      AsciiPrinter(Stderr).printVector(b);
      fprintf(Stderr,"\n");
      assert(a.size()==b.size());
    }
  int n=a.size();
  int nLoop=n;

  if(perturbationDegree>=0)nLoop=perturbationDegree;

  for(int i=0;i<nLoop;i++)
    {
      int index=(unsigned int)(i+largest)%(unsigned int)n;
      if((int64)scaleA*a[index]<(int64)scaleB*b[index])return true;
      if((int64)scaleA*a[index]>(int64)scaleB*b[index])return false;
    }
  return false;
}

/*bool LexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  int n=a.size();
  assert(b.size()==n);

  for(int i=0;i<n;i++)
    {
      if(b[i]>a[i])return true;
      if(b[i]<a[i])return false;
    }
  return false;
}*/


void LexicographicTermOrder::print(Printer &p)const
{
  p.printString("LexicographicTermOrder");
  p.printNewLine();
}


//-----------------------------------------
// LexicographicInvertedTermOrder
//-----------------------------------------

LexicographicInvertedTermOrder::LexicographicInvertedTermOrder()
{
}

int LexicographicInvertedTermOrder::rowDot(int row, const IntegerVector &v)const
{
  return v[v.size()-row-1];
}

bool LexicographicInvertedTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
//  assert(a.size()==b.size());
//  assert(perturbationDegree==-1);
	int n=a.size();

  for(int i=n-1;i>=0;i--)
    {
      if((int64)scaleA*a[i]<(int64)scaleB*b[i])return true;
      if((int64)scaleA*a[i]>(int64)scaleB*b[i])return false;
    }
  return false;
}

/*bool LexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  int n=a.size();
  assert(b.size()==n);

  for(int i=0;i<n;i++)
    {
      if(b[i]>a[i])return true;
      if(b[i]<a[i])return false;
    }
  return false;
}*/


void LexicographicInvertedTermOrder::print(Printer &p)const
{
  p.printString("LexicographicInvertedTermOrder");
  p.printNewLine();
}


//-----------------------------------------
// ReverseLexicographicTermOrder
//-----------------------------------------

ReverseLexicographicTermOrder::ReverseLexicographicTermOrder(int largest)
{
  this->largest=largest;
}


int ReverseLexicographicTermOrder::index(int row, const IntegerVector &a)const
{
    return (unsigned int)(-row+a.size()+largest-1)%(unsigned int)a.size(); //a>b>c //does largest work?
    //return (unsigned int)(row+largest)%(unsigned int)a.size();//a<b<c
}

int ReverseLexicographicTermOrder::rowDot(int row, const IntegerVector &v)const
{
  return -v[index(row,v)];
}


bool ReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  assert(a.size()==b.size());
  int n=a.size();
  int nLoop=n;

  if(perturbationDegree>=0)nLoop=perturbationDegree;

  for(int i=0;i<nLoop;i++)
    {
      int index=this->index(i,a);
      int64 A=(int64)scaleA*a[index];
      int64 B=(int64)scaleB*b[index];
      if(A>B)return true;
      if(A<B)return false;
    }
  return false;
}

/*
bool ReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  int n=a.size();
  assert(b.size()==n);

  for(int i=0;i<n;i++)
    {
      if(b[i]>a[i])return false;
      if(b[i]<a[i])return true;
    }
  return false;
}
*/

void ReverseLexicographicTermOrder::print(Printer &p)const
{
  p.printString("ReverseLexicographicTermOrder");
  p.printNewLine();
}


//-----------------------------------------
// ReverseLexicographicInvertedTermOrder
//-----------------------------------------

int ReverseLexicographicInvertedTermOrder::rowDot(int row, const IntegerVector &v)const
{
  return -v[row];
}


bool ReverseLexicographicInvertedTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  assert(a.size()==b.size());
  int n=a.size();

  for(int i=0;i<n;i++)
    {
      int64 A=(int64)scaleA*a.UNCHECKEDACCESS(i);
      int64 B=(int64)scaleB*b.UNCHECKEDACCESS(i);
      if(A!=B)
      {
    	  if(A>B)return true;
    	  if(A<B)return false;
      }
    }
  return false;
}


void ReverseLexicographicInvertedTermOrder::print(Printer &p)const
{
  p.printString("ReverseLexicographicInvertedTermOrder");
  p.printNewLine();
}


//-----------------------------------------
// StandardGradedLexicographicTermOrder
//-----------------------------------------

StandardGradedLexicographicTermOrder::StandardGradedLexicographicTermOrder(int largest)
{
  this->largest=largest;
}

int StandardGradedLexicographicTermOrder::rowDot(int row, const IntegerVector &v)const
{
  if(row==0)return v.sum();
  row--;

  return v[(unsigned int)(row+largest)%(unsigned int)v.size()];
}


bool StandardGradedLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  assert(a.size()==b.size());
  int n=a.size();

  if(perturbationDegree==0)return false;

  int64 difsum=(int64)scaleA*a.sum()-(int64)scaleB*b.sum();
  if(difsum<0)return true;
  if(difsum>0)return false;

  if(perturbationDegree==1)return false;

  int nLoop=n;
  if(perturbationDegree>=0)nLoop=perturbationDegree-1;

  for(int i=0;i<nLoop;i++)
    {
      int index=(unsigned int)(i+largest)%(unsigned int)n;
      if((int64)scaleA*a[index]<(int64)scaleB*b[index])return true;
      if((int64)scaleA*a[index]>(int64)scaleB*b[index])return false;
    }
  return false;
}

/*bool StandardGradedLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  int n=a.size();
  assert(b.size()==n);

  int difsum=a.sum()-b.sum();
  if(difsum<0)return true;
  if(difsum>0)return false;

  for(int i=0;i<n;i++)
    {
      if(b[i]>a[i])return true;
      if(b[i]<a[i])return false;
    }
  return false;
}*/

void StandardGradedLexicographicTermOrder::print(Printer &p)const
{
  p.printString("StandardGradedLexicographicTermOrder");
  p.printNewLine();
}

//-----------------------------------------
// StandardGradedReverseLexicographicTermOrder
//-----------------------------------------
/* to be fixed some day
StandardGradedReverseLexicographicTermOrder::StandardGradedReverseLexicographicTermOrder(int largest)
{
  this->largest=largest;
}

int StandardGradedLexicographicTermOrder::rowDot(int row, const IntegerVector &v)const
{
  if(row==0)return v.sum();
  row--;

  return v[(unsigned int)(row+largest)%(unsigned int)v.size()];
}


bool StandardGradedLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  assert(a.size()==b.size());
  int n=a.size();

  if(perturbationDegree==0)return false;

  int64 difsum=(int64)scaleA*a.sum()-(int64)scaleB*b.sum();
  if(difsum<0)return true;
  if(difsum>0)return false;

  if(perturbationDegree==1)return false;

  int nLoop=n;
  if(perturbationDegree>=0)nLoop=perturbationDegree-1;

  for(int i=0;i<nLoop;i++)
    {
      int index=(unsigned int)(i+largest)%(unsigned int)n;
      if((int64)scaleA*a[index]<(int64)scaleB*b[index])return true;
      if((int64)scaleA*a[index]>(int64)scaleB*b[index])return false;
    }
  return false;
}

void StandardGradedLexicographicTermOrder::print(Printer &p)const
{
  p.printString("StandardGradedLexicographicTermOrder");
  p.printNewLine();
}
*/
//-----------------------------------------
// WeightTermOrder
//-----------------------------------------

int WeightTermOrder::rowDot(int row, const IntegerVector &v)const
{
  if(row==0)return dot(v,weight);
  row--;

  return LexicographicTermOrder().rowDot(row,v);
}

bool WeightTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
{
  //  fprintf(Stderr,"Perturbation degree: %i\n",perturbationDegree);
  if(perturbationDegree==0)return false;

  int64 d=scaleA*dotLong(a,weight)-scaleB*dotLong(b,weight);
  //  int64 d=dotLong(scaleA*a-scaleB*b,weight);
  /*
  fprintf(Stderr,"sizeof:%i\n",sizeof(int64));
  AsciiPrinter(Stderr).printVector(a);
  AsciiPrinter(Stderr).printVector(b);
  fprintf(Stderr,"%x %x\n",(int)d,(int)(d>>32));
  fprintf(Stderr,"%x\n",(((int64)85569)*((int64)85569))/16);
  */
  if(d<0)return true;
  if(d>0)return false;
  return LexicographicTermOrder()(a,b,scaleA,scaleB,perturbationDegree-1);
}


/*bool WeightTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  int d=dot(a-b,weight);

  if(d<0)return true;
  if(d>0)return false;
  return LexicographicTermOrder()(a,b);
}*/


void WeightTermOrder::print(Printer &p)const
{
  p.printString("WeightTermOrder");
  p.printVector(weight);
  p.printNewLine();
}

//-----------------------------------------
// WeightReverseLexicographicTermOrder
//-----------------------------------------

IntegerVector WeightReverseLexicographicTermOrder::getWeight()const
{
  return weight;
}


int WeightReverseLexicographicTermOrder::rowDot(int row, const IntegerVector &v)const
{
  if(row==0)return dot(v,weight);
  row--;

  return ReverseLexicographicTermOrder().rowDot(row,v);
}

bool WeightReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
  //bool WeightReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
  if(perturbationDegree==0)return false;
  //int d=dot(scaleA*a-scaleB*b,weight);

  int64 d=(int64)scaleA*dotLong(a,weight)-(int64)scaleB*dotLong(b,weight);

  if(d<0)return true;
  if(d>0)return false;
  return ReverseLexicographicTermOrder()(a,b,scaleA,scaleB,perturbationDegree-1);
}


void WeightReverseLexicographicTermOrder::print(Printer &p)const
{
  p.printString("WeightReverseLexicographicTermOrder");
  p<<weight;
  p.printNewLine();
}


//-----------------------------------------
// MatrixTermOrder
//-----------------------------------------

int MatrixTermOrder::rowDot(int row, const IntegerVector &v)const
{
  int nrows=weights.size();//slow
  if(row<nrows)
	  {
		  IntegerVectorList::const_iterator i=weights.begin();
		  for(int j=0;j<nrows;j++)i++;
		  return dot(v,*i);
	  }
  row-=nrows;

  return ReverseLexicographicInvertedTermOrder().rowDot(row,v);
}

bool MatrixTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
  //bool WeightReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
	assert(perturbationDegree);
//  if(perturbationDegree==0)return false;
  //int d=dot(scaleA*a-scaleB*b,weight);

  for(IntegerVectorList::const_iterator i=weights.begin();i!=weights.end();i++)
  {
#if 0
	  int64 d=(int64)scaleA*dotLong(a,*i)-(int64)scaleB*dotLong(b,*i);
#else
	  int64 da=0;
	  int64 db=0;
	  int n=a.size();
	  for(int j=0;j<n;j++)
	  {
		  da+=((int64)a.UNCHECKEDACCESS(j))*i->UNCHECKEDACCESS(j);
		  db+=((int64)b.UNCHECKEDACCESS(j))*i->UNCHECKEDACCESS(j);
	  }
	  int64 d=scaleA*da-scaleB*db;
#endif
	  if(d<0)return true;
	  if(d>0)return false;
  }

  return ReverseLexicographicInvertedTermOrder()(a,b,scaleA,scaleB,perturbationDegree-1);
}


void MatrixTermOrder::print(Printer &p)const
{
  p.printString("MatrixTermOrder");
  p<<weights;
  p.printNewLine();
}

//-----------------------------------------
// TotalDegreeTieBrokenTermOrder
//-----------------------------------------

int TotalDegreeTieBrokenTermOrder::rowDot(int row, const IntegerVector &v)const
{
	if(row==0)return v.sum();
	return tieBreaker.rowDot(row-1,v);
}

bool TotalDegreeTieBrokenTermOrder::operator()(const IntegerVector &a, const IntegerVector &b, int scaleA, int scaleB, int perturbationDegree)const
  //bool WeightReverseLexicographicTermOrder::operator()(const IntegerVector &a, const IntegerVector &b)const
{
//	assert(perturbationDegree);

	int64 d=(int64)scaleA*(int64)a.sum()-(int64)scaleB*(int64)b.sum();

	  if(d<0)return true;
	  if(d>0)return false;

  return tieBreaker(a,b,scaleA,scaleB,perturbationDegree-1);
}


void TotalDegreeTieBrokenTermOrder::print(Printer &p)const
{
  p.printString("TotalDegreeTieBrokenTermOrder");
  tieBreaker.print(p);
  p.printNewLine();
}
