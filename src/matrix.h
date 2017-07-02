#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>
#include <list>
#include <assert.h>
#include <algorithm>
#include "vektor.h"
#include "printer.h"

using namespace std;

template <class typ> class Matrix{
  //public:
  int width,height;
//  vector<Vektor<typ> > rows;
  vector<typ> data;
public:
  inline int getHeight()const{return height;};
  inline int getWidth()const{return width;};
  Matrix(const Matrix &a):data(a.data),width(a.getWidth()),height(a.getHeight()){
  }
  Matrix(int height_, int width_):data(height_*width_),height(height_),width(width_){
    assert(height>=0);
    assert(width>=0);
    for(int i=0;i<height_*width_;i++)data[i]=0;
  };
  ~Matrix(){
  };
  Matrix():width(0),height(0){
  };
  Vektor<typ> column(int i)const
    {
      assert(i>=0);
      assert(i<getWidth());
      Vektor<typ> ret(getHeight());
      for(int j=0;j<getHeight();j++)ret[j]=data[j*getWidth()+i];
      return ret;
    }
  Matrix transposed()const
    {
      Matrix ret(getWidth(),getHeight());
      for(int i=0;i<getWidth();i++)
    	  for(int j=0;j<getHeight();j++)
    		  ret.data[i*getHeight()+j]=data[j*getWidth()+i];
      return ret;
    }
  static Matrix identity(int n)
    {
      Matrix m(n,n);
      for(int i=0;i<n;i++)m.data[i*(n+1)]=typ(1);
      return m;
    }
  void append(Matrix const &m)
    {
	  assert(getWidth()==m.getWidth());
	  assert(&m!=this);//If the two matrices are the same object the following will not work:
	  data.insert(data.end(),m.data.begin(),m.data.end());
      height+=m.height;
    }
  void appendRow(Vektor<typ> const &r)
    {
      assert(r.size()==width);
      int oldsize=data.size();
      data.resize(data.size()+width);
      for(int i=0;i<width;i++)data[oldsize+i]=r.v[i];
      height++;
    }
/*  void prependRow(Vektor<typ> const &r)
    {
      assert(r.size()==width);
      rows.front_back(r);
      height++;
    }*/
  void setRow(int i,IntegerVector const &v)
  {
	  assert(v.size()==getWidth());
	  assert(i>=0 && i<getHeight());
	  for(int j=0;j<getWidth();j++)data[i*getWidth()+j]=v.v[j];
  }
  IntegerVector getRow(int i)const
  {
	  IntegerVector ret(getWidth());
	  for(int j=0;j<getWidth();j++)ret.v[j]=data[i*getWidth()+j];
	  return ret;
  }
  IntegerVectorList getRows()const
    {
      IntegerVectorList ret;
      for(int i=0;i<height;i++)ret.push_back(getRow(i));
      return ret;
    }
  typ dotRow(IntegerVector const &v, int i)const
  {
	  assert(v.size()==getWidth());
	  typ s=0;
	  for(int j=0;j<getWidth();j++)
		  s+=data[j+i*getWidth()]*v.v[j];
	  return s;
  }
  typ dotRowLong(IntegerVector const &v, int i)const
  {
	  assert(v.size()==getWidth());
	  int64 s=0;
	  for(int j=0;j<getWidth();j++)
		  s+=((int64)(data[j+i*getWidth()]))*v.v[j];
	  return s;
  }
  IntegerVector vectormultiply(IntegerVector const &v)const
    {
      assert(v.size()==width);
      IntegerVector ret(height);
      for(int i=0;i<height;i++)
    	  ret[i]=dotRow(v,i);
      return ret;
    }
  /**
   * Decides if v is in the kernel of the matrix.
   */
  bool inKernel(IntegerVector const &v)const
    {
      assert(v.size()==width);
      for(int i=0;i<height;i++)
    	  if(dotRowLong(v,i)!=0)return false;
      return true;
    }
  inline friend Matrix operator*(typ s, const Matrix& q)
    {
      Matrix p=q;
      for(int i=0;i<q.data.size();i++)p.data[i]=s*q.data[i];
      return p;
    }
  friend Matrix tropicalProduct(const Matrix& a, const Matrix& b)
  {
    int neutral=-888888888;
    assert(a.width==b.height);
    Matrix ret(a.height,b.width);
    for(int i=0;i<a.height;i++)
      for(int j=0;j<b.width;j++)
        {
          int sum=neutral;
          for(int k=0;k<b.height;k++)
            {
              if(a[i][k]+b[k][j]>sum)sum=a[i][k]+b[k][j];
            }
          ret[i][j]=sum;
        }
    return ret;
  }
  friend Matrix operator*(const Matrix& a, const Matrix& b)
    {
      assert(a.width==b.height);
      Matrix ret(b.width,a.height);
      for(int i=0;i<b.width;i++)
        ret.setRow(i,a.vectormultiply(b.column(i)));
      return ret.transposed();
    }
  Matrix operator-()const
    {
      Matrix ret(height,width);
      for(int i=0;i<height*width;i++)
          ret.data[i]=-((*this).data[i]);
      return ret;
    }
  /*  template<class T>
    Matrix<T>(const Matrix<T>& c):v(c.size()){
    for(int i=0;i<size();i++)v[i]=typ(c[i]);}
  */

  /**
     Returns the specified submatrix. The endRow and endColumn are not included.
   */
  Matrix submatrix(int startRow, int startColumn, int endRow, int endColumn)const
  {
    assert(startRow>=0);
    assert(startColumn>=0);
    assert(endRow>=startRow);
    assert(endColumn>=startColumn);
    assert(endRow<=height);
    assert(endColumn<=width);
    Matrix ret(endRow-startRow,endColumn-startColumn);
    for(int i=startRow;i<endRow;i++)
      for(int j=startColumn;j<endColumn;j++)
//          ret[i-startRow][j-startColumn]=rows[i][j];
    	  ret.data[(i-startRow)*ret.width+j-startColumn]=data[i*width+j];
    return ret;
  }
  Matrix submatrixRows(int startRow, int endRow)const
  {
	  return submatrix(startRow,0,endRow,getWidth());
  }
  /**
     Returns the specified submatrix. The height is the height of *this, but the returned matrix only contains
     those columns i, for which subset[i]==1. All other entries of subset must be zero.
   */
  Matrix submatrixColumnSubsetBoolean(IntegerVector const &subset)const
  {
    assert(subset.size()==this->getWidth());
    Matrix ret(getHeight(),subset.sum());
    for(int i=0;i<getHeight();i++)
      {
        int J=0;
        for(int j=0;j<subset.size();j++)
          if(subset[j])ret[i][J++]=data[i*getWidth()+j];
      }
    return ret;
  }

  class RowRef;
  class const_RowRef{
    int rowNumM;
    Matrix const &matrix;
    friend class RowRef;
  public:
  inline const_RowRef(const Matrix  &matrix_, int rowNum_):
    rowNumM(rowNum_*matrix_.width),
      matrix(matrix_)
      {
      }
  inline typ const &operator[](int j)const
    {
	assert(j>=0);
	assert(j<matrix.width);
	return matrix.data[rowNumM+j];
    }
  inline typ const &UNCHECKEDACCESS(int j)const
    {
	return matrix.data[rowNumM+j];
    }
    const Vektor<typ> toVector()const
    {
      Vektor<typ> ret(matrix.width);
      for(int j=0;j<matrix.width;j++)
    	  ret[j]=matrix.data[rowNumM+j];
      return ret;
    }
//    operator const IntegerVector()const
    operator Vektor<typ>()const
		{
			return toVector();
		}
    bool operator==(Vektor<typ> const &b)const
		{
			return toVector()==b;
		}
    typ dot(Vektor<typ> const &b)const
		{
			return dot(toVector(),b);
		}
    Vektor<typ> operator-()const
    {
    	return -toVector();
    }
/*    bool isZero()const
    {
      for(int j=0;j<matrix.width;j++)if(!matrix.isZero(matrix.data[matrix.width*rowNum+j]))return false;
      return true;
    }*/
  };

  class RowRef{
    int rowNumM;
    Matrix &matrix;
  public:
  inline RowRef(Matrix &matrix_, int rowNum_):
    rowNumM(rowNum_*matrix_.width),
      matrix(matrix_)
      {
      }
    inline typ &operator[](int j)
      {
    	assert(j>=0);
    	assert(j<matrix.width);
    	return matrix.data[rowNumM+j];
      }
    inline typ &UNCHECKEDACCESS(int j)
      {
    	return matrix.data[rowNumM+j];
      }
    RowRef &operator=(Vektor<typ> const &v)
    {
        assert(v.size()==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v[j];

    	return *this;
    }
    RowRef &operator=(RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    RowRef &operator+=(Vektor<typ> const &v)
    {
        assert(v.size()==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.v[j];

    	return *this;
    }
    RowRef &operator+=(RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    RowRef &operator+=(const_RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    RowRef &operator=(const_RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    const Vektor<typ> toVector()const
    {
      Vektor<typ> ret(matrix.width);
      for(int j=0;j<matrix.width;j++)
    	  ret[j]=matrix.data[rowNumM+j];
      return ret;
    }
//    operator const IntegerVector()const
    operator Vektor<typ>()const
		{
			return toVector();
		}
    typ dot(Vektor<typ> const &b)const
		{
			return dot(toVector(),b);
		}
    /*    Vector toVector()
    {
      Vector ret(matrix.width);
      for(int j=0;j<matrix.width;j++)
	ret[j]=matrix.data[matrix.width*rowNum+j];
      return ret;
    }
    void set(Vector const &v)
    {
      assert(v.size()==matrix.width);
      for(int j=0;j<matrix.width;j++)
	matrix.data[matrix.width*rowNum+j]=v[j];
    }
    bool isZero()const
    {
      for(int j=0;j<matrix.width;j++)if(!matrix.isZero(matrix.data[matrix.width*rowNum+j]))return false;
      return true;
    }*/
  };

//  const Vektor<typ>& operator[](int n)const{assert(n>=0 && n<getHeight());return (rows[n]);}
//  Vektor<typ>& operator[](int n){assert(n>=0 && n<getHeight());return (rows[n]);}
// Bugfix for gcc4.5 (passing assertion to the above operator:
//  Vektor<typ>& operator[](int n){if(!(n>=0 && n<getHeight()))(*(const Matrix<typ>*)(this))[n];return (rows[n]);}

  inline RowRef operator[](int i)
  {
    assert(i>=0);
    assert(i<height);
    return RowRef(*this,i);
  }
  inline const_RowRef operator[](int i)const
  {
    assert(i>=0);
    assert(i<height);
    return const_RowRef(*this,i);
  }

  /**
     Takes two matrices with the same number of columns and construct
     a new matrix which has the rows of the matrix top on the top and
     the rows of the matrix bottom at the bottom. The return value is
     the constructed matrix.
   */
  friend Matrix combineOnTop(Matrix const &top, Matrix const &bottom)
  {
    assert(bottom.getWidth()==top.getWidth());
    Matrix ret(top.getHeight()+bottom.getHeight(),top.getWidth());
    for(int i=0;i<top.getHeight()*top.getWidth();i++)ret.data[i]=top.data[i];
    for(int i=0;i<bottom.getHeight()*top.getWidth();i++)ret.data[i+top.getHeight()*top.getWidth()]=bottom.data[i];

    return ret;
  }
  /**
     Takes two matrices with the same number of rows and construct
     a new matrix which has the columns of the matrix left on the left and
     the columns of the matrix right on the right. The return value is
     the constructed matrix.
   */
  friend Matrix combineLeftRight(Matrix const &left, Matrix const &right)
  {
    assert(left.getHeight()==right.getHeight());
    Matrix ret(left.getHeight(),left.getWidth()+right.getWidth());
    for(int i=0;i<left.getHeight();i++)
      {
        for(int j=0;j<left.getWidth();j++)ret.data[i*ret.getWidth()+j]=left.data[i*left.getWidth()+j];
        for(int j=0;j<right.getWidth();j++)ret.data[i*ret.getWidth()+j+left.getWidth()]=right.data[i*right.getWidth()+j];
      }
    return ret;
  }
  friend Printer& operator<<(Printer &p, Matrix const &m)
  {
	  Vektor<int> widths(m.getWidth());
	  for(int i=0;i<m.getHeight();i++)
		  for(int j=0;j<m.getWidth();j++)
		  {
			  stringstream s;
			  s<<m[i][j];
			  if(s.str().length()>widths[j])widths[j]=s.str().length();
		  }
	  stringstream s;
	  s<<"{";
	  for(int i=0;i<m.getHeight();i++)
	  {
		  if(i)s<<",";
		  s<<"\n";
		  s<<"(";
		  for(int j=0;j<m.getWidth();j++)
		  {
			  if(j)s<<",";
			  stringstream s2;
			  s2<<m[i][j];
			  for(int k=s2.str().length();k<widths[j];k++)s<<" ";
			  s<<s2.str();
		  }
		  s<<")";
	  }
	  s<<"}\n";

	  p<<s.str();
	  return p;
  }
  const typ& UNCHECKEDACCESS(int i,int j)const{
/*	    assert(i>=0);
	    assert(i<height);
	    assert(j>=0);
	    assert(j<width);*/
	  return (data[i*getWidth()+j]);}
  typ& UNCHECKEDACCESS(int i,int j){
/*	    assert(i>=0);
	    assert(i<height);
	    assert(j>=0);
	    assert(j<width);*/
	    return (data[i*getWidth()+j]);}
  typ* __restrict__ UNCHECKEDACCESSRESTRICT(int i,int j){return & (data[i*getWidth()+j]);}
  void swapColumns(int a, int b)
  {
	  if(a!=b)
		  for(int i=0;i<getHeight();i++)
			  swap(data[i*getWidth()+a],data[i*getWidth()+b]);
  }
};

typedef Matrix<int> IntegerMatrix;
typedef Matrix<double> FloatMatrix;

IntegerMatrix rowsToIntegerMatrix(IntegerVectorList const &rows, int width=-1);//width specifies the matrix width. If no width is specied the width is found by looking at the length of the rows. The function "asserts" if the length of the rows does not match the matrix size or if the width was not specified and could not be read off from the rows.
IntegerMatrix rowToIntegerMatrix(IntegerVector const &row);

FloatMatrix integerToFloatMatrix(IntegerMatrix const &m);
IntegerVector flattenMatrix(IntegerMatrix const &m);
int rank(IntegerMatrix const &m);

#endif
