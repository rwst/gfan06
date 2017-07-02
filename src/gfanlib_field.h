/*
 * gfanlib_field.h
 *
 *  Created on: Sep 2, 2015
 *      Author: anders
 */

#include <string>

#include <assert.h>
#include <memory>

#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>  //SHOULD BE REMOVED

#ifndef GFANLIB_FIELD_H_
#define GFANLIB_FIELD_H_

#define GFANLIB_FIELDELEMENTREFCOUNT 1
namespace gfan{


/** A FieldElement is an element of a Field in the algebraic sense. A
    FieldElement always knows its associated Field to which it
    belongs.  We may perform arithmetic operations on a FieldElement
    and a FieldElement is passed around by value. Thus in C++ it makes
    no sense to derive classes from the FieldElement class. The actual
    data is contained in an object of the class
    FieldElementImplementation and the FieldElement serves as a smart
    pointer without reference counting to avoid multi-threaded synchronisation.
 */
class FieldElement
{
// public:
#if GFANLIB_FIELDELEMENTREFCOUNT
	  std::shared_ptr<class FieldElementImplementation> implementingObject;
#else
  class FieldElementImplementation *implementingObject;
#endif
  void makeUnique()
  {
#if GFANLIB_FIELDELEMENTREFCOUNT
	  if(implementingObject.use_count()!=1)
	{
		  std::shared_ptr<class FieldElementImplementation> temp(implementingObject->copy());
		  swap(temp,implementingObject);
	}
#endif
  }
  FieldElement(class FieldElementImplementation *implementingObject_);
 public:
  class FieldImplementation *getField()const;
  /* Returns a pointer to the FieldImplementation which does NOT get
     its refCount increased. The pointer should immediately be cast to
     a Field object, there by increasing the refCount. This allows code like this
     void test(FieldElement a){Field b=a.getField();...}.
     Would it also work for this (a+a).getField(); or could the ref count get 0 before it is increased?
  */
#if 0
/*
 * The following is only supported for Z/pZ
 */
  int getIntegerRepresentation()const;
  int integerRepresentative()const; // Only allowed for Z/pZ

  /*
 * The following two are only supported for rationals.
 */
  mpq_t const *getGmpRationalTemporaryPointer()const;
  double floatingPointApproximation()const;
  bool isInteger()const;
  int pAdicValuation(int p)const; // Only allowed for the field Q.
  FieldElement pAdicRemainder(class Field const &ZModPZ)const; // Only allowed for the field Q.
#endif
  FieldElement one() const
  {
	  assert(implementingObject);
	  return FieldElement(implementingObject->one());
  }
  bool isZero()const
  {
	  assert(implementingObject);
	  return implementingObject->isZero();
  }

  bool isOne()const
  {
	  if(!implementingObject)return false;
	  return (*this-one()).isZero();
  }

  friend FieldElement operator+(const FieldElement &a,const FieldElement &b)
  {
	  assert(a.implementingObject);
	  assert(b.implementingObject);
	  return FieldElement(a.implementingObject->sum(*(b.implementingObject)));
  }

  friend FieldElement operator-(const FieldElement &a,const FieldElement &b)
  {
	  assert(a.implementingObject);
	  assert(b.implementingObject);
	  return FieldElement(a.implementingObject->difference(*(b.implementingObject)));
  }

  friend FieldElement operator-(const FieldElement &b)
  {
	  assert(b.implementingObject);
	  return FieldElement(b.implementingObject->negation());
  }

  FieldElement inverse()const
  {
	  assert(implementingObject);
	  return FieldElement(implementingObject->inverse());
  }

  int sign()const // Only allowed for ordered field. Asserts otherwise.
  {
	  assert(implementingObject);
	  return implementingObject->sign();
  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false ,bool latexMode=false) const;
  void operator*=(const FieldElement &a)
		{
	  assert(a.implementingObject);
	  assert(implementingObject);
	  makeUnique();
	  (*implementingObject)*=(*a.implementingObject);
		}
  void operator+=(const FieldElement &a)
		{
	  assert(a.implementingObject);
	  assert(implementingObject);
	  makeUnique();
	  (*implementingObject)+=(*a.implementingObject);
		}
  /** Adds a*b to the value of the object. */
  void madd(const FieldElement &a, const FieldElement &b)
  {
	  assert(a.implementingObject);
	  assert(b.implementingObject);
	  assert(implementingObject);
	  makeUnique();
	  (*implementingObject).madd(*a.implementingObject,*b.implementingObject);
  }
  friend FieldElement operator*(const FieldElement &a,const FieldElement &b)
  {
	  FieldElement c=a;
	  c*=b;
	  return c;
  }
  // Constructors:
  FieldElement(class Field const &f);//Initializes to zero

    FieldElement()              // This constructor causes a lot of trouble
    {                         // Having a legal FieldElement without an implementing object
//      implementingObject=0;   // we must check for each assignment (copy constructor/assignment operator)
    }                         // if the pointer is zero. - And some operation will be illegal on
                              // on this element ( one(),..... )

  FieldElement(const FieldElement &a);
  FieldElement& operator=(const FieldElement& a);
  virtual ~FieldElement();
};

/** The Field class describes an object which has as its value a field
    (in the algebraic sense). The value/object can be copied around as
    a value. The Field object serves as a smart pointer with reference
    counting to a FieldImplementation which is hidden for the user of
    the Field class. In C++ it makes no sense to derive classes from
    this class because the object must be passable by value. */

class Field
{
  friend class FieldElement;
  // protected:
 public:
//  class FieldImplementation *implementingObject;
  std::shared_ptr<class FieldImplementation> implementingObject;
 public:
  FieldElementImplementation *zHomomorphismImplementation(int n)const;
  /**
     @return The image of the integer n under the unique ring homomorphism from the integers Z to the Field taking 1 to the multiplicative neutral element.
   */
  FieldElement zHomomorphism(int n)const;
  int getCharacteristic()const;
  bool isRationals()const;

  const char *name();
  std::string toString()const;

  Field(Field const &a);//copy constructor
  Field(FieldImplementation *implObj);//constructor
  Field& operator=(const Field& a);//assignment
  ~Field();//destructor
};


class FieldElementImplementation
{
  static int numberOfLivingFieldElementImplementations;

  class FieldImplementation &theFieldImplementation;
  /* Ideally FieldElement would contain a Field object. However, since methods of a Field should be able to return FieldElements this seems impossible because of the stupid one-pass convention of C++. Instead FieldElement contains a FieldImplementation pointer and FieldElement must thus do the reference counting for this pointer on its own ???*/
 public:
  class FieldImplementation *getField()const;
  static int getNumberOfLivingFieldElementImplementations(){return numberOfLivingFieldElementImplementations;};
//  int refCount;
  FieldElementImplementation(FieldImplementation &a);//ref count = 1??
    virtual ~FieldElementImplementation();

  virtual FieldElementImplementation *one() const=0;

  virtual FieldElementImplementation *copy()const=0;

  virtual  bool isZero()const=0;
  virtual FieldElementImplementation *sum(const FieldElementImplementation &b)const=0;
  virtual FieldElementImplementation *difference(const FieldElementImplementation &b)const=0;
  virtual FieldElementImplementation *negation()const=0;
  virtual FieldElementImplementation *inverse()const=0;
  virtual std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false) const=0;
  virtual void operator*=(const FieldElementImplementation &a)=0;
  virtual void operator+=(const FieldElementImplementation &a)=0;
  virtual void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)=0;
  virtual int sign()const
  {
    assert(0);
    return 0;
  }
  virtual int pAdicValuation(int p)const
  {
    assert(0);
    return 0;
  }
  virtual FieldElement pAdicRemainder(Field const &ZModPZ)const
  {
    assert(0);
    return pAdicRemainder(ZModPZ);
  }
/*  virtual int integerRepresentative()const
  {
    assert(0);
    return 0;
  }*/
  virtual mpq_t const *getGmpRationalTemporaryPointer()const
    {
      fprintf(stderr,"*this object is not implemented using GMP\n");
      assert(0);

      return 0;
    }
  virtual int getIntegerRepresentation()const
  {
	  fprintf(stderr,"*this is not in Z/pZ.\n");
	  assert(0);
  }
  virtual bool isInteger()const
  {
    assert(0);
    return false;
  }

  Field& operator=(const Field& a)
    {
      assert(0);
    }//assignment
};


extern int FieldElementRationalsLiving;



class FieldImplementation
{
 public:
//  int refCount;
  //  static class Field *currentField;               // ??
  //  static void checkInitialized();                 // ??
 protected:
  //  class Field *next;
  //  static class Field *list;
  static int numberOfLivingFieldImplementations;
 public:
  virtual FieldElementImplementation *zHomomorphismImplementation(int n)=0;/* Creates FieldElementImplementation object with refcount1 */
 public:
  static int getNumberOfLivingFieldImplementations(){return numberOfLivingFieldImplementations;};
  FieldImplementation()//:refCount(0)
    {
      numberOfLivingFieldImplementations++;
    }
  virtual ~FieldImplementation()
    {
      numberOfLivingFieldImplementations--;
    }
  virtual bool isRationals()const
  {
    return false;
  }
    //  static Field *find(const char *name);
  //  static void printList(FILE *f);
    virtual FieldElement zHomomorphism(int n)=0;
    virtual int getCharacteristic()const=0;
    virtual const char *name()=0;
  virtual std::string toString()const=0;

  //  static void setField(Field *f);                 // Do we really want these
  //  static FieldElement staticZHomomorphism(int n); // two procedures?
};


}

#endif /* GFANLIB_FIELD_H_ */
