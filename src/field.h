#ifndef FIELD_H_INCLUDED
#define FIELD_H_INCLUDED

#define USESHAREDPTR 0

#include <string.h>
#include <string>
#include <stdio.h>

#include <assert.h>

#if USESHAREDPTR
#include <memory>
#endif

#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>  //SHOULD BE REMOVED

/** A FieldElement is an element of a Field in the algebraic sense. A
    FieldElement always knows its associated Field to which it
    belongs.  We may perform arithmetic operations on a FieldElement
    and a FieldElement is passed around by value. Thus in C++ it makes
    no sense to derive classes from the FieldElement class. The actual
    data is contained in an object of the class
    FieldElementImplementation and the FieldElement serves as a smart
    pointer with reference counting to this implementation.
 */
class FieldElement
{
 public:
#if USESHAREDPTR
  std::shared_ptr<class FieldElementImplementation> implementingObject;
#else
  class FieldElementImplementation *implementingObject;
#endif
  FieldElement(class FieldElementImplementation *implementingObject_);
  void makeUnique();
 public:
  class FieldImplementation *getField()const;
  /* Returns a pointer to the FieldImplementation which does NOT get
     its refCount increased. The pointer should immediately be cast to
     a Field object, there by incresing the refCount. This allows code like this
     void test(FieldElement a){Field b=a.getField();...}.
     Would it also work for this (a+a).getField(); or could the ref count get 0 before it is increased?
  */

/*
 * The following is only supported for Z/pZ
 */
  int getIntegerRepresentation()const;

  /*
 * The following two are only supported for rationals.
 */
  mpq_t const *getGmpRationalTemporaryPointer()const;
  double floatingPointApproximation()const;
  bool isInteger()const;
  virtual FieldElement multiplierGivingInteger()const;

#if !USESHAREDPTR
  void setImplementingObject(FieldElementImplementation *implementingObject_){implementingObject=implementingObject_;assert(0);}
#endif

  FieldElement one() const;
  bool isZero()const;
  virtual bool isOne()const;
  friend FieldElement operator+(const FieldElement &a,const FieldElement &b);
  friend FieldElement operator-(const FieldElement &a,const FieldElement &b);
  friend FieldElement operator-(const FieldElement &b);
  FieldElement inverse()const;
  int sign()const; // Only allowed for ordered field. Asserts otherwise.
  int pAdicValuation(int p)const; // Only allowed for the field Q.
  FieldElement pAdicRemainder(class Field const &ZModPZ)const; // Only allowed for the field Q.
  FieldElement modularRepresentative(class Field const &ZModPZ)const; // Only allowed for the field Q and only for integers in that field.
  int integerRepresentative()const; // Only allowed for Z/pZ
  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false ,bool latexMode=false) const;
  void operator*=(const FieldElement &a);
  void operator+=(const FieldElement &a);
  /** Adds a*b to the value of the object. */
  void madd(const FieldElement &a, const FieldElement &b);
  friend FieldElement operator*(const FieldElement &a,const FieldElement &b);

  // Constructors:
  FieldElement(class Field const &f);//Initializes to zero

    FieldElement()              // This constructor causes a lot of trouble
    {                         // Having a legal FieldElement without an implementing object
#if USESHAREDPTR
    	fprintf(stderr,"TROUBLE\n");
#else
    	      implementingObject=0;   // we must check for each assignment (copy constructor/assignment operator)
#endif
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
#if USESHAREDPTR
  std::shared_ptr<class FieldImplementation> implementingObject;
#else
  class FieldImplementation *implementingObject;
#endif
 public:
  FieldElementImplementation *zHomomorphismImplementation(int n)const;
  /**
     @return The image of the integer n under the unique ring homomorphism from the integers Z to the Field taking 1 to the multiplicative neutral element.
   */
  FieldElement zHomomorphism(int n)const;
  FieldElement random()const;
  int getCharacteristic()const;
  bool isRationals()const;

  const char *name();
  std::string toString()const;

  Field(Field const &a);//copy constructor
  //explicit
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
//  class Field getField2()const;
//  class FieldImplementation *getFieldImplementation()const;
//  class Field getField()const;
  static int getNumberOfLivingFieldElementImplementations(){return numberOfLivingFieldElementImplementations;};
#if !USESHAREDPTR
  int refCount;
#endif
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
  virtual FieldElement modularRepresentative(Field const &ZModPZ)const
  {
	  assert(0);
	  return modularRepresentative(ZModPZ);
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
  virtual FieldElementImplementation *multiplierGivingInteger()const
  {
	  assert(0);
	  return multiplierGivingInteger();
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
#if !USESHAREDPTR
	int refCount;
#endif
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
#if USESHAREDPTR
  FieldImplementation()//:refCount(0)
#else
  FieldImplementation():refCount(0)
#endif
  {
//	  fprintf(stderr,"FieldImplementationCreation\n");
      numberOfLivingFieldImplementations++;
    }
  virtual ~FieldImplementation()
    {
//	  fprintf(stderr,"FieldImplementationDestruction\n");
      numberOfLivingFieldImplementations--;
    }
  virtual bool isRationals()const
  {
    return false;
  }
    //  static Field *find(const char *name);
  //  static void printList(FILE *f);
  virtual FieldElement zHomomorphism(int n)=0;
  virtual FieldElement random()
  {
	  assert(0);
  }
    virtual int getCharacteristic()const=0;
    virtual const char *name()=0;
  virtual std::string toString()const=0;

  //  static void setField(Field *f);                 // Do we really want these
  //  static FieldElement staticZHomomorphism(int n); // two procedures?
};



#endif
