#include "field.h"

#include <memory>
#include <assert.h>
#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>

#include "printer.h"
#include "log.h"
#include "myassert.h"

// The following macro is used for logging of field creation. Note that the usual "log1" does not work since the logLevel might not be initialized at the point of field construction
#define flog1 if(0)
#define flog2 if(0)


FieldElement::FieldElement(Field const &f)
#if USESHAREDPTR
:
  //  theField(f)
implementingObject(f.zHomomorphismImplementation(0))
#endif
{
  assert(f.implementingObject);
#if !USESHAREDPTR
  implementingObject=f.zHomomorphismImplementation(0);//create object with refCount 1
#endif
  //  f.implementingObject->refCount++;
  flog2 fprintf(Stderr,"FieldElement - constructing1\n");
}

FieldElement::FieldElement(class FieldElementImplementation *implementingObject_):
  implementingObject(implementingObject_)
{
}

FieldElement::FieldElement(const FieldElement &a)
    //    theField(a.theField)
{
  //  assert(a.implementingObject);
  if(USESHAREDPTR||a.implementingObject)
    {
      implementingObject=a.implementingObject;
#if !USESHAREDPTR
      implementingObject->refCount++;
#endif
    }
  else
    { //this constructor is called implicitly in splitPoly
//	  assert(0);
	  implementingObject=0;
    }
  flog2 fprintf(Stderr,"FieldElement - constructing2\n");
}

/*int FieldElement::getIntegerRepresentation()const
{
	return implementingObject->getIntegerRepresentation();
}*/


mpq_t const *FieldElement::getGmpRationalTemporaryPointer()const
{
  return (implementingObject)->getGmpRationalTemporaryPointer();
}


double FieldElement::floatingPointApproximation()const
{
	mpq_t const *p=getGmpRationalTemporaryPointer();

	return mpq_get_d(*p);//LEAK??
}


bool FieldElement::isInteger()const
{
  return implementingObject->isInteger();
}


FieldElement FieldElement::multiplierGivingInteger()const
{
  return FieldElement(implementingObject->multiplierGivingInteger());
}


FieldElement::~FieldElement()
{
#if !USESHAREDPTR
  if(implementingObject && 0==(--(implementingObject->refCount)))
    {
      //implementingObject->refCount=0x8000;
      delete implementingObject;
      //      fprintf(Stderr,"DELETE\n");
      flog2 fprintf(Stderr,"Deleting implementing object\n");
    }
#endif

  flog2 fprintf(Stderr,"FieldElement - destructing\n");
}

FieldElement& FieldElement::operator=(const FieldElement& a)
{
#if !USESHAREDPTR
	if(this==&a)
    {
      flog2 fprintf(Stderr,"FieldElement - selfassignment\n");

      return *this;
    }
  if(implementingObject&& 0==(--(implementingObject->refCount)))
    {
      delete implementingObject;
      //      fprintf(Stderr,"DELETE\n");
      flog2 fprintf(Stderr,"FieldElement - deleting old implementing\n");
    }
  //assert(a.implementingObject);
  if(a.implementingObject)
    {
      implementingObject=a.implementingObject;
      implementingObject->refCount++;
    }
  else
    {
      implementingObject=0;
    }
#else
	implementingObject=a.implementingObject;
#endif
  flog2 fprintf(Stderr,"FieldElement - assignment\n");
  return *this;
}

class FieldImplementation *FieldElement::getField()const
{
//	debug<<"getField()\n";
	FieldImplementation *r=implementingObject->getField();
//	debug<<"getField()\n";
  return r;
}

FieldElement FieldElement::one() const//change this implementation
{
  assert(implementingObject);
  return FieldElement(implementingObject->one());
}


bool FieldElement::isZero()const
{
	  assert(implementingObject);
  return implementingObject->isZero();
}

bool FieldElement::isOne()const
{
  if(!implementingObject)return false;
  return (*this-one()).isZero();
}


FieldElement operator+(const FieldElement &a,const FieldElement &b)
{
  assert(a.implementingObject);
  assert(b.implementingObject);
  return FieldElement(a.implementingObject->sum(*(b.implementingObject)));
}


FieldElement operator-(const FieldElement &a,const FieldElement &b)
{
  assert(a.implementingObject);
  assert(b.implementingObject);
  return FieldElement(a.implementingObject->difference(*(b.implementingObject)));
}

FieldElement operator-(const FieldElement &b)
{
  assert(b.implementingObject);
  return FieldElement(b.implementingObject->negation());
}
FieldElement FieldElement::inverse()const
{
  assert(implementingObject);
  return FieldElement(implementingObject->inverse());
}

int FieldElement::sign()const
{
  assert(implementingObject);
  return implementingObject->sign();
}


int FieldElement::pAdicValuation(int p)const
{
  assert(implementingObject);
  return implementingObject->pAdicValuation(p);
}


FieldElement FieldElement::pAdicRemainder(class Field const &ZModPZ)const
{
  assert(implementingObject);
  return FieldElement(implementingObject->pAdicRemainder(ZModPZ));
}


FieldElement FieldElement::modularRepresentative(class Field const &ZModPZ)const
{
	  assert(implementingObject);
	  return FieldElement(implementingObject->modularRepresentative(ZModPZ));
}

int FieldElement::getIntegerRepresentation()const
{
  assert(implementingObject);
  return implementingObject->getIntegerRepresentation();
}

std::string FieldElement::toString(bool writeIfOne, bool alwaysWriteSign, bool latexMode) const
{
  assert(implementingObject);
  return implementingObject->toString(writeIfOne,alwaysWriteSign,latexMode);
}

#if USESHAREDPTR
void FieldElement::makeUnique()
{
//	if(implementingObject->refCount!=1)
		if(implementingObject.use_count()!=1)
	{
		  std::shared_ptr<class FieldElementImplementation> temp(implementingObject->copy());
		  swap(temp,implementingObject);
	}
}
#endif

void FieldElement::operator*=(const FieldElement &a)
{
  assert(a.implementingObject);
  assert(implementingObject);
#if !USESHAREDPTR
  if(implementingObject->refCount!=1)
    {
      implementingObject->refCount--;
      implementingObject=implementingObject->copy();
    }
#else
  makeUnique();
#endif
  (*implementingObject)*=(*a.implementingObject);
}


void FieldElement::operator+=(const FieldElement &a)
{
  assert(a.implementingObject);
  assert(implementingObject);
#if !USESHAREDPTR
  if(implementingObject->refCount!=1)
    {
      implementingObject->refCount--;
      implementingObject=implementingObject->copy();
    }
#else
  makeUnique();
#endif
  (*implementingObject)+=(*a.implementingObject);
}


void FieldElement::madd(const FieldElement &a, const FieldElement &b)
{
  assert(a.implementingObject);
  assert(b.implementingObject);
  assert(implementingObject);
#if !USESHAREDPTR
  if(implementingObject->refCount!=1)
    {
      implementingObject->refCount--;
      implementingObject=implementingObject->copy();
    }
#else
  makeUnique();
#endif
  (*implementingObject).madd(*a.implementingObject,*b.implementingObject);
}


FieldElement operator*(const FieldElement &a,const FieldElement &b)
{
  FieldElement c=a;
  c*=b;
  return c;
}







FieldElement Field::zHomomorphism(int n)const
{
  return implementingObject->zHomomorphism(n);
}


FieldElement Field::random()const
{
  return implementingObject->random();
}


int Field::getCharacteristic()const
{
  return implementingObject->getCharacteristic();
}



FieldElementImplementation::FieldElementImplementation(FieldImplementation &a):
#if !USESHAREDPTR
		refCount(1),
#endif
  theFieldImplementation(a)
{
#if !USESHAREDPTR
	a.refCount++;
#endif
  numberOfLivingFieldElementImplementations++;
};

FieldElementImplementation::~FieldElementImplementation()
{
  numberOfLivingFieldElementImplementations--;
#if !USESHAREDPTR
  if(0==(--(theFieldImplementation.refCount)))
    {
      delete &theFieldImplementation;//NOW THE REFERENCE IS INVALID....
      flog2 fprintf(Stderr,"Deleting field implementation\n");
    }
#endif
};

class FieldImplementation *FieldElementImplementation::getField()const
{
  return &theFieldImplementation;
}

/*class Field FieldElementImplementation::getField2()const
{
	return this->t
}*/

FieldElementImplementation *Field::zHomomorphismImplementation(int n)const
{
  return implementingObject->zHomomorphismImplementation(n);
}

const char *Field::name()
{
  return implementingObject->name();
}


bool Field::isRationals()const
{
  return implementingObject->isRationals();
}

std::string Field::toString()const{
  return implementingObject->toString();
}

Field::Field(FieldImplementation *implObj)
#if USESHAREDPTR
:
		implementingObject(implObj)
#endif
{
#if !USESHAREDPTR
	 implObj->refCount++;
	 implementingObject=implObj;
#endif
  flog1 fprintf(Stderr,"Constructing Field %s\n",toString().c_str());
  flog1 dumpStackTrace();
}

Field::~Field()
{
  assert(implementingObject);
#if !USESHAREDPTR
  implementingObject->refCount--;
  assert(implementingObject->refCount>=0);
  if(implementingObject->refCount==0)
    {
      flog1 fprintf(Stderr,"Deleting implementing object\n");
      delete implementingObject;
    }
  implementingObject=0;
#endif
  flog1   fprintf(Stderr,"Destructing Field\n");
}

Field::Field(Field const &a)//copy constructor
  :implementingObject(a.implementingObject)
{
#if !USESHAREDPTR
  implementingObject->refCount++;
#endif
  flog1 fprintf(Stderr,"Copying field\n");
}


Field& Field::operator=(const Field& a)
{
//	assert(0);
  flog1 fprintf(Stderr,"Assigning Field\n");
#if !USESHAREDPTR
  if(this==&a)
    {
      flog1 fprintf(Stderr,"---selfassigning\n");
      return *this;
    }

  if(implementingObject&& 0==(--(implementingObject->refCount)))
  {
	  debug<<"DELETING IMPLEMENTING OBJECT\n";
	  delete implementingObject;
  }
  //assert(a.implementingObject);
  if(a.implementingObject)
    {
      implementingObject=a.implementingObject;
      implementingObject->refCount++;
    }
  else
    implementingObject=0;
#else
  implementingObject=a.implementingObject;
#endif
  return *this;
}

//Field *Field::list;
//Field *Field::currentField;


//Field::Field()
//{
  /*  if(addToList)
    {
      next=list;
      list=this;
    }
  */
//}



/*Field *Field::find(const char *name)
{
   Field *l=list;
   while(l)
      {
  fprintf(Stderr,"testC\n");

         if(std::string(l->name())==std::string(name))break;
         l=l->next;
      }
  fprintf(Stderr,"testD\n");
   return l;
}
*/
 /*
void Field::printList(FILE *f)
{
   fprintf(f,"List of linked field classes:\n");
   Field *l=list;
   while(l)
      {
         fprintf(f," %s\n",l->name());
         l=l->next;
      }
}
 */

  /*void Field::checkInitialized()
{
  if(!currentField)
    {
  fprintf(Stderr,"test2323\n");
      setField(find("GmpRationals"));
    }
}
  */

   /*void Field::setField(Field *field)
{
  assert(field);
  currentField=field;
  fprintf(Stderr,"Field class being used: \"%s\".\n",currentField->name());
  }*/


    /*FieldElement Field::staticZHomomorphism(int n)
{
  fprintf(Stderr,"test");
  //  checkInitialized();
  return currentField->zHomomorphism(n);
  }*/









int FieldImplementation::numberOfLivingFieldImplementations;
int FieldElementImplementation::numberOfLivingFieldElementImplementations;
