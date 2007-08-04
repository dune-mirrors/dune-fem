#ifndef DUNE_OBJECTSTACK_HH
#define DUNE_OBJECTSTACK_HH

#include <vector>
#include <stack>

namespace Dune{

//! @ingroup HelperClasses
//! Stores pointers to a given class in a stack
//! used for local functions and for basefunctionsets
template <class ObjectFactoryImp> 
class ObjectStack 
{
  typedef ObjectFactoryImp ObjectFactoryType;
private:
  typedef ObjectStack<ObjectFactoryType> MyType;
  typedef typename ObjectFactoryType :: ObjectType ObjectType;

public:  
  typedef typename std::pair<ObjectType* , int* > StackStorageType;
private:
  std::stack < StackStorageType , std::vector<StackStorageType> > objStack_;
  const ObjectFactoryType & factory_;

  int numIssuedObjects_;

  public:
  //! constructor 
  ObjectStack (const ObjectFactoryType & factory) 
    : factory_(factory) , numIssuedObjects_(0) {}

  //! delete all objects on stack 
  ~ObjectStack ()
  {
    assert(numIssuedObjects_ == 0);

    while ( !objStack_.empty() )
    {
      StackStorageType obj = objStack_.top();
      objStack_.pop();
      delete obj.first; 
      obj.first = 0;
      delete obj.second; 
      obj.second = 0;
    }
  }
 
  //! get local function object
  StackStorageType getObject () 
  {
#ifndef NDEBUG
    ++numIssuedObjects_;
#endif
   
    if( objStack_.empty() )
    {
      // first pointer is the local function pointer 
      // ans second pointer is the reference counter initialized with 1  
      return StackStorageType ( factory_.newObject() , new int (1) );
    }
    else 
    {
      StackStorageType obj = objStack_.top();
      objStack_.pop();
      return obj;
    }
  }

  //! push local function to stack 
  void freeObject (StackStorageType & obj)
  {
#ifndef NDEBUG
    --numIssuedObjects_;
#endif
    objStack_.push(obj);
  }
  
private:
  //! prohibited methods 
  ObjectStack ( const MyType & c);
  MyType & operator = ( const MyType & c );
};

} // end namespace Dune
#endif
