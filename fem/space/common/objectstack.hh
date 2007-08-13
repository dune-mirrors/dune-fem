#ifndef DUNE_OBJECTSTACK_HH
#define DUNE_OBJECTSTACK_HH

#include <vector>
#include <stack>

#include <dune/fem/misc/debug.hh>

namespace Dune
{

  //! \ingroup HelperClasses
  //! Stores pointers to a given class in a stack
  //! used for local functions and for basefunctionsets
  template< class ObjectFactoryImp >
  class ObjectStack 
  {
    typedef ObjectFactoryImp ObjectFactoryType;

  private:
    typedef ObjectStack< ObjectFactoryType > ThisType;

  public:
    //! type of the stored objects
    typedef typename ObjectFactoryType :: ObjectType ObjectType;

    //! type of the storage objects (includes an additional reference counter)
    typedef typename std :: pair< ObjectType* , int* > StackStorageType;
    
  protected:
    std :: stack < StackStorageType, std::vector< StackStorageType > > objStack_;
    const ObjectFactoryType &factory_;

    DebugCounter<> numIssuedObjects_;

  public:
    //! constructor 
    ObjectStack ( const ObjectFactoryType &factory )
    : factory_( factory )
    {
    }

  private:
    // Disallow copying
    ObjectStack ( const ThisType &other );

  public:
    //! delete all objects on stack 
    ~ObjectStack ()
    {
      assert( numIssuedObjects_ == 0 );

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
    
  private:
    // Disallow copying
    ThisType &operator= ( const ThisType &other );

  public:
    //! get local function object
    StackStorageType getObject ()
    {
      ++numIssuedObjects_;
      if( objStack_.empty() )
      {
        // first pointer is the local function pointer 
        // ans second pointer is the reference counter initialized with 1  
        return StackStorageType( factory_.newObject() , new int( 1 ) );
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
      --numIssuedObjects_;
      objStack_.push(obj);
    }
  };

} // end namespace Dune

#endif
