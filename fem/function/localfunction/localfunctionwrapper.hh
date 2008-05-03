#ifndef DUNE_LOCALFUNCTIONWRAPPER_HH
#define DUNE_LOCALFUNCTIONWRAPPER_HH

//-s system includes 
#include <cassert>

//- Dune includes 
#include <dune/fem/storage/objectstack.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  template< class DFTraits >
  class DiscreteFunctionDefault;
//  template< class DiscreteFunctionSpaceType, class LocalFunctionImp >
//  class LocalFunctionDefault;

  template< class LocalFunctionStorage >
  class LocalFunctionWrapper;



  template< class LocalFunctionStorage >
  struct LocalFunctionWrapperTraits
  {
    typedef LocalFunctionStorage LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: ObjectType
      LocalFunctionImpType;
    typedef typename LocalFunctionImpType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef LocalFunctionWrapper< LocalFunctionStorageType >
      LocalFunctionUserType;
  };



  //**************************************************************************
  //
  //  --LocalFunctionWrapper 
  //
  //**************************************************************************
  //! Manages the getting and deleting of local function pointers and 
  //! acts like a local functions 
  template < class LocalFunctionStorage >
  class LocalFunctionWrapper
  : public LocalFunction< LocalFunctionWrapperTraits< LocalFunctionStorage > >
  {
  public:
    //! type of the traits
    typedef LocalFunctionWrapperTraits< LocalFunctionStorage > Traits;

    //! type of the local function storage
    typedef typename Traits :: LocalFunctionStorageType
      LocalFunctionStorageType;

    //! type of wrapped local function implementation 
    typedef typename Traits :: LocalFunctionImpType LocalFunctionImpType;

    //! type of discrete function space 
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
  private:
    typedef LocalFunctionWrapper< LocalFunctionStorageType > ThisType;
    typedef LocalFunction< Traits > BaseType;

    friend class EngineWrapper< LocalFunctionImpType, ThisType >;
    
  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

    enum {dimRange  = DiscreteFunctionSpaceType :: dimRange };
    enum {dimDomain = DiscreteFunctionSpaceType :: dimDomain };

  protected:
    // type of stack entry 
    typedef typename LocalFunctionStorageType :: PointerType
      LocalFunctionImpPtrType;
   
  protected:
    // pair storing pointer to local function and poiner to ref-counter 
    LocalFunctionImpPtrType lfptr_;

    // reference to local function 
    LocalFunctionImpType &lf_;

  public:
    //! Constructor initializing the underlying local function
    inline LocalFunctionWrapper ( const EntityType &entity,
                                  LocalFunctionStorageType &storage )
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
      // init real local function with entity
      asImp().init( entity );
    }

    //! Constructor creating empty local function 
    inline explicit LocalFunctionWrapper ( LocalFunctionStorageType &storage ) 
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
    }

    //! Constructor creating empty local function from given discrete
    //! function 
    template< class DiscreteFunctionType >
    inline explicit LocalFunctionWrapper ( DiscreteFunctionType &discreteFunction ) 
    : lfptr_( discreteFunction.localFunctionStorage().getObject() ),
      lf_( *lfptr_ )
    {
    }

    //! Copy constructor
    inline LocalFunctionWrapper ( const ThisType &other )
    : lfptr_( other.lfptr_ ),
      lf_( *lfptr_ )
    {
    }

  private:
    // prohibit assignment 
    inline ThisType &operator= ( const ThisType & );

  protected:
    const LocalFunctionImpType &asImp () const
    { 
      return lf_;
    } 

    LocalFunctionImpType &asImp () 
    {
      return lf_;
    } 
  }; // end LocalFunctionWrapper  

  

  template< class LocalFunctionFactoryImp >
  class LocalFunctionStack
  : public ObjectStack< LocalFunctionFactoryImp >
  {
  public:
    typedef LocalFunctionFactoryImp LocalFunctionFactoryType;

  private:
    typedef LocalFunctionStack< LocalFunctionFactoryType > ThisType;
    typedef ObjectStack< LocalFunctionFactoryImp > BaseType;

  public:
    typedef LocalFunctionWrapper< ThisType > LocalFunctionType;

  public:
    inline explicit LocalFunctionStack ( const LocalFunctionFactoryType &factory )
    : BaseType( factory )
    {
    }

  private:
    inline LocalFunctionStack ( const ThisType & );

    inline ThisType &operator= ( const ThisType & );

  public:
    LocalFunctionType localFunction ()
    {
      return LocalFunctionType( *this );
    }
    
    template< class EntityType >
    const LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return LocalFunctionType( entity, *this );
    }

    template< class EntityType >
    LocalFunctionType localFunction ( const EntityType &entity )
    {
      return LocalFunctionType( entity, *this );
    }
  };

} // end namespace Dune
#endif
