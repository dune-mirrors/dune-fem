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


  //**************************************************************************
  //
  //  --LocalFunctionWrapper 
  //
  //**************************************************************************
  //! Manages the getting and deleting of local function pointers and 
  //! acts like a local functions 
  template < class LocalFunctionStorage >
  class LocalFunctionWrapper
  : public LocalFunction
    < typename LocalFunctionStorage :: ObjectType :: DiscreteFunctionSpaceType,
      typename LocalFunctionStorage :: ObjectType,
      LocalFunctionWrapper< LocalFunctionStorage >
    >
  {
  public:
    //! type of the local function storage
    typedef LocalFunctionStorage LocalFunctionStorageType;

    //! type of wrapped local function implementation 
    typedef typename LocalFunctionStorageType :: ObjectType LocalFunctionImpType;

    //! type of discrete function space 
    typedef typename LocalFunctionImpType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
  private:
    typedef LocalFunctionWrapper< LocalFunctionStorageType > ThisType;
    typedef LocalFunction< DiscreteFunctionSpaceType, LocalFunctionImpType, ThisType >
      BaseType;

    friend class LocalFunction
      < DiscreteFunctionSpaceType, LocalFunctionImpType, ThisType >;
    
  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

#if 0
    enum { dimrange = DiscreteFunctionSpaceType::DimRange };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { DimRange = DiscreteFunctionSpaceType::DimRange };
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { DimDomain = DiscreteFunctionSpaceType::DimDomain };
#endif

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
