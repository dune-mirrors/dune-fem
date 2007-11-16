#ifndef DUNE_LOCALFUNCTIONWRAPPER_HH
#define DUNE_LOCALFUNCTIONWRAPPER_HH

//-s system includes 
#include <cassert>

//- Dune includes 
#include <dune/fem/storage/objectstack.hh>

namespace Dune
{

  template< class DFTraits >
  class DiscreteFunctionDefault;
  template< class DiscreteFunctionSpaceType, class LocalFunctionImp >
  class LocalFunctionDefault;


  //**************************************************************************
  //
  //  --LocalFunctionWrapper 
  //
  //**************************************************************************
  //! Manages the getting and deleting of local function pointers and 
  //! acts like a local functions 
  template < class LocalFunctionStorageImp >
  class LocalFunctionWrapper
  : public LocalFunctionDefault
    < typename LocalFunctionStorageImp :: ObjectType :: DiscreteFunctionSpaceType,
      LocalFunctionWrapper< LocalFunctionStorageImp >
    >
  {
  public:
    //! type of the local function storage
    typedef LocalFunctionStorageImp LocalFunctionStorageType;

    //! type of wrapped local function implementation 
    typedef typename LocalFunctionStorageType :: ObjectType WrappedLocalFunctionType;

    //! type of discrete function space 
    typedef typename WrappedLocalFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    
  private:
    typedef LocalFunctionWrapper< LocalFunctionStorageType > ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;

  public:
    //! type of base function set 
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    enum { dimrange = DiscreteFunctionSpaceType::DimRange };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { DimRange = DiscreteFunctionSpaceType::DimRange };
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { DimDomain = DiscreteFunctionSpaceType::DimDomain };

    // no docu here, then docu is copied from base class 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

  protected:
    // type of stack entry 
    typedef typename LocalFunctionStorageType :: PointerType
      WrappedLocalFunctionPtrType;
   
  private:
    // pair storing pointer to local function and poiner to ref-counter 
    WrappedLocalFunctionPtrType lfptr_;

    // reference to local function 
    WrappedLocalFunctionType &lf_;

  public:
    using BaseType :: axpy;
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! Constructor initializing the underlying local function 
    template< class EntityType > 
    inline LocalFunctionWrapper ( const EntityType &entity,
                                  LocalFunctionStorageType &storage )
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
      // init real local function with entity
      localFunction().init( entity );
    }

    //! Constructor creating empty local function 
    inline explicit LocalFunctionWrapper ( LocalFunctionStorageType &storage ) 
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
    }

    //! Constructor creating empty local function from given discrete
    //! function 
    template <class DiscreteFunctionImp> 
    inline explicit LocalFunctionWrapper ( DiscreteFunctionImp& discreteFunction ) 
    : lfptr_( discreteFunction.localFunctionStorage().getObject() ),
      lf_( *lfptr_ )
    {
    }

    //! Copy constructor
    inline LocalFunctionWrapper ( const LocalFunctionWrapper &org )
    : lfptr_( org.lfptr_ ),
      lf_( *lfptr_ )
    {
    }

    //! destructor pushing local funciton back to the stack (just forget the pointer)
    ~LocalFunctionWrapper () 
    { 
    }
    
  private:
    // prohibit assignment 
    inline ThisType &operator= ( const ThisType& );

  public:
    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) const */
    const RangeFieldType &operator[] ( const int num ) const
    {
      return localFunction()[ num ];
    }
    
    /** \copydoc Dune::LocalFunctionInterface::operator[](const int num) */
    RangeFieldType &operator[] ( const int num )
    {
      return localFunction()[ num ];
    }
 
    /** \copydoc LocalFunctionInterface::numDofs */
    int numDofs () const
    {
      return localFunction().numDofs();
    }
    
    /** \copydoc Dune::LocalFunctionInterface::evaluate(const PointType &x,RangeType &ret) const */
    template< class PointType >
    void evaluate ( const PointType &x,
                    RangeType &ret ) const
    {
      localFunction().evaluate( x , ret );
    }
    
#if 0
    /** \copydoc Dune::LocalFunctionInterface::evaluate( const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    void evaluate ( const QuadratureType &quadrature,
                    const int quadPoint,
                    RangeType &ret ) const
    {
      localFunction().evaluate( quadrature, quadPoint, ret );
    }
#endif
    
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const PointType &x,JacobianRangeType &ret) const */
    template< class PointType >
    void jacobian ( const PointType& x,
                    JacobianRangeType &ret ) const
    {
      localFunction().jacobian( x, ret ); 
    }
  
#if 0
    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType > 
    void jacobian ( const QuadratureType &quadrature,
                    const int quadPoint,
                    JacobianRangeType &ret ) const
    {
      localFunction().jacobian( quadrature, quadPoint, ret );
    }
#endif
   
    /** \brief update local function for given Entity
     */
    template <class EntityType > 
    void init ( const EntityType &en )
    { 
      localFunction().init(en);
    } 
    
    /** \copydoc Dune::LocalFunctionInterface::axpy(const PointType &x,const RangeType &factor) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor )
    {
      localFunction().axpy( x, factor );
    }

#if 0
    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const RangeType &factor) */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const RangeType &factor )
    {
      localFunction().axpy( quadrature, quadPoint, factor );
    }
#endif

    /** \copydoc Dune::LocalFunctionInterface::axpy(const PointType &x,const JacobianRangeType &factor) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const int quadPoint,
                       const JacobianRangeType &factor )
    {
      localFunction().axpy( x, factor );
    }
   
#if 0
    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &factor) */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const JacobianRangeType &factor )
    {
      localFunction().axpy( quadrature, quadPoint, factor );
    }
#endif
    
     /** \copydoc Dune::LocalFunctionInterface::axpy(const PointType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class PointType >
    inline void axpy ( const PointType &x,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 )
    {
      localFunction().axpy( x, factor1, factor2 );
    }
  
#if 0
    /** \copydoc Dune::LocalFunctionInterface::axpy(const QuadratureType &quadrture,const int quadPaint,const RangeType &factor1,const JacobianRangeType &factor2) */
    template< class QuadratureType >
    inline void axpy ( const QuadratureType &quadrature,
                       const int quadPoint,
                       const RangeType &factor1,
                       const JacobianRangeType &factor2 )
    {
      localFunction().axpy( quadrature, quadPoint, factor1, factor2 );
    }
#endif
    
    /** \copydoc LocalFunctionInterface::baseFunctionSet() const */
    const BaseFunctionSetType &baseFunctionSet() const 
    {
      return localFunction().baseFunctionSet();
    }

  private:
    const WrappedLocalFunctionType &localFunction () const 
    { 
      return lf_;
    } 

    WrappedLocalFunctionType &localFunction () 
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
