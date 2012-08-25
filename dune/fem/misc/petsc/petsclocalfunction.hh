// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCLOCALFUNCTION_HH
#define DUNE_FEM_PETSCLOCALFUNCTION_HH

#include <vector>

#include <dune/fem/misc/petsc/petsccommon.hh>

#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune 
{

  namespace Fem 
  {

    /*
     * forward declarations
     */
    template< typename DF > class PetscLocalFunction;

    /* ========================================
     * class PetscLocalFunctionFactory
     * =======================================
     */
    template< typename DF >
    class PetscLocalFunctionFactory
    {
      typedef PetscLocalFunctionFactory< DF > ThisType;

    public:
      typedef DF DiscreteFunctionType;
      typedef PetscLocalFunction< DF > ObjectType;

      explicit PetscLocalFunctionFactory ( DiscreteFunctionType &dFunction )
      : discreteFunction_( dFunction )
      {}

      // TODO: implement this!
      ObjectType* newObject () const { assert( false ); abort(); return 0; }

      DiscreteFunctionType& discreteFunction () { return discreteFunction_; }

      const DiscreteFunctionType& discreteFunction () const { return discreteFunction_; }

    private:
      const ThisType& operator= ( const ThisType& );

      // data fields
      DiscreteFunctionType &discreteFunction_;
    };

    /* ========================================
     * class PetscLocalFunction
     * =======================================
     */
    template< typename DF >
    class PetscLocalFunction
    : public LocalFunctionDefault< typename DF::DiscreteFunctionSpaceType, PetscLocalFunction< DF > >
    {

      typedef PetscLocalFunction ThisType;
      typedef LocalFunctionDefault< typename DF::DiscreteFunctionSpaceType, ThisType > BaseType;

    public:

      typedef DF DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    private:
      typedef typename DiscreteFunctionType::DofBlockType::DofProxy DofProxyType;
      typedef std::vector< DofProxyType > ProxyVectorType;

    public:

      enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
      enum { dimRange = DiscreteFunctionSpaceType::dimRange };

      explicit PetscLocalFunction ( DiscreteFunctionType &dFunction )
      : discreteFunction_( dFunction ),
        numDofs_( 0 ),
        proxyVector_( DiscreteFunctionSpaceType::localBlockSize * discreteFunction_.space().blockMapper().maxNumDofs() ),
        needCheckGeometry_( true ),
        entity_( 0 )
      {}

      PetscLocalFunction ( const PetscLocalFunction &other )
      : discreteFunction_( other.discreteFunction_ ),
        numDofs_( other.numDofs_ ),
        proxyVector_( other.proxyVector_ ),
        needCheckGeometry_( other.needCheckGeometry_ ),
        entity_( other.entity_ ),
        baseFunctionSet_( other.baseFunctionSet_ )
      {}

      const BaseFunctionSetType &baseFunctionSet () const 
      {
        return baseFunctionSet_;
      }

      void init ( const EntityType &entity )
      {
        typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
        typedef typename BlockMapperType::DofMapIteratorType DofMapIteratorType;
        enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

        typedef typename DiscreteFunctionType::DofBlockPtrType DofBlockPtrType;
        
        const DiscreteFunctionSpaceType &space = discreteFunction().space();
        const bool multipleBaseSets = space.multipleBaseFunctionSets();

        if( multipleBaseSets || needCheckGeometry_ )
        {
          // if multiple base sets skip geometry call
          bool updateBaseSet = true;
          if( !multipleBaseSets && ( entity_ != 0 ) )
            updateBaseSet = ( baseFunctionSet_.geometryType() != entity.type() );
          
          if( multipleBaseSets || updateBaseSet )
          {
            baseFunctionSet_ = space.baseFunctionSet( entity );

            // note, do not use baseFunctionSet() here, entity might no have been set
            numDofs_ = baseFunctionSet_.size();

            needCheckGeometry_ = space.multipleGeometryTypes();
          }
        }

        // cache entity
        entity_ = &entity;
        assert( baseFunctionSet_.geometryType() == entity.type() );

        assert( numDofs_ <= proxyVector_.size() );
        const BlockMapperType &mapper = space.blockMapper();
        const DofMapIteratorType end = mapper.end( entity );
        for( DofMapIteratorType it = mapper.begin( entity ); it != end; ++it )
        {
          assert( it.global() == int( mapper.mapToGlobal( entity, it.local() ) ) );
          
          DofBlockPtrType blockPtr = discreteFunction().block( it.global() );
          
          const unsigned int localBlock = it.local() * blockSize;
          for( unsigned int i = 0; i < blockSize; ++i )
            proxyVector_[ localBlock + i ].assign( (*blockPtr)[ i ] );
        }
      }

      int order() const
      {
        return discreteFunction_.space().order( entity() );
      }

      DofProxyType& operator[] ( unsigned int index )
      {
        assert( index < numDofs() );
        return proxyVector_[ index ];
      }

      const DofProxyType& operator[] ( unsigned int index ) const
      {
        assert( index < numDofs() );
        return proxyVector_[ index ];
      }

      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        baseFunctionSet().jacobianAll( x, entity().geometry().jacobianInverseTransposed( coordinate( x ) ), *this, ret);
      }

      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor )
      {
        baseFunctionSet().axpy( x, factor, *this );
      }
      template< class PointType >
      void axpy ( const PointType &x, const JacobianRangeType &factor)
      {
        baseFunctionSet().axpy( x, entity().geometry().jacobianInverseTransposed( coordinate( x ) ), factor, *this );
      }
      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
      {
        baseFunctionSet().axpy( x, entity().geometry().jacobianInverseTransposed( coordinate( x ) ), factor1, factor2, *this );
      }

      unsigned int numDofs () const { return numDofs_; }

      const EntityType& entity () const
      {
        assert( entity_ );
        return *entity_;
      }

    private:
      PetscLocalFunction ();  
      ThisType& operator= ( const ThisType& );


      // private methods
      DiscreteFunctionType& discreteFunction () { return discreteFunction_; }
      const DiscreteFunctionType& discreteFunction () const { return discreteFunction_; }
      
      /*
       * data fields
       */
      DiscreteFunctionType &discreteFunction_;
      unsigned int numDofs_;
      ProxyVectorType proxyVector_;
      bool needCheckGeometry_;
      const EntityType *entity_;
      BaseFunctionSetType baseFunctionSet_;
      
    };

    /* ========================================
     * class PetscLocalFunctionStack
     * =======================================
     */
    template< typename DF >
    class PetscLocalFunctionStack
    {
      typedef PetscLocalFunctionStack< DF > ThisType;

    public:
      typedef DF DiscreteFunctionType;
      typedef PetscLocalFunctionFactory< DiscreteFunctionType > LocalFunctionFactoryType;
      typedef PetscLocalFunction< DF > LocalFunctionType;

      explicit PetscLocalFunctionStack ( LocalFunctionFactoryType &factory )
      : factory_( factory )
      {}

      LocalFunctionType localFunction () const { return LocalFunctionType( *this ); }

      template< typename EntityType >
      const LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        LocalFunctionType ret( factory().discreteFunction() );
        ret.init( entity );
        return ret;
      }

      template< typename EntityType >
      LocalFunctionType localFunction ( const EntityType &entity )
      {
        LocalFunctionType ret( factory().discreteFunction() );
        ret.init( entity );
        return ret;
      }

    private:
      PetscLocalFunctionStack ();
      const ThisType& operator= ( const ThisType& );

      // private methods
      LocalFunctionFactoryType& factory() { return factory_; }

      const LocalFunctionFactoryType& factory() const { return factory_; }

      // data fields
      LocalFunctionFactoryType &factory_;
      
    };


  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_PETSCLOCALFUNCTION_HH
