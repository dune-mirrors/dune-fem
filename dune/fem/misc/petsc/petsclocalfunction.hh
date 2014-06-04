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
    template< class DiscreteFunctionTraits > class PetscLocalFunction;

    /* ========================================
     * class PetscLocalFunctionFactory
     * =======================================
     */
    template< class DiscreteFunctionTraits >
    class PetscLocalFunctionFactory
    {
      typedef PetscLocalFunctionFactory< DiscreteFunctionTraits > ThisType;

    public:
      typedef typename DiscreteFunctionTraits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionTraits :: DiscreteFunctionType      DiscreteFunctionType;
      typedef PetscLocalFunction< DiscreteFunctionType >                   ObjectType;

      explicit PetscLocalFunctionFactory ( DiscreteFunctionType &dFunction )
      : discreteFunction_( dFunction )
      {}

      // TODO: implement this!
      ObjectType* newObject () const { return new ObjectType( discreteFunction_ ); }

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
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

      //! type of local coordinates 
      typedef typename EntityType::Geometry::LocalCoordinate  LocalCoordinateType;

    private:
      // type of Dof, here it's a proxy object because of the access to PetscVec
      typedef typename DiscreteFunctionType::DofBlockType::DofProxy DofProxyType;
      typedef std::vector< DofProxyType > ProxyVectorType;

      struct AssignDofs;
    public:
      // type of Dofs returned by operator [] 
      typedef DofProxyType  DofType;

      enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
      enum { dimRange = DiscreteFunctionSpaceType::dimRange };

      explicit PetscLocalFunction ( DiscreteFunctionType &dFunction )
      : discreteFunction_( dFunction ),
        proxyVector_( DiscreteFunctionSpaceType::localBlockSize * discreteFunction_.space().blockMapper().maxNumDofs() ),
        entity_( 0 ),
        basisFunctionSet_(),
        numDofs_( 0 ),
        needCheckGeometry_( true )
      {}

      PetscLocalFunction ( const PetscLocalFunction &other )
      : discreteFunction_( other.discreteFunction_ ),
        proxyVector_( other.proxyVector_ ),
        entity_( other.entity_ ),
        basisFunctionSet_( other.basisFunctionSet_ ),
        numDofs_( other.numDofs_ ),
        needCheckGeometry_( other.needCheckGeometry_ )
      {}

      const BasisFunctionSetType &basisFunctionSet () const 
      {
        return basisFunctionSet_;
      }

      void init ( const EntityType &entity )
      {
        const DiscreteFunctionSpaceType &space = discreteFunction().space();

        // store basis function set and entity 
        basisFunctionSet_ = space.basisFunctionSet( entity );
        numDofs_ = basisFunctionSet_.size();
        entity_  = &entity;

        assert( basisFunctionSet_.type() == entity.type() );
        // get corresponding dofs 
        space.blockMapper().mapEach( entity, AssignDofs( discreteFunction_, proxyVector_ ) );
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
      ProxyVectorType proxyVector_;
      const EntityType *entity_;
      BasisFunctionSetType basisFunctionSet_;
      unsigned int numDofs_;
      bool needCheckGeometry_;
    };


    // PetscLocalFunctionImpl::AssignDofs
    // -------------------------------------

    template< class DiscreteFunction >
    struct PetscLocalFunction< DiscreteFunction >::AssignDofs
    {
      AssignDofs ( DiscreteFunctionType &discreteFunction, ProxyVectorType &values )
      : discreteFunction_( discreteFunction ), values_( values )
      {}

      template < class GlobalKey >
      void operator () ( const std::size_t local, const GlobalKey& globalKey )
      {
        typedef typename DiscreteFunctionType::DofBlockPtrType DofBlockPtrType;
        static const unsigned int blockSize = DiscreteFunctionSpaceType::localBlockSize;
      
        DofBlockPtrType blockPtr = discreteFunction_.block( globalKey );
        const unsigned int localBlock = local * blockSize;
        for( unsigned int i = 0; i < blockSize; ++i )
          values_[ localBlock + i ].assign( (*blockPtr)[ i ] );
      }

    private:
      DiscreteFunctionType &discreteFunction_;
      ProxyVectorType &values_;
    };

  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_PETSCLOCALFUNCTION_HH
