#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_GENERIC_INLINE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_GENERIC_INLINE_HH

#include "generic.hh"

namespace Dune
{

  namespace Fem
  {

    // GenericLocalFunctionImpl
    // -------------------------

    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::GenericLocalFunctionImpl ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction ),
      values_( DiscreteFunctionSpace::localBlockSize * discreteFunction_.space().blockMapper().maxNumDofs() ),
      basisFunctionSet_(),
      entity_( 0 ),
      numDofs_( 0 )
    {}
   


    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::GenericLocalFunctionImpl ( const ThisType &other )
    : discreteFunction_( other.discreteFunction_ ),
      values_( other.values_ ),
      basisFunctionSet_( other.basisFunctionSet_ ),
      entity_( other.entity_ ),
      numDofs_( other.numDofs_ )
    {}



    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline
    const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::DofType &
    GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::operator[] ( const int num ) const
    {
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }
    

    
    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline
    typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::DofType &
    GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::operator[] ( const int num )
    {
      assert( (num >= 0) && (num < numDofs()) );
      return *(values_[ num ]);
    }



    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline
    const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::BaseFunctionSetType &
    GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::basisFunctionSet () const
    {
      assert( entity_ != 0 );
      return basisFunctionSet_;
    }


    
    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline
    const typename GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::EntityType &
    GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::entity () const
    {
      assert( entity_ != 0 );
      return *entity_;
    }



    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline void GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >::init ( const EntityType &entity )
    {
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, DiscreteFunctionSpaceType :: localBlockSize > NonBlockMapper;
     
      const DiscreteFunctionSpaceType &space = discreteFunction_.space();
      
      basisFunctionSet_ = space.basisFunctionSet( entity );
      numDofs_ = basisFunctionSet_.size();
      entity_ = &entity;

      assert( numDofs_ <= values_.size() );
      
      const BlockMapperType &blockMapper = space.blockMapper();
      NonBlockMapper mapper( blockMapper );
      mapper.map( entity, indices_ );

      assert( indices_.size() == numDofs_ );
      for( unsigned int i = 0; i < numDofs_; ++i )
        values_[ i ] = &discreteFunction_.dof( indices_[ i ] );
    }



    template< class DiscreteFunction, class DiscreteFunctionSpace >
    inline int GenericLocalFunctionImpl< DiscreteFunction, DiscreteFunctionSpace >
      ::numDofs () const
    {
      return numDofs_;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_GENERIC_INLINE_HH
