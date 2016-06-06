#ifndef DUNE_FEM_SUBFUNCTION_HH
#define DUNE_FEM_SUBFUNCTION_HH

#include <memory>
#include <vector>

#include <dune/fem/space/combinedspace/combineddofstorage.hh>
#include <dune/fem/storage/subvector.hh>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>

namespace Dune
{

  namespace Fem
  {

    //! @ingroup SubDFunction
    //! A class for extracting sub functions from a
    //! discrete function containing pointbased combined data.
    template <class DiscreteFunctionImp>
    class SubFunctionStorage
    {
    protected:
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  SpaceType;
      enum { dimRange = SpaceType :: dimRange };
      typedef typename DiscreteFunctionType :: DofStorageType            DofStorageType;
    public:
      typedef typename SpaceType :: template ToNewDimRange < 1 > :: Type  SubSpaceType;

      typedef CombinedSubMapper< typename SubSpaceType :: MapperType , dimRange, PointBased >  SubMapperType;
      typedef Fem :: SubVector< DofStorageType, SubMapperType >                  SubDofVectorType;
      typedef VectorDiscreteFunction< SubSpaceType, SubDofVectorType >    SubDiscreteFunctionType;

      //! constructor storing the discrete function
      explicit SubFunctionStorage( DiscreteFunctionType& discreteFunction ) :
        discreteFunction_( discreteFunction ),
        space_( discreteFunction.space() ),
        subSpace_( space_.gridPart (),
                         space_.communicationInterface(),
                         space_.communicationDirection() ),
        subMapper_( dimRange, nullptr ),
        subVector_( dimRange, nullptr ),
        subDiscreteFunction_( dimRange, nullptr )
      {}

       SubFunctionStorage( const SubFunctionStorage& ) = delete;

      /** \brief return a SubDiscreteFunction repsenting only one
          component of the original discrete function
          \param component the component to be extracted

          \return reference to SubDiscreteFunction for given component
      */
      SubDiscreteFunctionType& subFunction(const size_t component) const
      {
        assert( component < dimRange );
        if( ! subDiscreteFunction_[ component ] )
        {
          subMapper_[ component ] = std::make_unique< SubMapperType >( subSpace_.mapper(), component );
          subVector_[ component ] = std::make_unique< SubDofVectorType >( discreteFunction_.dofStorage(),
                                                                          *subMapper_[component] );
          subDiscreteFunction_[ component ] =
            std::make_unique< SubDiscreteFunctionType >( std::string(discreteFunction_.name()+ "_sub"),
                                                         subSpace_, *( subVector_[ component ] ));
        }
        return *( subDiscreteFunction_[ component ] );
      }

    protected:
      DiscreteFunctionType& discreteFunction_;
      const SpaceType& space_;
      SubSpaceType subSpace_;
      mutable std::vector< std::unique_ptr< SubMapperType > > subMapper_;
      mutable std::vector< std::unique_ptr< SubDofVectorType > > subVector_;
      mutable std::vector< std::unique_ptr< SubDiscreteFunctionType > > subDiscreteFunction_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SUBFUNCTION_HH
