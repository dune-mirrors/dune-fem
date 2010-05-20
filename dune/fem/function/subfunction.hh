#ifndef DUNE_FEM_SUBFUNCTION_HH
#define DUNE_FEM_SUBFUNCTION_HH

#include <vector>

#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/storage/subarray.hh>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>

namespace Dune {

//  namespace Fem {

    //! @ingroup SubDFunction
    //! A class for extracting sub functions from a 
    //! discrete function containing pointbased combined data.
    template <class DiscreteFunctionImp>   
    class SubFunctionStorage
    {
      SubFunctionStorage( const SubFunctionStorage& );
    protected:  
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  SpaceType;
      enum { dimRange = SpaceType :: dimRange };
      typedef typename DiscreteFunctionType :: DofStorageType            DofStorageType;
    public:  
      typedef typename SpaceType :: template ToNewDimRange < 1 > :: Type  SubSpaceType;  

      typedef CombinedSubMapper< typename SubSpaceType :: MapperType , dimRange, PointBased >  SubMapperType;
      typedef SubVector< DofStorageType, SubMapperType >                  SubDofVectorType;
      typedef VectorDiscreteFunction< SubSpaceType, SubDofVectorType >    SubDiscreteFunctionType;

      //! constructor storing the discrete function 
      explicit SubFunctionStorage( DiscreteFunctionType& discreteFunction ) :
        discreteFunction_( discreteFunction ),
        space_( discreteFunction.space() ),
        subSpace_( space_.gridPart (), 
                         space_.communicationInterface(),
                         space_.communicationDirection() ),
        subMapper_( dimRange, (SubMapperType *) 0 ),
        subVector_( dimRange, (SubDofVectorType *) 0 ),
        subDiscreteFunction_( dimRange, (SubDiscreteFunctionType *) 0 )
      {}

      //! destructor 
      ~SubFunctionStorage() 
      {
        for(int i=0; i<dimRange; ++i)
        {
          delete subDiscreteFunction_[ i ]; subDiscreteFunction_[ i ] = 0;
          delete subVector_[ i ];           subVector_[ i ] = 0;
          delete subMapper_[ i ];           subMapper_[ i ] = 0;
        }
      }

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
          subMapper_[ component ] = new SubMapperType( subSpace_.mapper(), component );
          subVector_[ component ] = new SubDofVectorType( discreteFunction_.dofStorage(), 
                                                          *subMapper_[component] );
          subDiscreteFunction_[ component ] = 
            new SubDiscreteFunctionType( std::string(discreteFunction_.name()+ "_sub"),
                                         subSpace_, *( subVector_[ component ] ));
        }
        return *( subDiscreteFunction_[ component ] );
      }

    protected:
      DiscreteFunctionType& discreteFunction_;
      const SpaceType& space_;
      SubSpaceType subSpace_;
      mutable std::vector< SubMapperType * > subMapper_;
      mutable std::vector< SubDofVectorType * > subVector_;
      mutable std::vector< SubDiscreteFunctionType* > subDiscreteFunction_;
    };

//  }
}
#endif
