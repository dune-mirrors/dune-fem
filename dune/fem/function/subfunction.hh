#ifndef DUNE_FEM_SUBFUNCTION_HH
#define DUNE_FEM_SUBFUNCTION_HH

#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/storage/subarray.hh>
#include <dune/fem/function/vectorfunction/vectorfunction.hh>

namespace Dune {

//  namespace Fem {

    template <class DiscreteFunctionImp>   
    class SubFunctionStorage
    {
      SubFunctionStorage( const SubFunctionStorage& );
    protected:  
      typedef DiscreteFunctionImp DiscreteFunctionType;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  SpaceType;
      enum { dimRange = SpaceType :: dimRange };
      typedef typename SpaceType :: DofStorageType            DofStorageType;
    public:  
      typedef typename SpaceType :: template ToNewDimRange < 1 > :: Type  SubSpaceType;  
      typedef typename CombinedSpace< SubSpaceType, dimRange >            CombinedSpaceType;

      typedef CombinedSubMapper< typename SubSpaceType :: MapperType , dimRange, PointBased >  SubMapperType;
      typedef SubVector< DofStorageType, SubMapperType >                  SubDofVectorType;
      typedef VectorDiscreteFunction< SubSpaceType, SubDofVectorType >    SubDiscreteFunctionType;

      explicit SubFunctionStorage( const DiscreteFunctionType& discreteFunction ) 
        discreteFunction_( discreteFunction ),
        space_( discreteFunction.space() ),
        containedSpace_( space_.gridPart (), 
                         space_.communicationInterface(),
                         space_.communicationDirection() ),
        subMapper_( dimRange, (SubMapperType *) 0 ),
        subVector_( dimRange, (SubDofVectorType *) 0 ),
        subDiscreteFunction_( dimRange, (SubDiscreteFunctionType *) 0 )
      {}

      ~SubFunctionStorage() 
      {
        for(int i=0; i<dimRange; ++i)
        {
          delete subDiscreteFunction_[ i ]; subDiscreteFunction_[ i ] = 0;
          delete subVector_[ i ];           subVector_[ i ] = 0;
          delete subMapper_[ i ];           subMapper_[ i ] = 0;
        }
      }

      SubDiscreteFunctionType& subFunction(const size_t component) const 
      {
        assert( component < dimRange );
        if( ! subDiscreteFunction_[ component ] )
        {
          subMapper_[ component ] = new SubMapperType( containedSpace_.mapper(), component );
          subVector_[ component ] = new SubDofVectorType( discreteFunction_.dofStorage(), 
                                                          *subMapper_[component] );
          subDiscreteFunction_[ component ] = 
            new SubDiscreteFunctionType( std::string(discreteFunction_.name()+ "_sub"),
                                         containedSpace_, *( subVector_[ component ] ));
        }
        return *( subDiscreteFunction_[ component ] );
      }

    };



//  }
}
#endif
