namespace Dune
{

  ////////////////////////////////////////////////////////
  //- AdaptiveDiscreteFunction (specialisation)
  ////////////////////////////////////////////////////////
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  AdaptiveDiscreteFunction<
     CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveDiscreteFunction() 
  {
    for (int i=0;i<N;i++) {
      delete subDofVector_[i];
      delete subDofMapper_[i];
    }
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  typename AdaptiveDiscreteFunction<
    CombinedSpace<ContainedFunctionSpaceImp, N, p> >::SubDiscreteFunctionType
  AdaptiveDiscreteFunction<CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  subFunction(int component) 
  {
     const SubSpaceType& subSpace = this->space().containedSpace();
     return SubDiscreteFunctionType(std::string("Subfunction of ")+this->name(),
                                    subSpace,*(subDofVector_[component]));
  }

} // end namespace Dune
