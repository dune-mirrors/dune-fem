namespace Dune
{

  ////////////////////////////////////////////////////////
  //- AdaptiveDiscreteFunction (specialisation)
  ////////////////////////////////////////////////////////
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  inline AdaptiveDiscreteFunction<
     CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  ~AdaptiveDiscreteFunction() 
  {
    for (int i=0; i<N; ++i) 
    {
      delete subDiscFunc_[i]; subDiscFunc_[i] = 0;
      delete subDofVector_[i]; subDofVector_[i] = 0;
      delete subDofMapper_[i]; subDofMapper_[i] = 0;
    }
  }
  
  template <class ContainedFunctionSpaceImp, int N, DofStoragePolicy p>
  inline void AdaptiveDiscreteFunction<
     CombinedSpace<ContainedFunctionSpaceImp, N, p> >::
  initializeSubFunctions() 
  {
    const SubSpaceType& subSpace = this->space().containedSpace();
    for (int i=0;i<N;i++)
    {
      subDofMapper_[i] = new SubMapperType(this->spc_,i);
      subDofVector_[i] = new SubDofVectorType(this->dofStorage(), *subDofMapper_[i]);
      subDiscFunc_[i]  = new SubDiscreteFunctionType(
                               std::string(this->name() + "_sub"),
                               subSpace,*(subDofVector_[i]));
    }
  }
} // end namespace Dune
