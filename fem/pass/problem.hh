#ifndef DUNE_PROBLEM_HH
#define DUNE_PROBLEM_HH

namespace Dune {

  template <class ProblemImp, class FunctionSpaceImp>
  class ProblemInterface 
  {
  public:
    typedef ProblemImp ProblemType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    
  public:
    bool hasFlux() const { return asImp().hasFlux(); }
    bool hasSource() const { return asImp().hasSource(); }

    // Rest of the interface: what goes in and what comes out?
    void numericalFlux() { asImp().numericalFlux(); }

    void analyticalFlux() { asImp().analyticalFlux(); }

    void source() { asImp().source(); }

  private:
    ProblemType& asImp() { return static_cast<ProblemType&>(*this); }
    const ProblemType& asImp() const { return static_cast<const ProblemType&>(*this); }
  };
  
  template <class ProblemImp, class FunctionSpaceImp>
  class ProblemDefault : 
    public ProblemInterface<ProblemImp, FunctionSpaceImp> {
  public:
    // Here, implement standard methods to do nothing
    void numericalFlux() {}
    void analyticalFlux() {}
    void source() {}
  };

}  // end namespace Dune
#endif
