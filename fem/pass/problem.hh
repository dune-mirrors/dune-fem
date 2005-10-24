#ifndef DUNE_PROBLEM_HH
#define DUNE_PROBLEM_HH

namespace Dune {

  /**
   * @brief Interface class for problem definition in the LDG context.
   *
   * The problem interface prescribes the methods needed by the implementation
   * of local DG passes. Users need to derive from either this class or the
   * ProblemDefault class.
   */
  template <class ProblemImp, class FunctionSpaceImp>
  class ProblemInterface 
  {
  public:
    typedef ProblemImp ProblemType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    
  public:
    //! Returns true if problem has a flux contribution.
    //! If you program this method to return true, make sure to implement
    //! numericalFlux and analyticalFlux as well.
    bool hasFlux() const { return asImp().hasFlux(); }

    //! Returns true if problem has a source term.
    //! If you program this method to return true, make sure to implement
    //! method source well.
    bool hasSource() const { return asImp().hasSource(); }

    //! Computes the numerical flux at a cell interface.
    //! ... More to come ...
    template <
      class IntersectionIterator, class ArgumentTuple, class ResultType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const DomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         ResultType& gLeft,
                         ResultType& gRight)
    { return asImp().numericalFlux(); }

    //! Computes the analytical flux of the problem.
    //! ... More to come ...
    template <class Entity, class ArgumentTuple, class ResultType>
    void analyticalFlux(Entity& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, ResultType& f) 
    { asImp().analyticalFlux(); }

    //! Implements the source term of the problem.
    template <class Entity, class ArgumentTuple, class ResultType>
    void source(Entity& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, ResultType& s) 
    { asImp().source(); }

    // Need second source function with grad u as additional argument 

  private:
    ProblemType& asImp() { return static_cast<ProblemType&>(*this); }
    const ProblemType& asImp() const { return static_cast<const ProblemType&>(*this); }
  };
  
  //! Default implementation of the ProblemInterface where methods for 
  //! the fluxes and the source term do nothing, so that the user needn't
  //! implement them if not needed.
  template <class ProblemImp, class FunctionSpaceImp>
  class ProblemDefault : 
    public ProblemInterface<ProblemImp, FunctionSpaceImp> {
  public:
    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
   template <
      class IntersectionIterator, class ArgumentTuple, class ResultType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const DomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         ResultType& gLeft,
                         ResultType& gRight)
    { 
      assert(!this->hasFlux()); 
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class Entity, class ArgumentTuple, class ResultType>
    void analyticalFlux(Entity& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, ResultType& f)
    { assert(!this->hasFlux()); }

    //! Empty implementation that fails if problem claims to have a source 
    //! term.
    template <class Entity, class ArgumentTuple, class ResultType>
    void source(Entity& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, ResultType& s)
    { assert(!this->hasSource()); }
  };

}  // end namespace Dune
#endif
