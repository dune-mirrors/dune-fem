#ifndef DUNE_PROBLEM_HH
#define DUNE_PROBLEM_HH

#include <dune/common/fvector.hh>

namespace Dune {

  /**
   * @brief Interface class for problem definition in the LDG context.
   *
   * The problem interface prescribes the methods needed by the implementation
   * of local DG passes. Users need to derive from either this class or the
   * ProblemDefault class. The methods provided by this class are used by
   * LocalDGPass (and possibly other classes) to solve for the space 
   * discretisation of the equation d_t u + div f(u) = s(u), where f
   * represents a flux function and s a source term. Even non-conservative
   * fluxes can be treated (for details see paper by Dedner et al).
   * \note The definition of f(u) differs from the usual definition for 
   * conservation laws in the leading sign.
   */
  template <class ProblemTraits>
  class ProblemInterface 
  {
  public:
    //! Traits class defined by the user
    typedef ProblemTraits Traits;
    //! Implementation type for Barton-Nackman trick
    typedef typename Traits::ProblemType ProblemType;
    //! Function space (we take the discrete one even though the continuous one
    //! would suffice...)
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Coordinate type (world coordinates)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! Vector type of the function space's range
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! Matrix type of the function space's jacobian
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! Type of the grid
    typedef typename Traits::GridType GridType;
    //! Intersection iterator of the grid
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    //! Element (codim 0 entity) of the grid
    typedef typename GridType::template Codim<0>::Entity EntityType;

  public:
    //! Returns true if problem has a flux contribution.
    //! If you program this method to return true, make sure to implement
    //! numericalFlux and analyticalFlux as well.
    bool hasFlux() const { return asImp().hasFlux(); }

    //! Returns true if problem has a source term.
    //! If you program this method to return true, make sure to implement
    //! method source well.
    bool hasSource() const { return asImp().hasSource(); }

    //! \brief Computes the numerical flux at a cell interface.
    //! Interface method for calculating the numerical flux at cell interfaces.
    //! The method has distinct return values for the contribution to the inner
    //! and outer element, since the framework allows for non-conservative
    //! terms. The non-conservative terms don't get evaluated over separate
    //! methods, but need to be implemented using the flux and source methods.
    //! \param it Iterator of the face in consideration.
    //! \param time Global time.
    //! \param x Local coordinate of the point where the numerical flux gets
    //! evaluated.
    //! \param uLeft Tuple of the states on the inner side of the intersection.
    //! \param uRight Tuple of the states on the outer side of the
    //! intersection.
    //! \param gLeft The numerical flux contribution to the inner element 
    //! (return value).
    //! \param gRight The numerical flux contribution to the outer element
    //! (return value).
    //! \return Speed of the fastest wave. This information is used to compute
    //! the maximum admissible timestep size.
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { return asImp().numericalFlux(); }

    //! \brief Computes the flux at the boundary
    //! Special kind of numerical flux. The intersection iterator provides
    //! the necessary information to identify the type of boundary.
    //! \param it Iterator of the face in consideration.
    //! \param time Global time.
    //! \param x Local coordinate of the point where the numerical flux gets
    //! evaluated.
    //! \param uLeft Tuple of the states on the inner side of the intersection.
    //! \param gLeft The numerical flux contribution to the inner element 
    //! (return value).
    //! \return Speed of the fastest wave. This information is used to compute
    //! the maximum admissible timestep size.
    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    { return asImp().boundaryFlux(); }

    //! \brief Computes the analytical flux of the problem.
    //! Analytical flux of the problem.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param x Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param f The analytical flux (return value)
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f) 
    { asImp().analyticalFlux(); }

    //! \brief Implements the source term of the problem.
    //! Source term contribution. The source term can contain parts of the
    //! non-conservative contributions.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param x Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param jac Tuple of the gradients of the precedings passes' results 
    //! (needed for the non-conservative contributions)
    //! \param s The source contribution (return value).
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, 
                const JacobianTuple& jac, 
                RangeType& s) 
    { asImp().source(); }

  private:
    ProblemType& asImp() { return static_cast<ProblemType&>(*this); }
    const ProblemType& asImp() const { return static_cast<const ProblemType&>(*this); }
  };
  
  //! Default implementation of the ProblemInterface where methods for 
  //! the fluxes and the source term do nothing, so that the user needn't
  //! implement them if not needed.
  template <class ProblemTraits>
  class ProblemDefault : 
    public ProblemInterface<ProblemTraits> {
  public:
    typedef ProblemTraits Traits;
    typedef typename Traits::ProblemType ProblemType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename Traits::GridType GridType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;

  public:
    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      assert(!this->hasFlux()); 
      gLeft = 0.0;
      gRight = 0.0;
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      assert(!this->hasFlux());
      gLeft = 0.0;
      return 0.0;
    }

    //! Empty implementation that fails if problem claims to have a flux
    //! contribution.
    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      assert(!this->hasFlux()); 
      f = 0.0;
    }

    //! Empty implementation that fails if problem claims to have a source 
    //! term.
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, const JacobianTuple& jac, 
                RangeType& s)
    { 
      assert(!this->hasSource()); 
      s = 0.0;
    }
  };

}  // end namespace Dune
#endif
