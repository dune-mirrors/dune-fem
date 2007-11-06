#ifndef DUNE_ELLIPTICDISCRETEMODEL_HH
#define DUNE_ELLIPTICDISCRETEMODEL_HH


//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/common/bartonnackmanifcheck.hh>

//- local includes 
#include <dune/fem/misc/boundaryidentifier.hh>
#include "selection.hh"
#include "discretemodel.hh"

namespace Dune {

  /**
   * @brief Interface class for problem definition in the 
   * DG primal formulation context.
   *
   */
  template <class DiscreteModelTraits>
  class EllipticDiscreteModelInterface 
  {
  public:
    //! Traits class defined by the user
    typedef DiscreteModelTraits Traits;
    //! Implementation type for Barton-Nackman trick
    typedef typename Traits::DiscreteModelType DiscreteModelType;
    //! Function space (we take the discrete one even though the continuous one
    //! would suffice...)
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Coordinate type (world coordinates)
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! Vector type of the function space's range
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    //! Matrix type of the function space's jacobian
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    //! type of range field 
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    //! Type of GridPart 
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType; 
    //! Type of the grid
    typedef typename GridPartType::GridType GridType;

    //! dimension of grid 
    enum { dim = GridType :: dimension };

    //! type of domain for faces 
    typedef FieldVector<RangeFieldType,dim-1> FaceDomainType; 

    //! Intersection iterator of the grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    //! Element (codim 0 entity) of the grid
    typedef typename GridType::template Codim<0>::Entity EntityType;

    //! type of boundary identifier 
    typedef BoundaryIdentifier BoundaryIdentifierType; 

  public:
    //! Returns true if problem has a coefficient before laplace term.
    bool hasCoefficient () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().hasCoefficientFlux());
      return asImp().hasCoefficientFlux(); 
    }

    //! Returns true if problem has right hand side not qual to zero 
    bool hasRHS () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().hasRHS());
      return asImp().hasRHS(); 
    }

    //! \brief Computes the flux at the boundary
    //! Special kind of numerical flux. The intersection iterator provides
    //! the necessary information to identify the type of boundary.
    //! \param it Iterator of the face in consideration.
    //! \param time Global time.
    //! \param local local point where the boundaryvalue is evaluated
    //! \param u Tuple of the states of the previous passes and of the global argument.
    //! \param bndVal boundary Value (return value)
    template <class FaceDomType,
              class ArgumentTuple>
    BoundaryIdentifierType
    boundaryValue(const IntersectionIteratorType& it,
                  double time, const FaceDomType& local,
                  const ArgumentTuple& u,
                  RangeType& bndVal) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().boundaryValue(it,time,local,u,bndVal));
      return asImp().boundaryValue(it,time,local,u,bndVal);
    }

    template <class ArgumentTuple, class CoefficientType>
    void coefficientFace(const IntersectionIteratorType& it,
                         const double time,
                         const FaceDomainType& local,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         CoefficientType & coeffLeft,
                         CoefficientType & coeffRight) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().coefficientFace(it,time,local,
                        uLeft,uRight,coeffLeft,coeffRight));
    }

 
    //! \brief Implements the source term of the problem.
    //! Source term contribution. The source term can contain parts of the
    //! non-conservative contributions.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param local Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param jac Tuple of the gradients of the precedings passes' results 
    //! (needed for the non-conservative contributions)
    //! \param s The source contribution (return value).
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& local,
                const ArgumentTuple& u, 
                const JacobianTuple& jac, 
                RangeType& s) 
    { 
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().source(en,time,local,u,jac,s));
    }

    //! \brief Implements the right hand side of the problem.
    //! Right hand side term contribution.
    //! \param en Entity where the flux is to be evaluated
    //! \param time Global time.
    //! \param local Reference element coordinate where the flux is evaluated
    //! \param u Tuple of the states of the previous passes and of the global
    //! argument.
    //! \param jac Tuple of the gradients of the precedings passes' results 
    //! (needed for the non-conservative contributions)
    //! \param rhs The right hand side contribution (return value).
    template <class ArgumentTuple, class JacobianTuple>
    void rightHandSide(EntityType &en,
                       const double time,
                       const DomainType& local,
                       const ArgumentTuple& u,
                       const JacobianTuple& jac,
                       RangeType& rhs)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().rightHandSide(en,time,local,u,jac,rhs));
    }

  protected:
    // don't create instances from this class 
    EllipticDiscreteModelInterface() {}

  private:
    DiscreteModelType& asImp() { return static_cast<DiscreteModelType&>(*this); }
    const DiscreteModelType& asImp() const { return static_cast<const DiscreteModelType&>(*this); }
  };
  
}  // end namespace Dune

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION

#endif
