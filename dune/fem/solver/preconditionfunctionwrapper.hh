#ifndef DUNE_FEM_SOLVER_PRECONDITIONFUNCTIONWRAPPER_HH
#define DUNE_FEM_SOLVER_PRECONDITIONFUNCTIONWRAPPER_HH

#include <functional>

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{
  namespace Fem
  {
    /** \class PreconditionerFunctionWrapper
     *  \brief Wrapper for functions passed from Python side that implements a
     *  preconditioner.
     *
     *  \tparam DomainFunction argument function.
     *  \tparam RangeFunction  destination function.
     *
     */
    template <class DomainFunction, class RangeFunction = DomainFunction>
    class PreconditionerFunctionWrapper
      : public virtual Operator< DomainFunction, RangeFunction >
    {
    public:
      typedef DomainFunction  DomainFunctionType;
      typedef RangeFunction   RangeFunctionType;

      typedef std::reference_wrapper< const DomainFunctionType >  ConstDomainDFType;
      typedef std::reference_wrapper< RangeFunctionType >         RangeDFType;
      typedef std::function< void( ConstDomainDFType& , RangeDFType& ) > PreconditionerFunctionType;

    protected:
      const PreconditionerFunctionType& preconditioner_; // function given from Python side
    public:
      PreconditionerFunctionWrapper( const PreconditionerFunctionType& pre )
        : preconditioner_( pre ) {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &v ) const final override
      {
        // convert to reference_wrapper to avoid copying
        ConstDomainDFType uR( u );
        RangeDFType vR( v );

        // callback to python applying preconditioner
        preconditioner_( uR, vR );
      }
    };

  } // end namespace Fem

} // end namespace Dune
#endif
