#ifndef DUNE_FEM_EIGENOPERATOR_HH
#define DUNE_FEM_EIGENOPERATOR_HH

#ifdef HAVE_EIGEN

#include <dune/fem/operator/matrix/eigenmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // EigenLinearOperator
    // -----------------------

    template< class DomainFunction, class RangeFunction >
    class EigenLinearOperator
    : public EigenMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef EigenMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > Base;

    public:
      typedef typename Base::DomainSpaceType DomainSpaceType;
      typedef typename Base::RangeSpaceType RangeSpaceType;

      /** \copydoc Fem::Operator::assembled */
      static const bool assembled = true ;

      using Base::apply;

      EigenLinearOperator ( const std::string &name,
                                const DomainSpaceType &domainSpace,
                                const RangeSpaceType &rangeSpace,
                                const std::string &paramfile = "" )
      : Base( domainSpace, rangeSpace, paramfile )
      {}

      virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const
      {
        Base::apply( arg, dest );
      }

      const Base &systemMatrix () const
      {
        return *this;
      }

      void communicate () const
      {
      }
    };
  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef DUNE_FEM_SPLINEAR_HH
