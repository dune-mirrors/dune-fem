#ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
#define DUNE_FEM_ISTLLINEAROPERATOR_HH

#if HAVE_DUNE_ISTL

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/istlmatrix.hh>

namespace Dune
{

  namespace Fem
  {
    template< class DomainFunction, class RangeFunction >
    struct ISTLLinearOperator;

    //! ISTLMatrixOperator (any discrete function, ISTL spec below)
    template< class DomainFunction, class RangeFunction >
    struct ISTLLinearOperator
    : public ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef ISTLLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef ISTLMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;

      //! constructor
      //! \param domainSpace space defining domain of operator
      //! \param rangeSpace  space defining range of operator
      //! \param param ISTL matrix parameters for preconditioning
      //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
      //!         - Pre-iteration: number of iteration of preconditioner
      //!         - Pre-relaxation: relaxation factor
      ISTLLinearOperator( const std::string & ,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace,
                          const ISTLSolverParameter& param = ISTLSolverParameter() )
        : BaseType( domainSpace, rangeSpace, param )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      virtual void finalize() { BaseType::compress(); }

    };

    //! ISTLMatrixOperator
    template< class DomainFunctionSpace, class RangeFunctionSpace,
              class DomainBlock, class RangeBlock >
    struct ISTLLinearOperator< ISTLBlockVectorDiscreteFunction< DomainFunctionSpace, DomainBlock >,
                               ISTLBlockVectorDiscreteFunction< RangeFunctionSpace, RangeBlock > >
    : public ISTLMatrixObject< DomainFunctionSpace, RangeFunctionSpace, DomainBlock, RangeBlock >,
      public AssembledOperator< ISTLBlockVectorDiscreteFunction< DomainFunctionSpace, DomainBlock >,
                                ISTLBlockVectorDiscreteFunction< RangeFunctionSpace, RangeBlock > >
    {
      typedef DomainFunctionSpace DomainSpaceType;
      typedef RangeFunctionSpace  RangeSpaceType;

      typedef ISTLBlockVectorDiscreteFunction< DomainFunctionSpace, DomainBlock > DomainFunction;
      typedef ISTLBlockVectorDiscreteFunction< RangeFunctionSpace, RangeBlock >   RangeFunction;

      typedef ISTLLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef ISTLMatrixObject< DomainSpaceType, RangeSpaceType, DomainBlock, RangeBlock > BaseType;

      static constexpr bool assembled = true;

      using BaseType::apply;

      //! constructor
      //! \param domainSpace space defining domain of operator
      //! \param rangeSpace  space defining range of operator
      //! \param param ISTL matrix parameters for preconditioning
      //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
      //!         - Pre-iteration: number of iteration of preconditioner
      //!         - Pre-relaxation: relaxation factor
      ISTLLinearOperator( const std::string & ,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace,
                          const ISTLSolverParameter& param = ISTLSolverParameter() )
        : BaseType( domainSpace, rangeSpace, param )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      virtual void finalize() { BaseType::compress(); }

    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
