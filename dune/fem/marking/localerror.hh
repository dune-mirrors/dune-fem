#ifndef DUNE_FEM_MARKING_LOCALERROR_HH
#define DUNE_FEM_MARKING_LOCALERROR_HH

#include <dune/common/dynvector.hh>

#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionalError
    // --------------------

    template< class ErrorFunctional >
    class LocalFunctionalError
    {
      typedef LocalFunctionalError< ErrorFunctional > ThisType;

      static_assert( ErrorFunctional::dimRange == 1, "Error functionals must have dimRange == 1." );

    public:
      typedef ErrorFunctional ErrorFunctionalType;

      typedef typename ErrorFunctionalType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename ErrorFunctionalType::GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

      typedef ConstLocalFunction< ErrorFunctionalType > LocalErrorFunctionalType;

    private:
      struct LocalOne
      {
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

        typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

        static const int dimDomain = FunctionSpaceType::dimDomain;
        static const int dimRange = FunctionSpaceType::dimRange;

        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

        LocalOne ( const EntityType &entity, int order ) : entity_( entity ), order_( order ) {}

        template< class Point >
        void evaluate ( const Point &x, RangeType &value ) const
        {
          value[ 0 ] = 1;
        }

        template< class Quadrature, class Values >
        void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
        {
          for( const auto qp : quadrature )
            evaluate( qp, values[ qp.index() ] );
        }

        int order () const { return order_; }

        const EntityType &entity () const { return entity_; }

      private:
        const EntityType &entity_;
        int order_;
      };

    public:
      explicit LocalFunctionalError ( const ErrorFunctionalType &errorFunctional )
        : localErrorFunctional_( errorFunctional )
      {
        localIndicator_.reserve( space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize );
      }

      RangeFieldType operator() ( const ElementType &element ) const
      {
        localErrorFunctional_.init( element );
        localIndicator_.resize( localErrorFunctional_.basisFunctionSet().size() );

        const auto &interpolation = space().interpolation( element );
        interpolation( LocalOne( element, localErrorFunctional_.order() ), localIndicator_ );

        return localErrorFunctional_.localDofVector() * localIndicator_;
      }

      const DiscreteFunctionSpaceType &space () const { return localErrorFunctional_.discreteFunction().space(); }

    private:
      LocalErrorFunctionalType localErrorFunctional_;
      Dune::DynamicVector< RangeFieldType > localIndicator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MARKING_LOCALERROR_HH
