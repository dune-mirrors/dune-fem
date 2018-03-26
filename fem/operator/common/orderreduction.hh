#ifndef DUNE_FEM_DG_ORDERREDUCTION_HH
#define DUNE_FEM_DG_ORDERREDUCTION_HH

#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>

#include <dune/fem/function/tuplediscretefunction.hh>

namespace Dune {

  namespace Fem {

    /** \brief OrderReduction reduces the polynomial order of a discrete
     *         function by projecting from P_E -> P-1_E on each element.
     *
     *  \note  This is only implemented for hierarchic spaces currently.
     */
    template <class DF>
    class OrderReduction
    {
    public:
      typedef DF  DomainFunctionType;
      typedef DF  RangeFunctionType;

      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef DomainSpaceType RangeSpaceType;

      //static_assert(Dune::Fem::Capabilities::isHierarchic< DomainSpaceType > :: v, "OrderReduction is only implemented for hierarchic discrete function spaces");

      OrderReduction( const DomainSpaceType& )
      {
      }


      OrderReduction()
      {
      }

      void operator () (const DomainFunctionType& arg, RangeFunctionType& dest )
      {
        typedef typename DomainFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        const auto& basisSets = arg.space().basisFunctionSets();

        for( const auto& entity : arg.space() )
        {
          const auto hlf = arg.localFunction( entity );
          const int numDofs = hlf.numDofs();
          const int order   = hlf.order();

          const int lowerOrder = std::max( order-1, 0 );

          /*
          for( int o=0; o<order; ++o )
          {
            // get basis function set with one order lower
            const auto& base = basisSets.basisFunctionSet( entity, o );

            std::cout << "baseSet order = " << order << " " << o << std::endl;
            std::cout << "baseSet start = " << base.size() << " " << numDofs << std::endl;
          }
          */

          // get basis function set with one order lower
          const auto& lowerBase = basisSets.basisFunctionSet( entity, lowerOrder );

          auto llf = dest.localFunction( entity );
          // copy dofs first
          llf.assign( hlf );

          //std::cout << "baseSet order = " << order << " " << lowerOrder << std::endl;
          //std::cout << "baseSet start = " << lowerBase.size() << " " << numDofs << std::endl;

          // erase higher order moments
          for( int i = lowerBase.size(); i<numDofs; ++i )
          {
            llf[ i ] = 0;
          }
        }
      }
    };

    /** \brief Specialization of OrderReduction for TupleDiscreteFunction.
     *         The implementation default to the above for each tuple member.a
     *
     *  \note  This is only implemented for hierarchic spaces currently.
     */
    template < class... DFs >
    class OrderReduction< TupleDiscreteFunction< DFs... > >
    {
    public:
      typedef TupleDiscreteFunction< DFs... > DomainFunctionType;
      typedef DomainFunctionType  RangeFunctionType;

      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef DomainSpaceType RangeSpaceType;

      //static_assert(Dune::Fem::Capabilities::isHierarchic< DomainSpaceType > :: v, "OrderReduction is only implemented for hierarchic discrete function spaces");

      OrderReduction( const DomainSpaceType& )
      {
      }

      void operator () (const DomainFunctionType& arg, RangeFunctionType& dest )
      {
        Hybrid::forEach( typename DomainFunctionType::Sequence{}, [ & ]( auto i )
            {
              typedef typename DomainFunctionType :: template SubDiscreteFunction< i >::Type DF;
              OrderReduction< DF >()( arg.template subDiscreteFunction< i >(), dest.template subDiscreteFunction< i >() );
            } );
      }
    };

  } // end namespace Fem
} // end namespace Dune

#endif
