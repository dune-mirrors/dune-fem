#ifndef DUNE_FEM_L1NORM_HH
#define DUNE_FEM_L1NORM_HH

#include <dune/fem/quadrature/integrator.hh>

#include <dune/fem/misc/lpnorm.hh>

namespace Dune
{

  namespace Fem
  {

    // L1Norm
    // ------

    template< class GridPart >
    class L1Norm : public LPNormBase< GridPart, L1Norm< GridPart > >
    {
      typedef LPNormBase< GridPart, L1Norm< GridPart > > BaseType ;
      typedef L1Norm< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      using BaseType :: gridPart ;
      using BaseType :: comm ;

    protected:
      template< class Function >
      struct FunctionAbs;

      template< class UFunction, class VFunction >
      struct FunctionDistance;

      typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
      typedef typename GridIteratorType::Entity EntityType;
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

      const unsigned int order_;
    public:
      explicit L1Norm ( const GridPartType &gridPart, const unsigned int order = 0 );

      template< class DiscreteFunctionType >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u ) const;

      template< class UDiscreteFunctionType, class VDiscreteFunctionType >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

      template< class UDiscreteFunctionType,
                class VDiscreteFunctionType,
                class ReturnType >
      inline void
      distanceLocal ( const EntityType& entity, const unsigned int order,
                      const UDiscreteFunctionType &u,
                      const VDiscreteFunctionType &v,
                      ReturnType& sum ) const ;

      template< class UDiscreteFunctionType,
                class ReturnType >
      inline void
      normLocal ( const EntityType& entity, const unsigned int order,
                      const UDiscreteFunctionType &u,
                      ReturnType& sum ) const ;
    };



    // Implementation of L1Norm
    // ------------------------

    template< class GridPart >
    inline L1Norm< GridPart >::L1Norm ( const GridPartType &gridPart, const unsigned int order )
    : BaseType( gridPart ),
      order_( order )
    {}


    template< class GridPart >
    template< class DiscreteFunctionType >
    inline typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    L1Norm< GridPart >::norm ( const DiscreteFunctionType &u ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0), order_ );

      // return result, e.g. sum
      return comm().sum( sum[ 0 ] );
    }


    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    L1Norm< GridPart >
      ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0), order_ );

      // return result, e.g. sum
      return comm().sum( sum[ 0 ] );
    }

    template< class GridPart >
    template< class DiscreteFunctionType, class ReturnType >
    inline void
    L1Norm< GridPart >::normLocal ( const EntityType& entity, const unsigned int order,
                                    const DiscreteFunctionType &u,
                                    ReturnType& sum ) const
    {
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      Integrator< QuadratureType > integrator( order );

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionAbs< LocalFunctionType > ulocalAbs( ulocal );

      integrator.integrateAdd( entity, ulocalAbs, sum );
    }

    template< class GridPart >
    template< class UDiscreteFunctionType,
              class VDiscreteFunctionType,
              class ReturnType >
    inline void
    L1Norm< GridPart >::distanceLocal ( const EntityType& entity, const unsigned int order,
                                        const UDiscreteFunctionType &u,
                                        const VDiscreteFunctionType &v,
                                        ReturnType& sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
      typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

      Integrator< QuadratureType > integrator( order );

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
        LocalDistanceType;

      LocalDistanceType dist( ulocal, vlocal );
      FunctionAbs< LocalDistanceType > distAbs( dist );

      integrator.integrateAdd( entity, distAbs, sum );
    }


    template< class GridPart >
    template< class Function >
    struct L1Norm< GridPart >::FunctionAbs
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > RangeType;

      explicit FunctionAbs ( const FunctionType &function )
      : function_( function )
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret = phi.one_norm();
      }

    private:
      const FunctionType &function_;
    };


    template< class GridPart >
    template< class UFunction, class VFunction >
    struct L1Norm< GridPart >::FunctionDistance
    {
      typedef UFunction UFunctionType;
      typedef VFunction VFunctionType;

      typedef typename UFunctionType::RangeFieldType RangeFieldType;
      typedef typename UFunctionType::RangeType RangeType;
      typedef typename UFunctionType::JacobianRangeType JacobianRangeType;

      FunctionDistance ( const UFunctionType &u, const VFunctionType &v )
      : u_( u ), v_( v )
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        RangeType phi;
        u_.evaluate( x, ret );
        v_.evaluate( x, phi );
        ret -= phi;
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &ret ) const
      {
        JacobianRangeType phi;
        u_.jacobian( x, ret );
        v_.jacobian( x, phi );
        ret -= phi;
      }

    private:
      const UFunctionType &u_;
      const VFunctionType &v_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_L1NORM_HH
