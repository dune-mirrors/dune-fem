#ifndef DUNE_FEM_H1NORM_HH
#define DUNE_FEM_H1NORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart >
    class H1Norm
    : public LPNormBase< GridPart, H1Norm< GridPart> >
    {
      typedef H1Norm< GridPart > ThisType;
      typedef LPNormBase< GridPart, H1Norm< GridPart> > BaseType;

    public:
      typedef GridPart GridPartType;

      template< class Function >
      struct FunctionJacobianSquare;

    protected:
      typedef typename BaseType::GridIteratorType GridIteratorType;

      typedef typename GridIteratorType::Entity EntityType;

      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    public:
      typedef Integrator< QuadratureType > IntegratorType;


      using BaseType::gridPart;
      using BaseType::comm;

    public:
      explicit H1Norm ( const GridPartType &gridPart, const unsigned int order = 0 );
      H1Norm ( const ThisType &other );

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

    private:
      // prohibit assignment
      ThisType operator= ( const ThisType &other );

      const unsigned int order_;
    };



    // H1Norm::FunctionJacobianSquare
    // ------------------------------

    template< class GridPart >
    template< class Function >
    struct H1Norm< GridPart >::FunctionJacobianSquare
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< typename FunctionType::RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > RangeType;

    public:
      explicit FunctionJacobianSquare ( const FunctionType &function )
      : function_( function )
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        const int dimRange = FunctionType::RangeType::dimension;

        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret[ 0 ] = phi.two_norm2();

        typename FunctionType::JacobianRangeType grad;
        function_.jacobian( x, grad );
        for( int i = 0; i < dimRange; ++i )
          ret[ 0 ] += grad[ i ].two_norm2();
      }

    private:
      const FunctionType &function_;
    };



    // Implementation of H1 Norm
    // -------------------------

    template< class GridPart >
    inline H1Norm< GridPart >::H1Norm ( const GridPartType &gridPart, const unsigned int order )
    : BaseType( gridPart ),
      order_( order )
    {}



    template< class GridPart >
    inline H1Norm< GridPart >::H1Norm ( const ThisType &other )
    : BaseType( other ),
      order_( other.order_ )
    {}


    template< class GridPart >
    template< class DiscreteFunctionType >
    inline typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    H1Norm< GridPart >::norm ( const DiscreteFunctionType &u ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      ReturnType sum = BaseType :: forEach( u, ReturnType( 0 ), order_ );

      // return result, e.g. sqrt of calculated sum
      return sqrt( comm().sum( sum[ 0 ] ) );
    }

    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    H1Norm< GridPart >::distance ( const UDiscreteFunctionType &u,
                                   const VDiscreteFunctionType &v ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      ReturnType sum = BaseType :: forEach( u, v, ReturnType( 0 ), order_ );

      // return result, e.g. sqrt of calculated sum
      return sqrt( comm().sum( sum[ 0 ] ) );
    }

    template< class GridPart >
    template< class UDiscreteFunctionType,
              class VDiscreteFunctionType,
              class ReturnType >
    inline void
    H1Norm< GridPart >::distanceLocal ( const EntityType& entity, const unsigned int order,
                                        const UDiscreteFunctionType &u,
                                        const VDiscreteFunctionType &v,
                                        ReturnType& sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
      typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;
      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      typedef typename L2Norm< GridPart >::template FunctionDistance< ULocalFunctionType, VLocalFunctionType >
        LocalDistanceType;

      IntegratorType integrator( order );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionJacobianSquare< LocalDistanceType > dist2( dist );

      integrator.integrateAdd( entity, dist2, sum );
    }


    template< class GridPart >
    template< class DiscreteFunctionType, class ReturnType >
    inline void
    H1Norm< GridPart >::normLocal ( const EntityType& entity, const unsigned int order,
                                    const DiscreteFunctionType &u,
                                    ReturnType& sum ) const
    {
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      // evaluate norm locally

      IntegratorType integrator( order );

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionJacobianSquare< LocalFunctionType > ulocal2( ulocal );

      integrator.integrateAdd( entity, ulocal2, sum );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_H1NORM_HH
