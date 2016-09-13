#ifndef DUNE_FEM_H1NORM_HH
#define DUNE_FEM_H1NORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart >
    class H1Norm
    : public IntegralBase< GridPart, H1Norm< GridPart> >
    {
      typedef H1Norm< GridPart > ThisType;
      typedef IntegralBase< GridPart, H1Norm< GridPart> > BaseType;

    public:
      typedef GridPart GridPartType;

      template< class Function >
      struct FunctionJacobianSquare;

    protected:
      typedef typename BaseType::EntityType EntityType;

      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    public:
      typedef Integrator< QuadratureType > IntegratorType;


      using BaseType::gridPart;
      using BaseType::comm;

      /** \brief constructor
       *    \param gridPart     specific gridPart for selection of entities
       *    \param order        order of integration quadrature (default = 2*space.order())
       *    \param communicate  if true global (over all ranks) norm is computed (default = true)
       */
      explicit H1Norm ( const GridPartType &gridPart,
                        const unsigned int order = 0,
                        const bool communicate = true );

      //! || u ||_H1 on given set of entities (partition set)
      template< class DiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const;

      //! || u ||_H1 on interior partition entities
      template< class DiscreteFunctionType >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u ) const
      {
        return norm( u, Partitions::interior );
      }

      //! || u - v ||_H1 on given set of entities (partition set)
      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet& partitionSet ) const;

      //! || u - v ||_H1 on interior partition entities
      template< class UDiscreteFunctionType, class VDiscreteFunctionType >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
      {
        return distance( u, v, Partitions::interior );
      }

      template< class LocalFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const;

      template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const;

      ThisType operator= ( const ThisType& ) = delete;

    protected:
      const unsigned int order_;
      const bool communicate_;
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
    inline H1Norm< GridPart >::H1Norm ( const GridPartType &gridPart, const unsigned int order, const bool communicate )
    : BaseType( gridPart ),
      order_( order ),
      communicate_( communicate )
    {}


    template< class GridPart >
    template< class DiscreteFunctionType, class PartitionSet >
    inline typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    H1Norm< GridPart >::norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      ReturnType sum = BaseType :: forEach( u, ReturnType( 0 ), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum = comm().sum( sum[ 0 ] );
      }

      // return result, e.g. sqrt of calculated sum
      return sqrt( sum[ 0 ] );
    }

    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
    inline typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    H1Norm< GridPart >::distance ( const UDiscreteFunctionType &u,
                                   const VDiscreteFunctionType &v,
                                   const PartitionSet& partitionSet ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      ReturnType sum = BaseType :: forEach( u, v, ReturnType( 0 ), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum = comm().sum( sum[ 0 ] );
      }

      // return result, e.g. sqrt of calculated sum
      return sqrt( sum[ 0 ] );
    }

    template< class GridPart >
    template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
    inline void
    H1Norm< GridPart >::distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const
    {
      typedef typename L2Norm< GridPart >::template FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      IntegratorType integrator( order );

      LocalDistanceType dist( uLocal, vLocal );
      FunctionJacobianSquare< LocalDistanceType > dist2( dist );

      integrator.integrateAdd( entity, dist2, sum );
    }


    template< class GridPart >
    template< class LocalFunctionType, class ReturnType >
    inline void
    H1Norm< GridPart >::normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
    {
      // evaluate norm locally
      IntegratorType integrator( order );

      FunctionJacobianSquare< LocalFunctionType > uLocal2( uLocal );

      integrator.integrateAdd( entity, uLocal2, sum );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_H1NORM_HH
