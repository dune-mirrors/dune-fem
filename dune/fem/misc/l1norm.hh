#ifndef DUNE_FEM_L1NORM_HH
#define DUNE_FEM_L1NORM_HH

#include <dune/fem/quadrature/integrator.hh>

#include <dune/fem/misc/domainintegral.hh>

namespace Dune
{

  namespace Fem
  {

    // L1Norm
    // ------

    template< class GridPart >
    class L1Norm : public IntegralBase< GridPart, L1Norm< GridPart > >
    {
      typedef IntegralBase< GridPart, L1Norm< GridPart > > BaseType ;
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

      typedef typename BaseType::EntityType EntityType;
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

      const unsigned int order_;
      const bool communicate_;
    public:
      /** \brief constructor
       *    \param gridPart     specific gridPart for selection of entities
       *    \param order        order of integration quadrature (default = 2*space.order())
       *    \param communicate  if true global (over all ranks) norm is computed (default = true)
       */
      explicit L1Norm ( const GridPartType &gridPart,
                        const unsigned int order = 0,
                        const bool communicate = true );


      //! || u ||_L1 on given set of entities (partition set)
      template< class DiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const;

      //! || u ||_L1 on interior partition entities
      template< class DiscreteFunctionType >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u ) const
      {
        return norm( u, Partitions::interior );
      }

      //! || u - v ||_L2 on given set of entities (partition set)
      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet& partitionSet ) const;

      //! || u - v ||_L2 on interior partition entities
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
    };



    // Implementation of L1Norm
    // ------------------------

    template< class GridPart >
    inline L1Norm< GridPart >::L1Norm ( const GridPartType &gridPart, const unsigned int order, const bool communicate )
    : BaseType( gridPart ),
      order_( order ),
      communicate_( BaseType::checkCommunicateFlag( communicate ) )
    {}


    template< class GridPart >
    template< class DiscreteFunctionType, class PartitionSet >
    inline typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    L1Norm< GridPart >::norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum[ 0 ] = comm().sum( sum[ 0 ] );
      }

      return sum[ 0 ];
    }


    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
    inline typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    L1Norm< GridPart >
      ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet& partitionSet ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum[ 0 ] = comm().sum( sum[ 0 ] );
      }

      return sum[ 0 ];
    }

    template< class GridPart >
    template< class LocalFunctionType, class ReturnType >
    inline void
    L1Norm< GridPart >::normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
    {
      Integrator< QuadratureType > integrator( order );

      FunctionAbs< LocalFunctionType > uLocalAbs( uLocal );

      integrator.integrateAdd( entity, uLocalAbs, sum );
    }

    template< class GridPart >
    template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
    inline void
    L1Norm< GridPart >::distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const
    {
      Integrator< QuadratureType > integrator( order );

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      LocalDistanceType dist( uLocal, vLocal );
      FunctionAbs< LocalDistanceType > distAbs( dist );

      integrator.integrateAdd( entity, distAbs, sum );
    }


    template< class GridPart >
    template< class Function >
    struct L1Norm< GridPart >::FunctionAbs
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > RangeType;

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
