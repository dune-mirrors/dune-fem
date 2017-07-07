#ifndef DUNE_FEM_L2NORM_HH
#define DUNE_FEM_L2NORM_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

#include <dune/fem/misc/domainintegral.hh>

namespace Dune
{

  namespace Fem
  {

    // L2Norm
    // ------

    template< class GridPart >
    class L2Norm : public IntegralBase< GridPart, L2Norm< GridPart > >
    {
      typedef IntegralBase< GridPart, L2Norm< GridPart > > BaseType ;
      typedef L2Norm< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      using BaseType :: gridPart ;
      using BaseType :: comm ;

      template< class Function >
      struct FunctionSquare;

      template< class UFunction, class VFunction >
      struct FunctionDistance;

    protected:
      typedef typename BaseType::EntityType EntityType;
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
      typedef Integrator< QuadratureType > IntegratorType;

      const unsigned int order_;
      const bool communicate_;
    public:
      /** \brief constructor
       *    \param gridPart     specific gridPart for selection of entities
       *    \param order        order of integration quadrature (default = 2*space.order())
       *    \param communicate  if true global (over all ranks) norm is computed (default = true)
       */
      explicit L2Norm ( const GridPartType &gridPart,
                        const unsigned int order = 0,
                        const bool communicate = true );

      //! || u ||_L2 on given set of entities (partition set)
      template< class DiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const;

      //! || u ||_L2 on interior partition entities
      template< class DiscreteFunctionType >
      typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
      norm ( const DiscreteFunctionType &u ) const
      {
        return norm( u, Partitions::interior );
      }

      //! || u - v ||_L2 on given set of entities (partition set)
      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u,
                 const VDiscreteFunctionType &v,
                 const PartitionSet& partitionSet ) const;

      //! || u - v ||_L2 on interior partition entities
      template< class UDiscreteFunctionType, class VDiscreteFunctionType >
      typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
      {
        return distance( u, v, Partitions::interior );
      }

      template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const;

      template< class LocalFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const;
    };


    // Implementation of L2Norm
    // ------------------------

    template< class GridPart >
    inline L2Norm< GridPart >::L2Norm ( const GridPartType &gridPart, const unsigned int order, const bool communicate )
    : BaseType( gridPart ),
      order_( order ),
      communicate_( BaseType::checkCommunicateFlag( communicate ) )
    {
    }


    template< class GridPart >
    template< class DiscreteFunctionType, class PartitionSet >
    typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    L2Norm< GridPart >::norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const
    {
      typedef typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum[ 0 ] = comm().sum( sum[ 0 ] );
      }

      // return result, e.g. sqrt of calculated sum
      return std::sqrt( sum[ 0 ] );
    }


    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
    typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    L2Norm< GridPart >
      ::distance ( const UDiscreteFunctionType &u,
                   const VDiscreteFunctionType &v,
                   const PartitionSet& partitionSet ) const
    {
      typedef typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum[ 0 ] = comm().sum( sum[ 0 ] );
      }

      // return result, e.g. sqrt of calculated sum
      return std::sqrt( sum[ 0 ] );
    }

    template< class GridPart >
    template< class LocalFunctionType, class ReturnType >
    inline void L2Norm< GridPart >::normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
    {
      // evaluate norm locally
      IntegratorType integrator( order );

      FunctionSquare< LocalFunctionType > uLocal2( uLocal );

      integrator.integrateAdd( entity, uLocal2, sum );
    }

    template< class GridPart >
    template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
    inline void
    L2Norm< GridPart >::distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const
    {
      // evaluate norm locally
      IntegratorType integrator( order );

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      LocalDistanceType dist( uLocal, vLocal );
      FunctionSquare< LocalDistanceType > dist2( dist );

      integrator.integrateAdd( entity, dist2, sum );
    }


    template< class GridPart >
    template< class Function >
    struct L2Norm< GridPart >::FunctionSquare
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > RangeType;

      explicit FunctionSquare ( const FunctionType &function )
      : function_( function )
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret = phi.two_norm2();
      }

    private:
      const FunctionType &function_;
    };


    template< class GridPart >
    template< class UFunction, class VFunction >
    struct L2Norm< GridPart >::FunctionDistance
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



    // WeightedL2Norm
    // --------------

    template< class WeightFunction >
    class WeightedL2Norm
    : public L2Norm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType >
    {
      typedef WeightedL2Norm< WeightFunction > ThisType;
      typedef L2Norm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType > BaseType;

    public:
      typedef WeightFunction WeightFunctionType;

      typedef typename WeightFunctionType::DiscreteFunctionSpaceType WeightFunctionSpaceType;
      typedef typename WeightFunctionSpaceType::GridPartType GridPartType;

    protected:
      template< class Function >
      struct WeightedFunctionSquare;

      typedef ConstLocalFunction< WeightFunctionType > LocalWeightFunctionType;
      typedef typename WeightFunctionType::RangeType WeightType;

      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntegratorType IntegratorType;

      using BaseType::gridPart;
      using BaseType::comm;

    public:

      using BaseType::norm;
      using BaseType::distance;

      explicit WeightedL2Norm ( const WeightFunctionType &weightFunction, const unsigned int order = 0, const bool communicate = true );

      template< class LocalFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const;

      template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const;

    private:
      LocalWeightFunctionType wfLocal_;
    };




    // Implementation of WeightedL2Norm
    // --------------------------------

    template< class WeightFunction >
    inline WeightedL2Norm< WeightFunction >
      ::WeightedL2Norm ( const WeightFunctionType &weightFunction, const unsigned int order, const bool communicate )
    : BaseType( weightFunction.space().gridPart(), order, communicate ),
      wfLocal_( weightFunction )
    {
      static_assert( (WeightFunctionSpaceType::dimRange == 1),
                     "Wight function must be scalar." );
    }


    template< class WeightFunction >
    template< class LocalFunctionType, class ReturnType >
    inline void
    WeightedL2Norm< WeightFunction >::normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
    {
      // !!!! order !!!!
      IntegratorType integrator( order );

      wfLocal_.bind( entity );
      WeightedFunctionSquare< LocalFunctionType > uLocal2( wfLocal_, uLocal );
      integrator.integrateAdd( entity, uLocal2, sum );
      wfLocal_.unbind();
    }


    template< class WeightFunction >
    template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
    inline void
    WeightedL2Norm< WeightFunction >::distanceLocal ( const EntityType& entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType& sum ) const
    {
      typedef typename BaseType::template FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      // !!!! order !!!!
      IntegratorType integrator( order );

      wfLocal_.bind( entity );
      LocalDistanceType dist( uLocal, vLocal );
      WeightedFunctionSquare< LocalDistanceType > dist2( wfLocal_, dist );

      integrator.integrateAdd( entity, dist2, sum );
      wfLocal_.unbind();
    }


    template< class WeightFunction >
    template< class Function >
    struct WeightedL2Norm< WeightFunction >::WeightedFunctionSquare
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
      typedef FieldVector< RealType, 1 > RangeType;

      WeightedFunctionSquare ( const LocalWeightFunctionType &weightFunction,
                               const FunctionType &function )
      : weightFunction_( weightFunction ),
        function_( function )
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        WeightType weight;
        weightFunction_.evaluate( x, weight );

        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret = weight[ 0 ] * (phi * phi);
      }

    private:
      const LocalWeightFunctionType &weightFunction_;
      const FunctionType &function_;
    };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_L2NORM_HH
