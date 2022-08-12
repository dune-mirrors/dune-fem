#ifndef DUNE_FEM_INTEGRAL_HH
#define DUNE_FEM_INTEGRAL_HH

#include <numeric>
#include <type_traits>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

#include <dune/fem/function/common/gridfunctionadapter.hh>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem/common/bindguard.hh>

namespace Dune
{

  namespace Fem
  {
    // IntegralBase
    // ----------

    template< class GridPart, class NormImplementation >
    class IntegralBase
      : public BartonNackmanInterface< IntegralBase< GridPart, NormImplementation >,
                                       NormImplementation >
    {
      typedef BartonNackmanInterface< IntegralBase< GridPart, NormImplementation >,
                                      NormImplementation >  BaseType ;
      typedef IntegralBase< GridPart, NormImplementation > ThisType;

    public:
      typedef GridPart GridPartType;

    protected:
      using BaseType :: asImp ;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      template <bool uDiscrete, bool vDiscrete>
      struct ForEachCaller
      {
      private:
        template <class IteratorRange, class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType>
        static ReturnType forEachLocal ( const ThisType &norm, const IteratorRange& iterators,
                                         const UDiscreteFunctionType &u, const VDiscreteFunctionType &v,
                                         const ReturnType &initialValue,
                                         unsigned int order )
        {
          static_assert( uDiscrete && vDiscrete, "Distance can only be calculated between GridFunctions" );

          ReturnType sum( 0 );
          {
            ConstLocalFunction< UDiscreteFunctionType > uLocal( u );
            ConstLocalFunction< VDiscreteFunctionType > vLocal( v );
            for( const EntityType &entity : iterators )
            {
              auto uGuard = bindGuard( uLocal, entity );
              auto vGuard = bindGuard( vLocal, entity );

              const unsigned int uOrder = uLocal.order();
              const unsigned int vOrder = vLocal.order();
              const unsigned int orderLocal = (order == 0 ? 2*std::max( uOrder, vOrder ) : order);
              norm.distanceLocal( entity, orderLocal, uLocal, vLocal, sum );
            }
          }
          return sum;
        }

      public:
        template <class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType, class PartitionSet>
        static ReturnType forEach ( const ThisType &norm, const UDiscreteFunctionType &u, const VDiscreteFunctionType &v,
                                    const ReturnType &initialValue,
                                    const PartitionSet& partitionSet,
                                    unsigned int order )
        {
          const int nThreads = MPIManager::numThreads();
          if( nThreads == 1 )
            return forEachLocal( norm, elements( norm.gridPart_, partitionSet ), u, v, initialValue, order );

          std::vector< ReturnType > sums( nThreads, ReturnType(0) );

          ThreadIterator< GridPartType, PartitionSet::partitionIterator() >
            iterators( norm.gridPart_ );

          auto doIntegrate = [ &norm, &iterators, &sums, &u, &v, &initialValue, &order ] ()
          {
            sums[ MPIManager::thread() ] = forEachLocal( norm, iterators, u, v, initialValue, order );
          };

          try {
            // run threaded
            MPIManager::run( doIntegrate );
          }
          catch ( const SingleThreadModeError& e )
          {
            // return single threaded variant in this case
            return forEachLocal( norm, elements( norm.gridPart_, partitionSet ), u, v, initialValue, order );
          }

          return std::accumulate( sums.begin(), sums.end(), ReturnType(0) );
        }
      };

      // this specialization creates a grid function adapter
      template <bool vDiscrete>
      struct ForEachCaller<false, vDiscrete>
      {
        template <class F,
                  class VDiscreteFunctionType,
                  class ReturnType,
                  class PartitionSet>
        static
        ReturnType forEach ( const ThisType& norm,
                             const F& f, const VDiscreteFunctionType& v,
                             const ReturnType& initialValue,
                             const PartitionSet& partitionSet,
                             const unsigned int order )
        {
          typedef GridFunctionAdapter< F, GridPartType>  GridFunction ;
          GridFunction u( "Integral::adapter" , f , v.space().gridPart(), v.space().order() );

          return ForEachCaller< true, vDiscrete > ::
            forEach( norm, u, v, initialValue, partitionSet, order );
        }
      };

      // this specialization simply switches arguments
      template <bool uDiscrete>
      struct ForEachCaller<uDiscrete, false>
      {
        template <class UDiscreteFunctionType,
                  class F,
                  class ReturnType,
                  class PartitionSet>
        static
        ReturnType forEach ( const ThisType& norm,
                             const UDiscreteFunctionType& u,
                             const F& f,
                             const ReturnType& initialValue,
                             const PartitionSet& partitionSet,
                             const unsigned int order )
        {
          return ForEachCaller< false, uDiscrete > ::
            forEach( norm, f, u, initialValue, partitionSet, order );
        }
      };

      template <class IteratorRange, class UDiscreteFunctionType, class ReturnType>
      ReturnType forEachLocal ( const IteratorRange& iterators,
                                const UDiscreteFunctionType &u,
                                const ReturnType &initialValue,
                                unsigned int order ) const
      {
        ReturnType sum( 0 );
        {
          ConstLocalFunction< UDiscreteFunctionType > uLocal( u );
          for( const EntityType &entity : iterators )
          {
            auto uGuard = bindGuard( uLocal, entity );

            const unsigned int uOrder = uLocal.order();
            const unsigned int orderLocal = (order == 0 ? 2*uOrder : order);
            normLocal( entity, orderLocal, uLocal, sum );
          }
        }
        return sum;
      }

      template< class DiscreteFunctionType, class ReturnType, class PartitionSet >
      ReturnType forEach ( const DiscreteFunctionType &u, const ReturnType &initialValue,
                           const PartitionSet& partitionSet,
                           unsigned int order = 0 ) const
      {
        static_assert( (std::is_base_of<Fem::HasLocalFunction, DiscreteFunctionType>::value),
                            "Norm only implemented for quantities implementing a local function!" );
        const int nThreads = MPIManager::numThreads();
        if( nThreads == 1 )
          return forEachLocal( elements( gridPart_, partitionSet ), u, initialValue, order );

        std::vector< ReturnType > sums( nThreads, ReturnType(0) );

        ThreadIterator< GridPartType, PartitionSet::partitionIterator() >
          iterators( gridPart_ );

        auto doIntegrate = [ this, &iterators, &sums, &u, &initialValue, &order ] ()
        {
          sums[ MPIManager::thread() ] = forEachLocal( iterators, u, initialValue, order );
        };

        try {
          // run threaded
          MPIManager::run( doIntegrate );
        }
        catch ( const SingleThreadModeError& e )
        {
          // return single threaded variant in this case
          return forEachLocal( elements( gridPart_, partitionSet ), u, initialValue, order );
        }

        return std::accumulate( sums.begin(), sums.end(), ReturnType(0) );
      }

      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType, class PartitionSet >
      ReturnType forEach ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v,
                           const ReturnType &initialValue, const PartitionSet& partitionSet,
                           unsigned int order = 0 ) const
      {
        enum { uDiscrete = std::is_convertible<UDiscreteFunctionType, HasLocalFunction>::value };
        enum { vDiscrete = std::is_convertible<VDiscreteFunctionType, HasLocalFunction>::value };

        // call forEach depending on which argument is a grid function,
        // i.e. has a local function
        return ForEachCaller< uDiscrete, vDiscrete > ::
                  forEach( *this, u, v, initialValue, partitionSet, order );
      }

    public:
      explicit IntegralBase ( const GridPartType &gridPart )
        : gridPart_( gridPart )
      {}

    protected:
      template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().distanceLocal( entity, order, uLocal, vLocal, sum ) );
      }

      template< class LocalFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().normLocal( entity, order, uLocal, sum ) );
      }

      const GridPartType &gridPart () const { return gridPart_; }

      const typename GridPartType::CommunicationType& comm () const
      {
        return gridPart().comm();
      }

      bool checkCommunicateFlag( bool communicate ) const
      {
#ifndef NDEBUG
        bool commMax = communicate;
        assert( communicate == comm().max( commMax ) );
#endif
        return communicate;
      }

    private:
      const GridPartType &gridPart_;

    };

    // Integral
    // --------

    template< class GridPart >
    class Integral : public IntegralBase< GridPart, Integral< GridPart > >
    {
      typedef IntegralBase< GridPart, Integral< GridPart > > BaseType ;
      typedef Integral< GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      using BaseType :: gridPart ;
      using BaseType :: comm ;

    protected:

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
      explicit Integral ( const GridPartType &gridPart,
                          const unsigned int order = 0,
                          const bool communicate = true );


      //! || u ||_L1 on given set of entities (partition set)
      template< class DiscreteFunctionType, class PartitionSet >
      typename DiscreteFunctionType::RangeType
      norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const;

      //! || u ||_L1 on interior partition entities
      template< class DiscreteFunctionType >
      typename DiscreteFunctionType::RangeType
      norm ( const DiscreteFunctionType &u ) const
      {
        return norm( u, Partitions::interior );
      }

      //! || u - v ||_L2 on given set of entities (partition set)
      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
      typename UDiscreteFunctionType::RangeType
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet& partitionSet ) const;

      //! || u - v ||_L2 on interior partition entities
      template< class UDiscreteFunctionType, class VDiscreteFunctionType >
      typename UDiscreteFunctionType::RangeType
      distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
      {
        return distance( u, v, Partitions::interior );
      }

      template< class LocalFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const;

      template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const;
    };



    // Implementation of Integral
    // ------------------------

    template< class GridPart >
    inline Integral< GridPart >::Integral ( const GridPartType &gridPart, const unsigned int order, const bool communicate )
    : BaseType( gridPart ),
      order_( order ),
      communicate_( BaseType::checkCommunicateFlag( communicate ) )
    {}


    template< class GridPart >
    template< class DiscreteFunctionType, class PartitionSet >
    inline typename DiscreteFunctionType::RangeType
    Integral< GridPart >::norm ( const DiscreteFunctionType &u, const PartitionSet& partitionSet ) const
    {
      typedef typename DiscreteFunctionType::RangeType RangeType;
      typedef RangeType ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum = comm().sum( sum );
      }

      return sum;
    }


    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
    inline typename UDiscreteFunctionType::RangeType
    Integral< GridPart >
      ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet& partitionSet ) const
    {
      typedef typename UDiscreteFunctionType::RangeType RangeType;
      typedef RangeType ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0), partitionSet, order_ );

      // communicate_ indicates global norm
      if( communicate_ )
      {
        sum = comm().sum( sum );
      }

      return sum;
    }

    template< class GridPart >
    template< class LocalFunctionType, class ReturnType >
    inline void
    Integral< GridPart >::normLocal ( const EntityType &entity, unsigned int order, const LocalFunctionType &uLocal, ReturnType &sum ) const
    {
      Integrator< QuadratureType > integrator( order );

      integrator.integrateAdd( entity, uLocal, sum );
    }

    template< class GridPart >
    template< class ULocalFunctionType, class VLocalFunctionType, class ReturnType >
    inline void
    Integral< GridPart >::distanceLocal ( const EntityType &entity, unsigned int order, const ULocalFunctionType &uLocal, const VLocalFunctionType &vLocal, ReturnType &sum ) const
    {
      Integrator< QuadratureType > integrator( order );

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      LocalDistanceType dist( uLocal, vLocal );

      integrator.integrateAdd( entity, dist, sum );
    }

    template< class GridPart >
    template< class UFunction, class VFunction >
    struct Integral< GridPart >::FunctionDistance
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

#endif // #ifndef DUNE_FEM_INTEGRAL_HH
