#ifndef DUNE_FEM_LPNORM_HH
#define DUNE_FEM_LPNORM_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

#include <dune/fem/function/common/gridfunctionadapter.hh>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
  {

    // LPNormBase
    // ----------

    template< class GridPart, class NormImplementation >
    class LPNormBase
      : public BartonNackmanInterface< LPNormBase< GridPart, NormImplementation >,
                                       NormImplementation >
    {
      typedef BartonNackmanInterface< LPNormBase< GridPart, NormImplementation >,
                                      NormImplementation >  BaseType ;
      typedef LPNormBase< GridPart, NormImplementation > ThisType;

    public:
      typedef GridPart GridPartType;

    protected:
      using BaseType :: asImp ;

      typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
      typedef typename GridIteratorType::Entity EntityType;

      template < class ReturnType,
                 class UDiscreteFunctionType,
                 class VDiscreteFunctionType >
      struct NormOnEntityFunctor
      {
        const ThisType& norm_;
        const UDiscreteFunctionType& u_;
        const VDiscreteFunctionType* v_;
        const unsigned int order_ ;
        ReturnType sum_;

        unsigned int getOrder( const UDiscreteFunctionType& u,
                               const VDiscreteFunctionType& v ) const
        {
          const unsigned int uorder = u.space().order();
          const unsigned int vorder = v.space().order();
          const unsigned int order = 2 * std::max( uorder, vorder );
          return order;
        }

        NormOnEntityFunctor( const ThisType& norm,
                      const UDiscreteFunctionType& u,
                      const unsigned int order = 0 )
          : norm_( norm ),
            u_( u ),
            v_( 0 ),
            order_( ( order == 0 ) ? ( getOrder( u, u ) ) : order ),
            sum_( 0 )
        {}

        NormOnEntityFunctor( const ThisType& norm,
                      const UDiscreteFunctionType& u,
                      const VDiscreteFunctionType& v,
                      const unsigned int order = 0 )
          : norm_( norm ),
            u_( u ),
            v_( &v ),
            order_( ( order == 0 ) ? ( getOrder( u, v ) ) : order ),
            sum_( 0 )
        {
        }

        // apply norm operation on entity
        void operator() ( const EntityType& entity )
        {
          // if v is not zero apply distance
          if( v_ )
            norm_.distanceLocal( entity, order_, u_, *v_, sum_ );
          else
            norm_.normLocal( entity, order_, u_, sum_ );
        }

        const ReturnType& result() const { return sum_; }
      };

      template <bool uDiscrete, bool vDiscrete>
      struct ForEachCaller
      {
        template <class UDiscreteFunctionType,
                  class VDiscreteFunctionType,
                  class ReturnType>
        static
        ReturnType forEach ( const ThisType& norm,
                             const UDiscreteFunctionType& u,
                             const VDiscreteFunctionType& v,
                             const ReturnType& initialValue,
                             const unsigned int order )
        {
          static_assert( uDiscrete && vDiscrete, "Distance can only be calculated between GridFunctions" );

          typedef NormOnEntityFunctor< ReturnType, UDiscreteFunctionType, VDiscreteFunctionType >
            ForEachFunctorType;

          // create functor
          ForEachFunctorType normOnEntity( norm, u, v, order );

          // do grid part traversal and apply action on entity
          u.space().forEach( normOnEntity );

          return normOnEntity.result();
        }
      };

      // this specialization creates a grid function adapter
      template <bool vDiscrete>
      struct ForEachCaller<false, vDiscrete>
      {
        template <class F,
                  class VDiscreteFunctionType,
                  class ReturnType>
        static
        ReturnType forEach ( const ThisType& norm,
                             const F& f, const VDiscreteFunctionType& v,
                             const ReturnType& initialValue,
                             const unsigned int order )
        {
          typedef GridFunctionAdapter< F, GridPartType>  GridFunction ;
          GridFunction u( "LPNorm::adapter" , f , v.space().gridPart(), v.space().order() );

          return ForEachCaller< true, vDiscrete > ::
            forEach( norm, u, v, initialValue, order );
        }
      };

      // this specialization simply switches arguments
      template <bool uDiscrete>
      struct ForEachCaller<uDiscrete, false>
      {
        template <class UDiscreteFunctionType,
                  class F,
                  class ReturnType>
        static
        ReturnType forEach ( const ThisType& norm,
                             const UDiscreteFunctionType& u,
                             const F& f,
                             const ReturnType& initialValue,
                             const unsigned int order )
        {
          return ForEachCaller< false, uDiscrete > ::
            forEach( norm, f, u, initialValue, order );
        }
      };

      template <class DiscreteFunctionType, class ReturnType>
      ReturnType
      forEach ( const DiscreteFunctionType& u,
                const ReturnType& initialValue,
                const unsigned int order = 0 ) const
      {
        static_assert( (IsBaseOf<Fem::HasLocalFunction, DiscreteFunctionType>::value),
                            "Norm only implemented for quantities implementing a local function!" );

        typedef NormOnEntityFunctor< ReturnType, DiscreteFunctionType, DiscreteFunctionType >
          ForEachFunctorType;

        // create functor
        ForEachFunctorType normOnEntity( *this, u, order );

        // do grid part traversal and apply action on entity
        u.space().forEach( normOnEntity );

        return normOnEntity.result();
      }

      template <class UDiscreteFunctionType,
                class VDiscreteFunctionType,
                class ReturnType>
      ReturnType
      forEach ( const UDiscreteFunctionType& u,
                const VDiscreteFunctionType& v,
                const ReturnType& initialValue,
                const unsigned int order = 0 ) const
      {
        enum { uDiscrete = Conversion<UDiscreteFunctionType, HasLocalFunction>::exists };
        enum { vDiscrete = Conversion<VDiscreteFunctionType, HasLocalFunction>::exists };

        // call forEach depending on which argument is a grid function,
        // i.e. has a local function
        return ForEachCaller< uDiscrete, vDiscrete > ::
                  forEach( *this, u, v, initialValue, order );
      }

    public:
      explicit LPNormBase ( const GridPartType &gridPart )
        : gridPart_( gridPart )
      {}

    protected:
      template< class UDiscreteFunctionType,
                class VDiscreteFunctionType,
                class ReturnType >
      inline void
      distanceLocal ( const EntityType& entity, const unsigned int order,
                      const UDiscreteFunctionType &u,
                      const VDiscreteFunctionType &v,
                      ReturnType& sum ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().distanceLocal( entity, order, u, v, sum ) );
      }

      template< class UDiscreteFunctionType,
                class ReturnType >
      inline void
      normLocal ( const EntityType& entity, const unsigned int order,
                  const UDiscreteFunctionType &u,
                  ReturnType& sum ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().normLocal( entity, order, u, sum ) );
      }

      const GridPartType &gridPart () const { return gridPart_; }

      typename GridPartType::CollectiveCommunicationType comm () const
      {
        return gridPart().comm();
      }

    private:
      const GridPartType &gridPart_;
    };



    // TODO weighte LP norm might be adapted later
    // LPNorm
    //
    // !!!!! It is not cleared which quadrature order have to be applied for p > 2!!!!
    // !!!!! For p = 1 this norm does not work !!!!
    // ------


    //! Quadrature Order Interface
    struct OrderCalculatorInterface
    {
      virtual int operator() (const double p)=0;
    };

    //! default Quadrature Order Calculator
    //  can be re-implemented in order to use
    //  a different type of calculation
    //  which can be sepcified in the second template argument of LPNorm
    struct DefaultOrderCalculator : public OrderCalculatorInterface
    {
      int operator() (const double p)
      {
        int ret=0;
        const double q = p / (p-1);
        double max = std::max(p,q);
        assert(max < std::numeric_limits<int>::max()/2. );
        ret = max +1;
        return ret;
      }
    };

    template< class GridPart, class OrderCalculator = DefaultOrderCalculator >
    class LPNorm : public LPNormBase < GridPart, LPNorm< GridPart, OrderCalculator > >
    {
      typedef LPNorm< GridPart, OrderCalculator > ThisType;
      typedef LPNormBase< GridPart, ThisType > BaseType ;

    public:
      typedef GridPart GridPartType;

      using BaseType::gridPart;
      using BaseType::comm;

    protected:
      template< class Function >
      struct FunctionMultiplicator;

      template< class UFunction, class VFunction >
      struct FunctionDistance;

      typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
      typedef typename GridIteratorType::Entity EntityType;
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
      typedef Integrator< QuadratureType > IntegratorType;

    public:
      explicit LPNorm ( const GridPartType &gridPart, const double p );

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

      int order ( const int spaceOrder ) const ;

    protected:
      double p_ ;
      OrderCalculator *orderCalc_;
    };




    // WeightedLPNorm
    // --------------

    template< class WeightFunction, class OrderCalculator = DefaultOrderCalculator >
    class WeightedLPNorm
    : public LPNorm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType,
                     OrderCalculator >
    {
      typedef WeightedLPNorm< WeightFunction, OrderCalculator > ThisType;
      typedef LPNorm< typename WeightFunction::DiscreteFunctionSpaceType::GridPartType, OrderCalculator> BaseType;

    public:
      typedef WeightFunction WeightFunctionType;

      typedef typename WeightFunctionType::DiscreteFunctionSpaceType WeightFunctionSpaceType;
      typedef typename WeightFunctionSpaceType::GridPartType GridPartType;

    protected:
      template< class Function >
      struct WeightedFunctionMultiplicator;

      typedef typename WeightFunctionType::LocalFunctionType LocalWeightFunctionType;
      typedef typename WeightFunctionType::RangeType WeightType;

      typedef typename BaseType::GridIteratorType GridIteratorType;
      typedef typename BaseType::IntegratorType IntegratorType;

      typedef typename GridIteratorType::Entity EntityType;

      using BaseType::gridPart;
      using BaseType::comm;

    public:
      using BaseType::norm;
      using BaseType::distance;

      explicit WeightedLPNorm ( const WeightFunctionType &weightFunction, const double p );

      template< class UDiscreteFunctionType, class ReturnType >
      void normLocal ( const EntityType &entity, const int order,
              const UDiscreteFunctionType &u,
              ReturnType& sum ) const;

      template< class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType >
      void distanceLocal ( const EntityType &entity, const int order,
              const UDiscreteFunctionType &u,
              const VDiscreteFunctionType &v,
              ReturnType& sum ) const;

    private:
      const WeightFunctionType &weightFunction_;
      const double p_;
    };


    // Implementation of LPNorm
    // ------------------------

    template< class GridPart, class OrderCalculator >
    inline LPNorm< GridPart, OrderCalculator >::LPNorm ( const GridPartType &gridPart, const double p )
      :  BaseType( gridPart ),
         p_(p),
         orderCalc_( new OrderCalculator() )
    {
    }

    template< class GridPart, class OrderCalculator>
    inline int LPNorm< GridPart, OrderCalculator>::order(const int spaceOrder) const
    {
      return spaceOrder * orderCalc_->operator() (p_);
    }


    template< class GridPart, class OrderCalculator>
    template< class DiscreteFunctionType >
    inline typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type
    LPNorm< GridPart, OrderCalculator >::norm ( const DiscreteFunctionType &u ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0) );

      // return result
      return std::pow ( comm().sum( sum[ 0 ] ), (1.0 / p_) );
    }

    template< class GridPart, class OrderCalculator >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename Dune::FieldTraits< typename UDiscreteFunctionType::RangeFieldType >::real_type
    LPNorm< GridPart, OrderCalculator >
      ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0) );

      // return result
      return std::pow( comm().sum( sum[ 0 ] ), (1.0/p_) );
    }

    template< class GridPart, class OrderCalculator >
    template< class DiscreteFunctionType, class ReturnType >
    inline void
    LPNorm< GridPart, OrderCalculator >::
    normLocal ( const EntityType& entity,
                const unsigned int order,
                const DiscreteFunctionType &u,
                ReturnType& sum ) const
    {
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

      IntegratorType integrator( order );

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionMultiplicator< LocalFunctionType > ulocalp( ulocal, p_ );

      integrator.integrateAdd( entity, ulocalp, sum );
    }

    template< class GridPart, class OrderCalculator >
    template< class UDiscreteFunctionType,
              class VDiscreteFunctionType,
              class ReturnType >
    inline void
    LPNorm< GridPart, OrderCalculator >::
    distanceLocal ( const EntityType& entity,
                    const unsigned int order,
                    const UDiscreteFunctionType &u,
                    const VDiscreteFunctionType &v,
                    ReturnType& sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
      typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType >
        LocalDistanceType;

      IntegratorType integrator( order );

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      FunctionMultiplicator< LocalDistanceType > distp( dist, p_ );

      integrator.integrateAdd( entity, distp, sum );
    }


    template< class GridPart, class OrderCalculator >
    template< class Function >
    struct LPNorm< GridPart, OrderCalculator >::FunctionMultiplicator
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > RangeType ;

      explicit FunctionMultiplicator ( const FunctionType &function, double p )
      : function_( function ),
        p_(p)
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret = std :: pow ( phi.two_norm(), p_);
      }

    private:
      const FunctionType &function_;
      double p_;
    };


    template< class GridPart, class OrderCalculator >
    template< class UFunction, class VFunction >
    struct LPNorm< GridPart, OrderCalculator >::FunctionDistance
    {
      typedef UFunction UFunctionType;
      typedef VFunction VFunctionType;

      typedef typename UFunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
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


    // Implementation of WeightedL2Norm
    // --------------------------------

    template< class WeightFunction, class OrderCalculator >
    inline WeightedLPNorm< WeightFunction, OrderCalculator >
      ::WeightedLPNorm ( const WeightFunctionType &weightFunction, double p )
    : BaseType( weightFunction.space().gridPart(), p ),
      weightFunction_( weightFunction ),
      p_(p)
    {
      static_assert( (WeightFunctionSpaceType::dimRange == 1),
                          "Weight function must be scalar." );
    }


    template< class WeightFunction, class OrderCalculator >
    template< class UDiscreteFunctionType, class ReturnType >
    inline void
    WeightedLPNorm< WeightFunction, OrderCalculator >
      ::normLocal ( const EntityType &entity, const int order,
          const UDiscreteFunctionType &u,
          ReturnType &sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType LocalFunctionType;

      // !!!! order !!!!
      IntegratorType integrator( order );

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      LocalFunctionType ulocal = u.localFunction( entity );

      WeightedFunctionMultiplicator< LocalFunctionType > ulocal2( wflocal, ulocal, p_ );

      integrator.integrateAdd( entity, ulocal2, sum );
    }


    template< class WeightFunction,class OrderCalculator >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType >
    inline void
    WeightedLPNorm< WeightFunction, OrderCalculator >
      ::distanceLocal ( const EntityType &entity, const int order,
          const UDiscreteFunctionType &u,
          const VDiscreteFunctionType &v,
          ReturnType &sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
      typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

      typedef typename BaseType::template FunctionDistance
        < ULocalFunctionType, VLocalFunctionType >
        LocalDistanceType;

      // !!!! order !!!!
      IntegratorType integrator( order );

      LocalWeightFunctionType wflocal = weightFunction_.localFunction( entity );
      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      LocalDistanceType dist( ulocal, vlocal );
      WeightedFunctionMultiplicator< LocalDistanceType > dist2( wflocal, dist );

      integrator.integrateAdd( entity, dist2, sum );
    }


    template< class WeightFunction, class OrderCalculator>
    template< class Function >
    struct WeightedLPNorm< WeightFunction, OrderCalculator>::WeightedFunctionMultiplicator
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;
      typedef FieldVector< real_type, 1 > RangeType;

      WeightedFunctionMultiplicator ( const LocalWeightFunctionType &weightFunction,
                                      const FunctionType &function,
                                      double p )
      : weightFunction_( weightFunction ),
        function_( function ),
        p_(p)
      {}

      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        WeightType weight;
        weightFunction_.evaluate( x, weight );

        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        ret = weight[ 0 ] * std::pow ( phi.two_norm(), p_);
      }

    private:
      const LocalWeightFunctionType &weightFunction_;
      const FunctionType &function_;
      double p_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LPNORM_HH
