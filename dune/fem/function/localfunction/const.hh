#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH

#include <algorithm>
#include <type_traits>
#include <utility>

#include <dune/common/dynvector.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/common/intersectionside.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declerations
    // -----------------------------

    template< class >
    struct DiscreteFunctionTraits;

    class HasLocalFunction;
    class IsDiscreteFunction;
    struct BindableFunction;

    // BasicConstLocalFunction
    // -----------------------

    template < class BasisFunctionSet, class LocalDofVector >
    class BasicConstLocalFunction
    : public LocalFunction< BasisFunctionSet, LocalDofVector >
    {
      typedef BasicConstLocalFunction< BasisFunctionSet, LocalDofVector >  ThisType;
      typedef LocalFunction< BasisFunctionSet, LocalDofVector > BaseType;

    public:
      //! type of Dof
      typedef typename BaseType::DofType DofType;

      //! type of Entity
      typedef typename BaseType :: EntityType EntityType;

      //! type of BasisFunctionSet
      typedef typename BaseType :: BasisFunctionSetType BasisFunctionSetType;

      //! type of LocalDofVector
      typedef typename BaseType :: LocalDofVectorType LocalDofVectorType;

      //! type of SizeType
      typedef typename BaseType::SizeType SizeType;

      //! default ctor
      BasicConstLocalFunction () {}

      explicit BasicConstLocalFunction ( const BasisFunctionSetType & basisFunctionSet ) : BaseType( basisFunctionSet ) {}

      explicit BasicConstLocalFunction ( const LocalDofVectorType &localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, const LocalDofVectorType &localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      explicit BasicConstLocalFunction ( LocalDofVectorType &&localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, LocalDofVectorType &&localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      BasicConstLocalFunction ( const BaseType &other ) : BaseType( other ) {}

      BasicConstLocalFunction ( const ThisType &other ) : BaseType( static_cast<const BaseType &>( other ) ) {}
      BasicConstLocalFunction ( ThisType && other ) : BaseType( static_cast<BaseType&&>(other) ) {}

      const DofType &operator[] ( SizeType i ) const { return static_cast< const BaseType & >( *this )[ i ]; }
      const DofType &operator[] ( SizeType i ) { return static_cast< const BaseType & >( *this )[ i ]; }

      using BaseType::localDofVector;

   protected:
      using BaseType::clear;
      using BaseType::assign;
      using BaseType::operator +=;
      using BaseType::operator -=;
      using BaseType::axpy;
    };

    /** \ingroup LocalFunction
        \class ConstLocalDiscreteFunction
        \brief A constant local function carrying values for one entity

        A ConstLocalDiscreteFunction is a LocalFunction which is basically doing the same as the
        LocalFunction of a discrete function. The difference is that the local dofs
        are not kept as references but are copied to a local storage.
        Therefore, this is a const local function and any modification of dofs is not
        allowed.

        \note Local DoF numbers correspond directly to array indices. Hence it
        may be more cache efficient to generate a ConstLocalFunction when only a
        const access to the local function is needed.

        \param DiscreteFunction type of the discrete function, the
                                local function shall belong to
     */
    template< class DiscreteFunction >
    class ConstLocalDiscreteFunction
    : public BasicConstLocalFunction<
      typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DiscreteFunctionSpaceType::BasisFunctionSetType,
      Dune::DynamicVector< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType,
        typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::LocalDofVectorAllocatorType
      :: template rebind< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > > ::DofType > ::other > >
    {
      typedef ConstLocalDiscreteFunction< DiscreteFunction > ThisType;
      typedef BasicConstLocalFunction< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DiscreteFunctionSpaceType::BasisFunctionSetType,
              Dune::DynamicVector< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType,
              typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > > :: LocalDofVectorAllocatorType
              :: template rebind< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType >::other  > >
          BaseType;

    public:
      typedef std::remove_const_t< DiscreteFunction > DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

      typedef DiscreteFunctionType GridFunctionType;

      typedef typename BaseType::DofType DofType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;
      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      /** \brief constructor creating a local function without binding it to an
                 entity

          Creates the local function without initializing the fields depending on
          the current entity.

          \note Before using the local function it must be initilized by
          \code
          localFunction.init( entity );
          \endcode

          \param[in] df discrete function the local function shall belong to
       */
      explicit ConstLocalDiscreteFunction ( const DiscreteFunctionType &df )
      : BaseType( LocalDofVectorType( df.localDofVectorAllocator() ) )
      , discreteFunction_( &df )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {}

      //! cast a MutableLocalFunction into this one !!! expensive !!!
      ConstLocalDiscreteFunction ( const typename DiscreteFunctionType::LocalFunctionType &localFunction )
      : BaseType( localFunction.basisFunctionSet(), LocalDofVectorType( localFunction.size(), localFunction.discreteFunction().localDofVectorAllocator() ) )
      , discreteFunction_( &localFunction.discreteFunction() )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {
        std::copy( localFunction.localDofVector().begin(), localFunction.localDofVector().end(), localDofVector().begin() );
      }

      /** \brief constructor creating a local function and binding it to an
                 entity

          Creates the local function and initilizes the fields depending on the
          current entity. It is not necessary, though allowed, to call init
          before using the discrete function.

          \note The degrees of freedom are not initialized by this function.

          \param[in] df      discrete function the local function shall
                             belong to
          \param[in] entity  entity for initialize the local function to
       */
      ConstLocalDiscreteFunction ( const DiscreteFunctionType &df, const EntityType &entity )
      : BaseType( df.space().basisFunctionSet( entity ), LocalDofVectorType( df.localDofVectorAllocator() )  )
      , discreteFunction_( &df )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }
      ConstLocalDiscreteFunction ( const EntityType &entity, const DiscreteFunctionType &df )
      : BaseType( df.space().basisFunctionSet( entity ), LocalDofVectorType( df.localDofVectorAllocator() )  )
      , discreteFunction_( &df )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

      //! copy constructor
      ConstLocalDiscreteFunction ( const ThisType &other )
      : BaseType( static_cast<const BaseType &>( other ) )
      , discreteFunction_( other.discreteFunction_ )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {}

      //! move constructor
      ConstLocalDiscreteFunction ( ThisType &&other )
      : BaseType( static_cast< BaseType &&>( other ) )
      , discreteFunction_( other.discreteFunction_ )
#ifdef TESTTHREADING
      , thread_(-1)
#endif
      {}

      using BaseType::localDofVector;

      using BaseType::evaluate;
      using BaseType::jacobian;
      using BaseType::hessian;

      /** \brief evaluate the local function
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \returns          value of the function in the given point
       */
      template< class Point >
      RangeType evaluate ( const Point &p ) const
      {
        RangeType val;
        evaluate( p, val );
        return val;
      }

      /** \brief evaluate Jacobian of the local function
       *
       *  \note Though the Jacobian is evaluated on the reference element, the
       *        return value is the Jacobian with respect to the actual entity.
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \returns          Jacobian of the function in the evaluation point
       */
      template< class Point >
      JacobianRangeType jacobian ( const Point &p ) const
      {
        JacobianRangeType jac;
        jacobian( p, jac );
        return jac;
      }

      /** \brief evaluate Hessian of the local function
       *
       *  \note Though the Hessian is evaluated on the reference element, the
       *        return value is the Hessian with respect to the actual entity.
       *
       *  \param[in]   x        evaluation point in local coordinates
       *  \returns              Hessian of the function in the evaluation point
       */
      template< class Point >
      HessianRangeType hessian ( const Point &p ) const
      {
        HessianRangeType h;
        hessian( p, h );
        return h;
      }

      /** \copydoc Dune::Fem::LocalFunction :: init */
      void init ( const EntityType &entity )
      {
#ifdef TESTTHREADING
        if (thread_ == -1) thread_ = MPIManager::thread();
        if (thread_ != MPIManager::thread())
        {
          std::cout << "wrong thread number\n";
          assert(0);
          std::abort();
        }
#endif
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

      void bind ( const EntityType &entity ) { init( entity ); }
      void unbind ()
      {
#ifdef TESTTHREADING
        if (thread_ != MPIManager::thread())
        {
          std::cout << "wrong thread number\n";
          assert(0);
          std::abort();
        }
#endif
      }
      void bind(const IntersectionType &intersection, IntersectionSide side)
      {
        // store local copy to avoid problems with casting to temporary types
        const EntityType entity = side==IntersectionSide::in? intersection.inside(): intersection.outside();
        bind( entity );
      }

      const DiscreteFunctionType &discreteFunction() const { return *discreteFunction_; }
      const GridFunctionType &gridFunction() const { return discreteFunction(); }

    protected:
      const DiscreteFunctionType* discreteFunction_;
#ifdef TESTTHREADING
      int thread_;
#endif
    };



    // ConstLocalFunction
    // ------------------

    namespace Impl
    {

      template< class GF, class = void >
      struct ConstLocalFunction;

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value && std::is_base_of< Fem::IsDiscreteFunction, GF >::value > >
      {
        typedef ConstLocalDiscreteFunction< GF > Type;
      };

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value && !std::is_base_of< Fem::IsDiscreteFunction, GF >::value  && std::is_class< typename GF::LocalFunctionType >::value > >
      {
        struct Type
          : public GF::LocalFunctionType
        {
          typedef GF GridFunctionType;
          typedef typename GridFunctionType::LocalFunctionType::EntityType EntityType;

          typedef typename GF::LocalFunctionType::DomainType DomainType;
          typedef typename GF::LocalFunctionType::RangeType RangeType;
          typedef typename GF::LocalFunctionType::JacobianRangeType JacobianRangeType;
          typedef typename GF::LocalFunctionType::HessianRangeType HessianRangeType;

          explicit Type ( const GridFunctionType &gridFunction )
            : GridFunctionType::LocalFunctionType( gridFunction ),
              gridFunction_( gridFunction ) {}
          explicit Type ( const EntityType &entity, const GridFunctionType &gridFunction )
            : GridFunctionType::LocalFunctionType( gridFunction ),
              gridFunction_( gridFunction )
          { bind(entity); }

          using GF::LocalFunctionType::evaluate;
          using GF::LocalFunctionType::jacobian;
          using GF::LocalFunctionType::hessian;
          using GF::LocalFunctionType::init;
          using GF::LocalFunctionType::entity;

          //! evaluate local function
          template< class Point >
          RangeType evaluate ( const Point &p ) const
          {
            RangeType val;
            evaluate( p, val );
            return val;
          }

          //! jacobian of local function
          template< class Point >
          JacobianRangeType jacobian ( const Point &p ) const
          {
            JacobianRangeType jac;
            jacobian( p, jac );
            return jac;
          }

          //! hessian of local function
          template< class Point >
          HessianRangeType hessian ( const Point &p ) const
          {
            HessianRangeType h;
            hessian( p, h );
            return h;
          }

          void bind ( const EntityType &entity ) { init( entity ); }
          void unbind () {}
          template <class IntersectionType>
          void bind(const IntersectionType &intersection, IntersectionSide side)
          {
            // store local copy to avoid problems with casting to temporary types
            const EntityType entity = side==IntersectionSide::in? intersection.inside(): intersection.outside();
            bind( entity );
          }

          const GridFunctionType &gridFunction () const { return gridFunction_; }

        private:
          const GridFunctionType &gridFunction_;
        };
      };

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::BindableFunction, std::decay_t<GF> >::value && !std::is_base_of< Fem::IsDiscreteFunction, std::decay_t<GF> >::value > >
      {
        static_assert( !std::is_reference<GF>::value );
        struct Type
        {
          typedef GF GridFunctionType;
          typedef std::decay_t<GF> GridFunctionDecayType;
          typedef typename GridFunctionDecayType::GridPartType GridPartType;
          typedef typename GridFunctionDecayType::EntityType EntityType;
          typedef typename GridFunctionDecayType::RangeFieldType RangeFieldType;
          typedef typename GridFunctionDecayType::DomainType DomainType;
          typedef typename GridFunctionDecayType::RangeType RangeType;
          typedef typename GridFunctionDecayType::JacobianRangeType JacobianRangeType;
          typedef typename GridFunctionDecayType::HessianRangeType HessianRangeType;
          typedef typename GridFunctionDecayType::FunctionSpaceType FunctionSpaceType;

          template<class Arg, std::enable_if_t<std::is_constructible<GF, Arg>::value, int> = 0>
          explicit Type ( Arg&& gridFunction )
            : gridFunction_( std::forward<Arg>(gridFunction) )
#ifdef TESTTHREADING
            , thread_(-1)
#endif
          {
            // if (MPIManager::thread()==0 || MPIManager::thread()==1) std::cout << "[" << MPIManager::thread() << "]: CLF " << &gridFunction_ << std::endl;
          }
          template<class Arg, std::enable_if_t<std::is_constructible<GF, Arg>::value, int> = 0>
          explicit Type ( const EntityType &entity, Arg&& gridFunction )
            :  gridFunction_( std::forward<Arg>(gridFunction) )
#ifdef TESTTHREADING
            ,  thread_(-1)
#endif
          {
            bind(entity);
          }

          template <class Point>
          void evaluate(const Point &x, RangeType &ret) const
          {
            gridFunction().evaluate(x,ret);
          }
          template <class Point>
          void jacobian(const Point &x, JacobianRangeType &ret) const
          {
            gridFunction().jacobian(x,ret);
          }
          template <class Point>
          void hessian(const Point &x, HessianRangeType &ret) const
          {
            gridFunction().hessian(x,ret);
          }
          unsigned int order() const { return gridFunction().order(); }

          //! evaluate local function
          template< class Point >
          RangeType evaluate ( const Point &p ) const
          {
            RangeType val;
            evaluate( p, val );
            return val;
          }

          //! jacobian of local function
          template< class Point >
          JacobianRangeType jacobian ( const Point &p ) const
          {
            JacobianRangeType jac;
            jacobian( p, jac );
            return jac;
          }

          //! hessian of local function
          template< class Point >
          HessianRangeType hessian ( const Point &p ) const
          {
            HessianRangeType h;
            hessian( p, h );
            return h;
          }

          template< class Quadrature, class ... Vectors >
          void evaluateQuadrature ( const Quadrature &quad, Vectors & ... values ) const
          {
            static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
            evaluateFullQuadrature( PriorityTag<42>(), quad, values... );
          }
          template< class Quadrature, class Jacobians >
          void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
          { jacobianQuadrature(quadrature,jacobians, PriorityTag<42>() ); }
          template< class Quadrature, class Hessians >
          void hessianQuadrature ( const Quadrature &quadrature, Hessians &hessians ) const
          { hessianQuadrature(quadrature,hessians, PriorityTag<42>() ); }

          void bind ( const EntityType &entity )
          {
#ifdef TESTTHREADING
            if (thread_ == -1) thread_ = MPIManager::thread();
            if (thread_ != MPIManager::thread())
            {
              std::cout << "wrong thread number\n";
              assert(0);
              std::abort();
            }
#endif
            gridFunction_.bind( entity );
          }
          void unbind ()
          {
#ifdef TESTTHREADING
            if (thread_ != MPIManager::thread())
            {
              std::cout << "wrong thread number\n";
              assert(0);
              std::abort();
            }
#endif
            gridFunction_.unbind();
          }
          template <class IntersectionType>
          void bind(const IntersectionType &intersection, IntersectionSide side)
          {
#ifdef TESTTHREADING
            if (thread_ == -1) thread_ = MPIManager::thread();
            if (thread_ != MPIManager::thread())
            {
              std::cout << "wrong thread number\n";
              assert(0);
              std::abort();
            }
#endif
            defaultIntersectionBind(gridFunction_,intersection, side);
          }

          const EntityType& entity() const
          {
            return gridFunction_.entity();
          }

          const GridFunctionDecayType &gridFunction () const { return gridFunction_; }

        private:
          template< class Quadrature, class ... Vectors, class GF_=GridFunctionDecayType >
          auto evaluateFullQuadrature ( PriorityTag<1>, const Quadrature &quad, Vectors & ... values ) const
          -> std::enable_if_t< std::is_void< decltype( std::declval< const GF_& >().evaluateQuadrature(quad,values...))>::value >
          { gridFunction().evaluateQuadrature(quad,values...); }
          template< class Quadrature, class ... Vectors >
          void evaluateFullQuadrature ( PriorityTag<0>, const Quadrature &quad, Vectors & ... values ) const
          { std::ignore = std::make_tuple( ( evaluateSingleQuadrature( quad, values ), 1 ) ... ); }

          template< class Quadrature, class Jacobians, class GF_=GridFunctionDecayType>
          auto jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians, PriorityTag<1> ) const
          -> std::enable_if_t< std::is_void< decltype( std::declval< const GF_& >().jacobianQuadrature(quadrature,jacobians))>::value >
          { gridFunction().jacobianQuadrature(quadrature,jacobians); }
          template< class Quadrature, class Jacobians >
          void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians, PriorityTag<0> ) const
          {
            for( const auto qp : quadrature )
              jacobians[ qp.index() ] = jacobian( qp );
          }

          template< class Quadrature, class Hessians, class GF_=GridFunctionDecayType >
          auto hessianQuadrature ( const Quadrature &quadrature, Hessians &hessians, PriorityTag<1> ) const
          -> std::enable_if_t< std::is_void< decltype( std::declval< const GF_& >().hessianQuadrature(quadrature,hessians))>::value >
          { gridFunction().hessianQuadrature(quadrature,hessians); }
          template< class Quadrature, class Hessians >
          void hessianQuadrature ( const Quadrature &quadrature, Hessians &hessians, PriorityTag<0> ) const
          {
            for( const auto qp : quadrature )
              hessians[ qp.index() ] = hessian( qp );
          }

          template< class Quadrature, class Vector >
          auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
          -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, RangeType >::value >
          {
            for( const auto qp : quad )
              v[ qp.index() ] = evaluate( qp );
          }
          template< class Quadrature, class Vector >
          auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
          -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, JacobianRangeType >::value >
          { jacobianQuadrature(quad,v); }
          template< class Quadrature, class Vector >
          auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
          -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, HessianRangeType >::value >
          { hessianQuadrature(quad,v); }

          GridFunctionType gridFunction_;
#ifdef TESTTHREADING
          int thread_;
#endif
        };
      };
    } // namespace Impl


    template< class GridFunction >
      using ConstLocalFunction = typename Impl::ConstLocalFunction< GridFunction >::Type;
    /**@internal Default FalseType.*/
    template<class T, class SFINAE = void>
      struct IsConstLocalFunction
      : std::false_type
      {};

    /**@internal Forward to decay_t.*/
    template<class T>
      struct IsConstLocalFunction<T, std::enable_if_t<!std::is_same<T, std::decay_t<T> >{}> >
      : IsConstLocalFunction<std::decay_t<T> >
      {};

    /**TrueType if a T can be wrapped into a Fem::ConstLocalFunction.*/
    template<class T>
      struct IsConstLocalFunction<
      T,
      std::enable_if_t<(std::is_same<T, std::decay_t<T> >{}
          && std::is_same<T, Fem::ConstLocalFunction<typename T::GridFunctionType> >{}
          )> >
        : std::true_type
        {};


    /**Wrap an F into a Fem::ConstLocalFunction if directly allowed
     * by Fem::ConstLocalFunction.
     */
    template<class F, std::enable_if_t<!IsConstLocalFunction<F>::value, int> = 0>
    constexpr auto constLocalFunction(F&& f)
    {
      return Fem::ConstLocalFunction<std::decay_t<F> >(std::forward<F>(f));
    }

    /**@internal Forward a Fem::ConstLocalFunction as is.*/
    template<class F, std::enable_if_t<IsConstLocalFunction<F>::value, int> = 0>
    constexpr decltype(auto) constLocalFunction(F&& f)
    {
      return std::forward<F>(f);
    }

    template<class F, class Entity>
    constexpr auto constLocalFunction(F&& f, const Entity &entity)
    {
      return Dune::Fem::ConstLocalFunction<std::decay_t<F> >(entity,std::forward<F>(f));
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
