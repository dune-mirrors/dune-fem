#ifndef DUNE_FEM_GRIDFUNCTIONADAPTER_HH
#define DUNE_FEM_GRIDFUNCTIONADAPTER_HH

#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>

//- local includes
#include <dune/fem/version.hh>
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/function/common/discretefunction.hh>

// for compatibility
#include <dune/fem/function/common/localfunctionadapter.hh>

namespace Dune
{

  namespace Fem
  {

    /** @addtogroup GridFunctionAdapter

        To plug an \ref Fem::Function "analytical function"
        into a operator taking
        \ref DiscreteFunctionInterface "discrete functions",
        i.e., expecting \ref LocalFunction "local functions"
        a wrapper can be applied to the analytical function
        instance.
        The resulting class is still a \ref Fem::Function "Function"
        but with the property \ref HasLocalFunction "\" has local function \"" added.

        @{
    */

    //- forward declaration
    template <class FunctionImp, class GridPartImp>
    class BasicGridFunctionAdapter;

    //! traits of GridFunctionAdapter
    template <class FunctionImp, class GridPartImp>
    struct BasicGridFunctionAdapterTraits
    {
      typedef typename std::decay_t< FunctionImp >::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      typedef DiscreteFunctionSpaceAdapter<FunctionSpaceType,GridPartImp> DiscreteFunctionSpaceType;

      typedef GridPartImp GridPartType;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridPartType :: template Codim<0> :: EntityType EntityType;
      typedef typename GridPartType :: IntersectionType IntersectionType;
      //! type of iterator
      typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
      //! type of IndexSet
      typedef typename GridPartType :: IndexSetType IndexSetType;

      typedef BasicGridFunctionAdapter<FunctionImp,GridPartImp> DiscreteFunctionType;
    };



    // BasicGridFunctionAdapter
    // ------------------------

    /** \brief BasicGridFunctionAdapter provides local functions for a Function.
     */
    template< class FunctionImp, class GridPartImp >
    class BasicGridFunctionAdapter
    : public Function< typename std::decay_t< FunctionImp >::FunctionSpaceType,
                       BasicGridFunctionAdapter< FunctionImp, GridPartImp > >,
      public HasLocalFunction
    {
      typedef BasicGridFunctionAdapter< FunctionImp, GridPartImp > ThisType;
      typedef Function< typename std::decay_t< FunctionImp >::FunctionSpaceType, ThisType > BaseType;

      // Make sure the function is not a discrete functon
      static_assert( !(std::is_convertible< FunctionImp, HasLocalFunction >::value),
                     "FunctionType may not be a discrete function type." );

    public:
      //! type of traits
      typedef BasicGridFunctionAdapterTraits< FunctionImp, GridPartImp > Traits;

      //! type of function
      typedef std::decay_t< FunctionImp > FunctionType;

      //! type of grid part
      typedef GridPartImp GridPartType;

      //! type of discrete function space
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      // type of discrete function space
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;

      //! type of grid
      typedef typename DiscreteFunctionSpaceType::GridType GridType;

      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
      //! jacobian type (from function space)
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

      //! type of codim 0 entity
      typedef typename Traits :: EntityType  EntityType;
      typedef typename Traits :: IntersectionType IntersectionType;

    private:
      class LocalFunction;

    public:
      //! type of local function to export
      typedef LocalFunction LocalFunctionType;

      // reference to function this local belongs to
      BasicGridFunctionAdapter ( std::string name, FunctionImp f, const GridPartType &gridPart, unsigned int order = DiscreteFunctionSpaceType::polynomialOrder )
      : space_( gridPart, order ),
        function_( std::move( f ) ),
        name_( std::move( name ) )
      {}

      // reference to function this local belongs to
      BasicGridFunctionAdapter ( const ThisType &other )
      : space_( other.space_ ),
        function_( other.function_ ),
        name_( other.name_ )
      {}

      //! evaluate function on local coordinate local
      void evaluate ( const DomainType &global, RangeType &result ) const
      {
        function_.evaluate( global, result );
      }

      //! evaluate function on local coordinate local
      void jacobian ( const DomainType &global, JacobianRangeType &result ) const
      {
        function_.jacobian(global,result);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      LocalFunctionType localFunction ( const EntityType &entity )
      {
        return LocalFunctionType( entity, *this );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */
      const LocalFunctionType localFunction ( const EntityType &entity ) const
      {
        return LocalFunctionType( entity, *this );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name() const */
      const std::string &name () const
      {
        return name_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::space() const */
      const DiscreteFunctionSpaceType &space () const
      {
        return space_;
      }

      const GridPartType &gridPart () const
      {
        return space().gridPart();
      }

      //! return true, probably
      inline int order () const
      {
        return space().order();
      }

      //! return true, probably
      inline bool continuous () const
      {
        return space().continuous();
      }

    protected:
      DiscreteFunctionSpaceType space_;
      FunctionImp function_;
      const std::string name_;
    };



    // BasicGridFunctionAdapter::LocalFunction
    // ---------------------------------------

    template< class Function, class GridPart >
    class BasicGridFunctionAdapter< Function, GridPart >::LocalFunction
    {
      typedef LocalFunction ThisType;
      typedef BasicGridFunctionAdapter< Function, GridPart > DiscreteFunctionType;

    public:
      //! function space type
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;

      //! domain field type (from function space)
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      //! range field type (from function space)
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      //! domain dimension (from function space)
      static const int dimDomain = GridPart::GridType::dimensionworld;
      //! range dimension (from function space)
      static const int dimRange = FunctionSpaceType::dimRange;

      //! domain type (from function space)
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! range type (from function space)
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! jacobian type (from function space)
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! hessian type (from function space)
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! entity type
      typedef typename Traits::EntityType EntityType;
      typedef typename Traits::IntersectionType IntersectionType;
      //! local coordinate type
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
      //! local dimension
      static const int dimLocal = LocalCoordinateType::dimension;

      //! constructor initializing local function
      LocalFunction ( const EntityType &entity, const DiscreteFunctionType &df )
      : function_( &df.function_ ),
        entity_( &entity ),
        order_( df.space().order() )
      {}

      LocalFunction ( const DiscreteFunctionType &df )
      : function_( &df.function_ ),
        entity_( 0 ),
        order_( df.space().order() )
      {}

      //! copy constructor
      LocalFunction ( const LocalFunction &other ) = default;

      //! evaluate local function
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const
      {
        const auto geometry = entity().geometry();
        auto global = geometry.global( coordinate( x ) );
        function().evaluate( global, ret );
      }
      template< class PointType >
      RangeType operator() ( const PointType &x ) const
      {
        RangeType ret;
        evaluate(x,ret);
        return ret;
      }

      //! jacobian of local function
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        const auto geometry = entity().geometry();
        auto global = geometry.global( coordinate( x ) );
        function().jacobian( global, ret );

        if( dimLocal != dimDomain )
        {
          // This computes the projection to the tangential space
          // (i.e. the hyperplane this entity is contained in). This
          // is done in a generic way by first projecting to the local
          // tangential space of the reference elment, and then
          // projecting back to the ambient space.

          const auto gjt = geometry.jacobianTransposed( coordinate( x ) );
          const auto gjit = geometry.jacobianInverseTransposed( coordinate( x ) );

          FieldVector< RangeFieldType, dimLocal > tmp;
          for( auto i = 0; i < dimRange; ++i )
          {
            gjit.mtv( ret[ i ], tmp );
            gjt.mtv( tmp, ret[ i ] );
          }
        }
      }

      //! hessian of local function
      template< class PointType >
      void hessian ( const PointType &x, HessianRangeType &ret ) const
      {
        DUNE_THROW( NotImplemented, "Method hessian() not implemented yet" );
      }

      //! evaluate function or jacobian of function for given quadrature
      template < class QuadratureType, class ... Vectors >
      void evaluateQuadrature( const QuadratureType& quadrature, Vectors& ... values ) const
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore = std::make_tuple( ( evaluateQuadratureImp( quadrature, values ), 1 ) ... );
      }

      int order () const { return order_; }

      //! init local function
      void init ( const EntityType &entity )
      {
        entity_ = &entity;
      }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

    protected:
      template < class QuadratureType, class VectorType >
      auto evaluateQuadratureImp( const QuadratureType& quadrature, VectorType& values ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(values[ 0 ] ) >, RangeType >::value >
      {
        for( auto qp : quadrature )
          evaluate( qp, values[ qp.index() ] );
      }

      template < class QuadratureType, class VectorType >
      auto evaluateQuadratureImp( const QuadratureType& quadrature, VectorType& values ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(values[ 0 ] ) >, JacobianRangeType >::value >
      {
        for( auto qp : quadrature )
          jacobian( qp, values[ qp.index() ] );
      }

      const FunctionType &function () const
      {
        return *function_;
      }

      const FunctionType *function_;
      const EntityType *entity_;
      int order_;
    };



    // GridFunctionAdapter
    // -------------------

    template< class Function, class GridPart >
    using GridFunctionAdapter = BasicGridFunctionAdapter< const Function &, GridPart >;



    // gridFunctionAdapter
    // -------------------

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  name      name of the grid function
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only lvalue references. The grid function only
     *       references the original function.
     **/
    template< class Function, class GridPart >
    inline static GridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( std::string name, const Function &function, const GridPart &gridPart, unsigned int order )
    {
      return GridFunctionAdapter< Function, GridPart >( std::move( name ), function, gridPart, order );
    }

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only lvalue references. The grid function only
     *       references the original function.
     **/
    template< class Function, class GridPart >
    inline static GridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( const Function &function, const GridPart &gridPart, unsigned int order )
    {
      return GridFunctionAdapter< Function, GridPart >( std::string(), function, gridPart, order );
    }

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  name      name of the grid function
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only lvalue references. The grid function only
     *       references the original function.
     **/
    template< class Function, class GridPart >
    inline static GridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( std::string name, Function &function, const GridPart &gridPart, unsigned int order )
    {
      const Function& cf = function;
      return gridFunctionAdapter( name, cf, gridPart, order );
    }

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only lvalue references. The grid function only
     *       references the original function.
     **/
    template< class Function, class GridPart >
    inline static GridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( Function &function, const GridPart &gridPart, unsigned int order )
    {
      const Function& cf = function;
      return gridFunctionAdapter( cf, gridPart, order );
    }

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  name      name of the grid function
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only rvalue references. The original function
     *       is move-constructed into the grid function.
     **/
    template< class Function, class GridPart >
    inline static BasicGridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( std::string name, Function &&function, const GridPart &gridPart, unsigned int order )
    {
      return BasicGridFunctionAdapter< Function, GridPart >( std::move( name ), std::move( function ), gridPart, order );
    }

    /**
     * \brief convert a function to a grid function
     *
     * \param[in]  function  function to convert
     * \param[in]  gridPart  grid part to restrict the domain to
     * \param[in]  order     polynomial order to report
     *
     * \note This version accepts only rvalue references. The original function
     *       is move-constructed into the grid function.
     **/
    template< class Function, class GridPart >
    inline static BasicGridFunctionAdapter< Function, GridPart >
    gridFunctionAdapter ( Function &&function, const GridPart &gridPart, unsigned int order )
    {
      return BasicGridFunctionAdapter< Function, GridPart >( std::string(), std::move( function ), gridPart, order );
    }



    namespace
    {
      template <class FunctionImp,class GridPartType,bool>
      struct ConvertDFTypeHelper;

      template <class FunctionImp,class GridPartType>
      struct ConvertDFTypeHelper<FunctionImp,GridPartType,true>
      {
        typedef ConvertDFTypeHelper<FunctionImp,GridPartType,true> ThisType;
        enum {compatible = std::is_convertible<GridPartType,typename FunctionImp::DiscreteFunctionSpaceType::GridPartType>::value};
        typedef FunctionImp FunctionType;
        typedef typename FunctionType::DiscreteFunctionSpaceType DFSType;
        ConvertDFTypeHelper(const std::string& name,const FunctionImp& func,const GridPartType& gp) :
          func_(func)
        {}
        ConvertDFTypeHelper(const ConvertDFTypeHelper& other) :
          func_(other.func_)
        {}
        const FunctionType& function() const
        {
          return func_;
        }
        const DFSType& space() const
        {
          return func_.space();
        }
        private:
        const FunctionImp& func_;
      };

      template <class FunctionImp,class GridPartType>
      struct ConvertDFTypeHelper<FunctionImp,GridPartType,false>
        : GridFunctionAdapter<FunctionImp,GridPartType>
      {
        typedef ConvertDFTypeHelper<FunctionImp,GridPartType,false> ThisType;
        typedef GridFunctionAdapter<FunctionImp,GridPartType> BaseType;
        typedef BaseType FunctionType;
        typedef typename FunctionType::DiscreteFunctionSpaceType DFSType;
        ConvertDFTypeHelper(const std::string& name,const FunctionImp& func,const GridPartType& gp) :
          BaseType(name,func,gp)
        {}
        ConvertDFTypeHelper(const ConvertDFTypeHelper& other) :
          BaseType(other)
        {}
        const FunctionType& function() const
        {
          return *this;
        }
        const DFSType& space() const
        {
          return BaseType::space();
        }
      };
    }

    template< class FunctionImp, class GridPartImp >
    class ConvertToGridFunction
    : public Function< typename FunctionImp::FunctionSpaceType,
                       ConvertToGridFunction< FunctionImp, GridPartImp > >,
      public HasLocalFunction
    {
      typedef ConvertToGridFunction< FunctionImp, GridPartImp > ThisType;
      typedef Function< typename FunctionImp::FunctionSpaceType, ThisType > BaseType;
      static const bool hasLocalFunction = std::is_convertible< FunctionImp, HasLocalFunction >::value;
      typedef ConvertDFTypeHelper< FunctionImp, GridPartImp, hasLocalFunction > Helper;
      typedef typename Helper::FunctionType ConvertedType;

    public:
      typedef FunctionImp FunctionType;
      typedef GridPartImp GridPartType;

      //! type of discrete function space
      typedef typename ConvertedType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      // type of discrete function space
      typedef typename ConvertedType::FunctionSpaceType FunctionSpaceType;

      //! type of grid
      typedef typename DiscreteFunctionSpaceType::GridType GridType;

      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType ;
      //! domain type (from function space)
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType ;
      //! range type (from function space)
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType ;
      //! jacobian type (from function space)
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

      //! type of codim 0 entity
      typedef typename GridPartType :: template Codim<0> :: EntityType  EntityType;

      //! type of local function to export
      typedef typename ConvertedType::LocalFunctionType LocalFunctionType;

      //! constructor
      ConvertToGridFunction ( const std::string &name,
                              const FunctionImp &function,
                              const GridPartType &gridPart )
      : name_( name ),
        helper_( name, function, gridPart )
      {}

      ConvertToGridFunction ( const ThisType &other )
      : name_( other.name_ ),
        helper_( other.helper_ )
      {}

      //! evaluate function on local coordinate local
      void evaluate ( const DomainType &global, RangeType &result ) const
      {
        helper_.function().evaluate(global,result);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */
      const LocalFunctionType localFunction( const EntityType &entity ) const
      {
        return helper_.function().localFunction(entity);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
      LocalFunctionType localFunction( const EntityType &entity )
      {
        return helper_.function().localFunction(entity);
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::name */
      const std::string &name() const
      {
        return name_;
      }

      const DiscreteFunctionSpaceType &space() const
      {
        return helper_.space();
      }

    private:
      const std::string name_;
      Helper helper_;
    };

    template< class Function, class GridPart >
    inline ConvertToGridFunction< Function, GridPart >
    convertToGridFunction ( const std::string &name,
                            const Function &function,
                            const GridPart &gridPart )
    {
      return ConvertToGridFunction< Function, GridPart >( name, function, gridPart );
    }

  } // namespace Fem

} // namespace Dune

//@}

#endif // #ifndef DUNE_DISCRETEFUNCTIONADAPTER_HH
