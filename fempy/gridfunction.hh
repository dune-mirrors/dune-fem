#ifndef DUNE_FEMPY_GRIDFUNCTION_HH
#define DUNE_FEMPY_GRIDFUNCTION_HH

#include <limits>
#include <string>
#include <utility>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

namespace Dune
{

  namespace FemPy
  {

    // GridFunctionBase
    // ----------------

    template< class GridPart >
    struct GridFunctionBase
      : public Fem::HasLocalFunction
    {};



    // RestrictProlongInterface
    // ------------------------

    template< class Grid >
    struct RestrictProlongInterface
    {
      typedef typename Grid::ctype ctype;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometry;

      virtual ~RestrictProlongInterface () = default;

      virtual void setFatherChildWeight ( const ctype &weight ) = 0;

      virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) = 0;
      virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) = 0;

      virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) = 0;
      virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) = 0;

      virtual void addToList ( Fem::CommunicationManagerList &commList ) = 0;
      virtual void removeFromList ( Fem::CommunicationManagerList &commList ) = 0;
    };



    // RestrictProlongImpl
    // -------------------

    template< class Grid, class Impl >
    class RestrictProlongImpl final
      : public RestrictProlongInterface< Grid >
    {
      typedef RestrictProlongInterface< Grid > Base;

    public:
      typedef typename Base::ctype ctype;

      typedef typename Base::Element Element;
      typedef typename Base::LocalGeometry LocalGeometry;

      template< class ... Args >
      RestrictProlongImpl ( Args && ... args )
        : impl_( std::forward< Args >( args ) ... )
      {}

      virtual void setFatherChildWeight ( const ctype &weight ) override
      {
        impl_.setFatherChildWeight( weight );
      }

      virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) override
      {
        impl_.restrictLocal( father, child, initialize );
      }

      virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) override
      {
        impl_.restrictLocal( father, child, geometryInFather, initialize );
      }

      virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) override
      {
        impl_.prolongLocal( father, child, initialize );
      }

      virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) override
      {
        impl_.prolongLocal( father, child, geometryInFather, initialize );
      }

      virtual void addToList ( Fem::CommunicationManagerList &commList ) override
      {
        //  impl_.addToList( commList );
      }

      virtual void removeFromList ( Fem::CommunicationManagerList &commList ) override
      {
        // impl_.removeFromList( commList );
      }

    private:
      Impl impl_;
    };



    // makeDefaultRestrictProlong
    // --------------------------

    template< class DiscreteFunction >
    inline static std::shared_ptr< RestrictProlongInterface< typename DiscreteFunction::GridPartType::GridType > >
    makeDefaultRestrictProlong ( DiscreteFunction &discreteFunction )
    {
      typedef RestrictProlongImpl< typename DiscreteFunction::GridPartType::GridType, Fem::RestrictProlongDefault< DiscreteFunction > > RP;
      return std::make_shared< RP >( discreteFunction );
    }



    // AdaptiveGridFunction
    // --------------------

    template< class Grid >
    struct AdaptiveGridFunction
    {
      typedef double DofType;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      typedef RestrictProlongInterface< Grid > RestrictProlong;

      explicit AdaptiveGridFunction ( std::shared_ptr< RestrictProlong > restrictProlong )
        : restrictProlong_( std::move( restrictProlong ) )
      {}

      virtual ~AdaptiveGridFunction () = default;

      virtual void enableDofCompression () = 0;

      virtual std::size_t numLocalDofs ( const Element &element ) const = 0;

      virtual DofType *getLocalDofs ( const Element &element, DofType *localDofs ) const = 0;
      virtual const DofType *setLocalDofs ( const Element &element, const DofType *localDofs ) = 0;

      RestrictProlong &restrictProlong () { assert( restrictProlong_ ); return *restrictProlong_; }

    protected:
      void resetRestrictProlong ( std::shared_ptr< RestrictProlong > restrictProlong )
      {
        restrictProlong_ = std::move( restrictProlong );
      }

    private:
      std::shared_ptr< RestrictProlong > restrictProlong_;
    };



    // GridFunction
    // ------------

    template< class GridPart, int dimR, class RangeField = double >
    struct GridFunction
      : public GridFunctionBase< GridPart >,
        public std::enable_shared_from_this< GridFunction< GridPart, dimR, RangeField > >
    {
      typedef RangeField RangeFieldType;
      typedef GridPart GridPartType;
      static const int dimRange = dimR;
      typedef Fem::FunctionSpace< double, RangeFieldType,
                                        GridPart::dimensionworld, dimRange > FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPartType::IntersectionType IntersectionType;
      typedef typename EntityType::Geometry::LocalCoordinate LocalDomainType;
      typedef typename Fem::CachingQuadrature< GridPartType, 0 >::QuadraturePointWrapperType Point;
      typedef typename Fem::CachingQuadrature< GridPartType, 1 >::QuadraturePointWrapperType IntersectionPoint;
      typedef typename Fem::ElementQuadrature< GridPartType, 0 >::QuadraturePointWrapperType ElementPoint;
      typedef typename Fem::ElementQuadrature< GridPartType, 1 >::QuadraturePointWrapperType ElementIntersectionPoint;

      virtual ~GridFunction () = default;

      virtual const std::string &name () const = 0;

      template< class Point >
      const DomainType xLocal ( const Point &point ) const
      {
        return Fem::coordinate( point );
      }

      template< class Point >
      const DomainType xGlobal ( const EntityType &entity, const Point &point ) const
      {
        return entity.geometry().global( xLocal( point ) );
      }

      struct LocalFunction
      {
        typedef RangeField RangeFieldType;
        LocalFunction ( const std::shared_ptr< const GridFunction< GridPart, dimRange > > &df )
          : df_( df ), object_( df_->attach( this ) )
        {}

        LocalFunction ( const std::shared_ptr< const GridFunction< GridPart, dimRange > > &df, const EntityType &entity )
          : df_( df ), object_( df_->attach( this ) )
        {
          init( entity );
        }

        LocalFunction ( const GridFunction< GridPart, dimRange > &df )
          : df_( df.shared_from_this()), object_( df_->attach( this ) )
        {}

        LocalFunction ( const GridFunction< GridPart, dimRange > &df, const EntityType &entity )
          : df_( df.shared_from_this()), object_( df_->attach( this ) )
        {
          init( entity );
        }

        ~LocalFunction ()
        {
          df_->detach( this );
        }

        void init ( const EntityType &entity )
        {
          entity_ = &entity;
          df_->init( this );
        }

        template< class Point >
        void evaluate ( const Point &point, RangeType &ret ) const
        {
          df_->evaluate( this, point, ret );
        }

        template< class Point >
        void jacobian ( const Point &point, JacobianRangeType &ret ) const
        {
          df_->jacobian( this, point, ret );
        }

        int order () const { return order_; }

        const EntityType &entity () const { return *entity_; }

        void *object () const { return object_; }

        std::shared_ptr< const GridFunction< GridPart, dimRange > > df_;
        mutable void *object_;
        mutable const EntityType  *entity_;
        int order_ = std::numeric_limits< int >::max();
      };

      LocalFunction localFunction ( const EntityType &entity ) const
      {
        return LocalFunction( this->shared_from_this(), entity );
      }

      virtual void evaluate( const LocalFunction *lf, const Point &x, RangeType &ret ) const = 0;
      virtual void evaluate( const LocalFunction *lf, const IntersectionPoint &x, RangeType &ret ) const = 0;
      virtual void evaluate( const LocalFunction *lf, const ElementPoint &x, RangeType &ret ) const = 0;
      virtual void evaluate( const LocalFunction *lf, const ElementIntersectionPoint &x, RangeType &ret ) const = 0;
      virtual void evaluate( const LocalFunction *lf, const LocalDomainType &x, RangeType &ret ) const = 0;

      template< class PointType >
      void evaluate ( const LocalFunction *lf, const PointType &x, RangeType &ret ) const
      {
        evaluate( lf, Fem::coordinate( x ), ret );
      }

      virtual void jacobian( const LocalFunction *lf, const Point &x, JacobianRangeType &ret ) const = 0;
      virtual void jacobian( const LocalFunction *lf, const IntersectionPoint &x, JacobianRangeType &ret ) const = 0;
      virtual void jacobian( const LocalFunction *lf, const ElementPoint &x, JacobianRangeType &ret ) const = 0;
      virtual void jacobian( const LocalFunction *lf, const ElementIntersectionPoint &x, JacobianRangeType &ret ) const = 0;
      virtual void jacobian( const LocalFunction *lf, const LocalDomainType &x, JacobianRangeType &ret ) const = 0;

      template< class PointType >
      void jacobian ( const LocalFunction *lf, const PointType &x, JacobianRangeType &ret ) const
      {
        jacobian( lf, Fem::coordinate( x ), ret );
      }

      // use this method to cache some data needed for a specific local function
      virtual void *attach ( LocalFunction *lf ) const { return nullptr; }
      virtual void detach ( LocalFunction *lf ) const {}
      virtual void init ( LocalFunction *lf ) const { lf->order_ = std::numeric_limits< int >::max(); }

      typedef LocalFunction LocalFunctionType;

      int getDimRange () { return dimRange; }

      struct Space
      {
        Space ( const int &order ) : order_( order ) {}
        int order () const { return order_; }
        const int order_;
      };

      Space space () const { return Space( order_ ); }

      int order_ = std::numeric_limits< int >::max();
    };


    template< class Container, class DiscreteFunction >
    struct GridFunctionDiscrete
      : public GridFunction< typename DiscreteFunction::GridPartType, DiscreteFunction::DiscreteFunctionSpaceType::dimRange, typename DiscreteFunction::DiscreteFunctionSpaceType::RangeFieldType >,
        public AdaptiveGridFunction< typename DiscreteFunction::GridPartType::GridType >
    {
      typedef DiscreteFunction DiscreteFunctionType;
      typedef typename DiscreteFunctionType::GridPartType GridPartType;
      static const int dimRange = DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
      typedef GridFunction< GridPartType, dimRange, RangeFieldType > Base;
      typedef typename Base::FunctionSpaceType FunctionSpaceType;
      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;
      typedef typename Base::JacobianRangeType JacobianRangeType;
      typedef typename Base::HessianRangeType HessianRangeType;
      typedef typename Base::DomainFieldType DomainFieldType;
      typedef typename Base::EntityType EntityType;
      typedef typename Base::IntersectionType IntersectionType;
      typedef typename Base::LocalDomainType LocalDomainType;
      typedef typename Base::Point Point;
      typedef typename Base::IntersectionPoint IntersectionPoint;
      typedef typename Base::ElementPoint ElementPoint;
      typedef typename Base::ElementIntersectionPoint ElementIntersectionPoint;
      typedef typename Base::LocalFunction LocalFunctionType;

      typedef typename AdaptiveGridFunction< typename DiscreteFunction::GridPartType::GridType >::DofType DofType;

      typename DiscreteFunctionType::LocalFunctionType &object ( const LocalFunctionType *lf ) const
      {
        return *( reinterpret_cast< typename DiscreteFunctionType::LocalFunctionType * >( lf->object()) );
      }

      GridFunctionDiscrete ( const std::shared_ptr< Container > &container, const DiscreteFunctionType &df )
        : AdaptiveGridFunction< typename GridPartType::GridType >( makeDefaultRestrictProlong( const_cast< DiscreteFunction & >( df ) ) ),
        container_( container ), df_( df )
      {
        Base::order_ = df_.space().order();
      }

      ~GridFunctionDiscrete ()
      {
        std::cout << "In GridFunctionDiscrete destructor" << std::endl;
        // Hack: ensure our restrict-prolong object is destroyed before the container
        AdaptiveGridFunction< typename GridPartType::GridType >::resetRestrictProlong( nullptr );
      }

      virtual const std::string &name () const { return df_.name(); }

      virtual void *attach ( LocalFunctionType *lf ) const
      {
        return (void *)( new typename DiscreteFunctionType::LocalFunctionType( df_ ) );
      }

      virtual void detach ( LocalFunctionType *lf ) const {}

      virtual void init ( LocalFunctionType *lf ) const
      {
        object( lf ).init( lf->entity() );
        lf->order_ = object( lf ).order();
      }

      virtual void evaluate ( const LocalFunctionType *lf, const Point &x, RangeType &ret ) const
      {
        object( lf ).evaluate( x, ret );
      }

      virtual void evaluate ( const LocalFunctionType *lf, const IntersectionPoint &x, RangeType &ret ) const
      {
        object( lf ).evaluate( x, ret );
      }

      virtual void evaluate ( const LocalFunctionType *lf, const ElementPoint &x, RangeType &ret ) const
      {
        object( lf ).evaluate( x, ret );
      }

      virtual void evaluate ( const LocalFunctionType *lf, const ElementIntersectionPoint &x, RangeType &ret ) const
      {
        object( lf ).evaluate( x, ret );
      }

      virtual void evaluate ( const LocalFunctionType *lf, const LocalDomainType &x, RangeType &ret ) const
      {
        object( lf ).evaluate( x, ret );
      }

      virtual void jacobian ( const LocalFunctionType *lf, const Point &x, JacobianRangeType &ret ) const
      {
        object( lf ).jacobian( x, ret );
      }

      virtual void jacobian ( const LocalFunctionType *lf, const IntersectionPoint &x, JacobianRangeType &ret ) const
      {
        object( lf ).jacobian( x, ret );
      }

      virtual void jacobian ( const LocalFunctionType *lf, const ElementPoint &x, JacobianRangeType &ret ) const
      {
        object( lf ).jacobian( x, ret );
      }

      virtual void jacobian ( const LocalFunctionType *lf, const ElementIntersectionPoint &x, JacobianRangeType &ret ) const
      {
        object( lf ).jacobian( x, ret );
      }

      virtual void jacobian ( const LocalFunctionType *lf, const LocalDomainType &x, JacobianRangeType &ret ) const
      {
        object( lf ).jacobian( x, ret );
      }

      RangeType eval ( const EntityType &entity, const LocalDomainType &x ) const
      {
        RangeType ret;
        df_.localFunction( entity ).evaluate( x, ret );
        return ret;
      }

      RangeType deriv ( int i, const EntityType &entity, const LocalDomainType &x ) const
      {
        RangeType ret;
        JacobianRangeType jac;
        df_.localFunction( entity ).jacobian( x, jac );
        for( unsigned int j = 0; j < ret.size(); ++j )
          ret[ j ] = jac[ j ][ i ];
        return ret;
      }

      virtual void enableDofCompression () final override { const_cast< DiscreteFunction & >( df_ ).enableDofCompression(); }

      virtual std::size_t numLocalDofs ( const EntityType &entity ) const final override
      {
        return df_.space().blockMapper().numDofs( entity ) * DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize;
      }

      virtual DofType *getLocalDofs ( const EntityType &entity, DofType *localDofs ) const final override
      {
        using Fem::dofBlockFunctor;

        Fem::AssignFunctor< DofType * > assignFunctor( localDofs );
        df_.space().blockMapper().mapEach( entity, dofBlockFunctor( df_.dofVector(), assignFunctor ) );

        return localDofs + numLocalDofs( entity );
      }

      virtual const DofType *setLocalDofs ( const EntityType &entity, const DofType *localDofs ) final override
      {
        using Fem::dofBlockFunctor;

        Fem::LeftAssign< const DofType * > assignFunctor( localDofs );
        df_.space().blockMapper().mapEach( entity, dofBlockFunctor( const_cast< DiscreteFunction & >( df_ ).dofVector(), assignFunctor ) );

        return localDofs + numLocalDofs( entity );
      }

      std::shared_ptr< Container > container_;
      const DiscreteFunctionType &df_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRIDFUNCTION_HH
