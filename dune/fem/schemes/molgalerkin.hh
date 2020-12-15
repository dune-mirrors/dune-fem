#ifndef DUNE_FEM_SCHEMES_MOLGALERKIN_HH
#define DUNE_FEM_SCHEMES_MOLGALERKIN_HH

// fem includes
#include <dune/fem/schemes/galerkin.hh>

namespace Dune
{

  namespace Fem
  {

    // GalerkinOperator
    // ----------------

    template< class Integrands, class DomainFunction, class RangeFunction = DomainFunction >
    struct MOLGalerkinOperator
      : public virtual Operator< DomainFunction, RangeFunction >
    {
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;


      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      typedef typename RangeFunctionType::GridPartType GridPartType;

      typedef Impl::GalerkinOperator< Integrands > GalerkinOperatorImplType;
      typedef typename RangeFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

      typedef typename GalerkinOperatorImplType::template QuadratureSelector<
            DiscreteFunctionSpaceType > :: InteriorQuadratureType    InteriorQuadratureType;

      typedef LocalMassMatrix< DiscreteFunctionSpaceType, InteriorQuadratureType >  LocalMassMatrixType ;

      template< class... Args >
      explicit MOLGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : impl_( gridPart, std::forward< Args >( args )... ),
          communicate_( true )
      {
        // disable communicate in Impl::GalerkinOperator
        // since applyInverseMass has to be applied first
        impl_.setCommunicate( false );
      }

      void setCommunicate( const bool communicate ) { communicate_ = communicate; }
      void setQuadratureOrders(unsigned int interior, unsigned int surface) { impl_.setQuadratureOrders(interior,surface); }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        evaluate( u, w );
      }

      const GridPartType &gridPart () const { return impl_.gridPart(); }

      typedef Integrands ModelType;
      typedef Integrands DirichletModelType;
      ModelType &model() const { return impl_.model(); }

    protected:
      void applyInverseMass( RangeFunctionType& w ) const
      {
        // get-set local contribution
        Dune::Fem::SetSelectedLocalContribution< RangeFunctionType > wLocal( w );

        LocalMassMatrixType localMassMatrix( w.space(), impl_.interiorQuadratureOrder( w.space().order() ) );

        // iterate over all elements
        for( const auto& entity : elements( gridPart(), Partitions::interiorBorder ) )
        {
          // fill local contribution
          auto guard = bindGuard( wLocal, entity );

          // apply inverse mass matrix
          // TODO: add mass term if needed (from ufl expression)
          localMassMatrix.applyInverse( entity, wLocal );
        }
      }

      template< class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        // Impl::GalerkinOperator::evaluate without communicate
        impl_.evaluate( u, w );

        // method of lines
        applyInverseMass( w );

        // synchronize data
        if( communicate_ )
          w.communicate();
      }

      // GalerkinOperator implementation (see galerkin.hh)
      GalerkinOperatorImplType impl_;
      bool communicate_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Integrands, class JacobianOperator >
    class MOLDifferentiableGalerkinOperator
      : public MOLGalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType >,
        public DifferentiableOperator< JacobianOperator >
    {
     typedef MOLGalerkinOperator< Integrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType > BaseType;

    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit MOLDifferentiableGalerkinOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                   const RangeDiscreteFunctionSpaceType &rSpace,
                                                   Args &&... args )
        : BaseType( rSpace.gridPart(), std::forward< Args >( args )... ),
          dSpace_(dSpace),
          rSpace_(rSpace)
      {}

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        // assemble Jacobian, same as GalerkinOperator
        impl_.assemble( u, jOp );
        // apply inverse mass
        applyInverseMass( jOp, impl_.model().hasSkeleton() );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        // assemble Jacobian, same as GalerkinOperator
        impl_.assemble( u, jOp );
        // apply inverse mass
        applyInverseMass( jOp, impl_.model().hasSkeleton() );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return rSpace_;
      }

      using BaseType::gridPart;

    protected:
      void applyInverseMass ( JacobianOperatorType &jOp, const bool hasSkeleton ) const
      {
        typedef typename BaseType::LocalMassMatrixType  LocalMassMatrixType;

        LocalMassMatrixType localMassMatrix( jOp.rangeSpace(), impl_.interiorQuadratureOrder( jOp.rangeSpace().order() ) );

        Dune::Fem::SetSelectedLocalContribution< JacobianOperatorType > jOpLocal( jOp );

        // multiply with inverse mass matrix
        for( const auto& inside : elements( gridPart(), Partitions::interiorBorder ) )
        {
          // scale diagonal
          {
            auto guard = bindGuard( jOpLocal, inside, inside );
            localMassMatrix.leftMultiplyInverse( jOpLocal );
          }

          if( hasSkeleton )
          {
            for( const auto &intersection : intersections( gridPart(), inside ) )
            {
              // scale off-diagonal
              if( intersection.neighbor() )
              {
                const auto& outside = intersection.outside();
                auto guard = bindGuard( jOpLocal, outside, inside );
                localMassMatrix.leftMultiplyInverse( jOpLocal );
              }
            }
          }
        }
      }

      using BaseType::impl_;
      const DomainDiscreteFunctionSpaceType &dSpace_;
      const RangeDiscreteFunctionSpaceType &rSpace_;
    };



    // AutomaticDifferenceGalerkinOperator
    // -----------------------------------

    template< class Integrands, class DomainFunction, class RangeFunction >
    class MOLAutomaticDifferenceGalerkinOperator
      : public MOLGalerkinOperator< Integrands, DomainFunction, RangeFunction >,
        public AutomaticDifferenceOperator< DomainFunction, RangeFunction >
    {
      typedef MOLGalerkinOperator< Integrands, DomainFunction, RangeFunction > BaseType;
      typedef AutomaticDifferenceOperator< DomainFunction, RangeFunction > AutomaticDifferenceOperatorType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit MOLAutomaticDifferenceGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : BaseType( gridPart, std::forward< Args >( args )... ), AutomaticDifferenceOperatorType()
      {}
    };



    // ModelDifferentiableGalerkinOperator
    // -----------------------------------

    template < class LinearOperator, class ModelIntegrands >
    struct MOLModelDifferentiableGalerkinOperator
      : public MOLDifferentiableGalerkinOperator< ModelIntegrands, LinearOperator >
    {
      typedef MOLDifferentiableGalerkinOperator< ModelIntegrands, LinearOperator > BaseType;

      typedef typename ModelIntegrands::ModelType ModelType;

      typedef typename LinearOperator::DomainFunctionType RangeFunctionType;
      typedef typename LinearOperator::RangeSpaceType DiscreteFunctionSpaceType;

      MOLModelDifferentiableGalerkinOperator ( ModelType &model, const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( dfSpace.gridPart(), model )
      {}

      template< class GridFunction >
      void apply ( const GridFunction &u, RangeFunctionType &w ) const
      {
        (*this)( u, w );
      }

      template< class GridFunction >
      void apply ( const GridFunction &u, LinearOperator &jOp ) const
      {
        (*this).jacobian( u, jOp );
      }
    };


    // MethodOfLinesScheme
    // -------------------

    template< class Integrands, class LinearOperator, class InverseOperator, bool addDirichletBC>
    using MethodOfLinesScheme = Impl::GalerkinSchemeImpl< Integrands, LinearOperator, InverseOperator, addDirichletBC,
                                                          MOLDifferentiableGalerkinOperator >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_MOLGALERKIN_HH
