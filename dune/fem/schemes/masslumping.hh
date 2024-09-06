#ifndef DUNE_FEM_SCHEMES_MASSLUMPING_HH
#define DUNE_FEM_SCHEMES_MASSLUMPING_HH

// fem includes
#include <dune/fem/schemes/galerkin.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

namespace Dune
{

  namespace Fem
  {

#if HAVE_PETSC
    //forward declaration of PetscLinearOperator
    template <typename DomainFunction, typename RangeFunction >
    class PetscLinearOperator ;
#endif

    // MassLumpingOperator
    // -------------------

    template< class Integrands, class MassIntegrands, class DomainFunction, class RangeFunction = DomainFunction >
    struct MassLumpingOperator
      : public virtual Operator< DomainFunction, RangeFunction >
    {
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;


      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      typedef typename RangeFunctionType::GridPartType GridPartType;
      typedef ThreadIterator< GridPartType >  ThreadIteratorType;

      typedef Impl::GalerkinOperator< GridPartType >    GalerkinOperatorImplType;
      typedef Impl::LocalGalerkinOperator< Integrands > LocalGalerkinOperatorImplType;

      template <class Space>
      struct LumpingQuadratureSelector
      {
        typedef typename Space :: GridPartType GridPartType;
        typedef CachingLumpingQuadrature< GridPartType, 0 > InteriorQuadratureType;

        // not needed, since lumping is only on element terms
        typedef CachingQuadrature< GridPartType, 1, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > SurfaceQuadratureType;
      };

      // define Galerkin operator with different quadratures (i.e. LumpingQuadrature)
      typedef Impl::LocalGalerkinOperator< MassIntegrands, LumpingQuadratureSelector > MassOperatorImplType;

      typedef typename RangeFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

      template< class... Args >
      explicit MassLumpingOperator ( const GridPartType &gridPart,
                                     const Integrands& integrands,
                                     const MassIntegrands& massIntegrands,
                                     Args &&... args )
        : iterators_( gridPart ),
          opImpl_( gridPart ),
          localOp_( gridPart, std::move( integrands ), std::forward< Args >( args )... ),
          mass_( gridPart, std::move( massIntegrands ), std::forward< Args >( args )... ),
          gridSizeInterior_( 0 ),
          communicate_( true )
      {
        if( mass().hasSkeleton() )
          DUNE_THROW(InvalidStateException,"MassLumpingOperator: Mass model cannot have a skeleton term");

        // volume order (same as space for lumping)
        auto interiorOrder = [] (const int order) { return order; };
        // surface order does not matter for this part
        auto surfaceOrder  = [] (const int order) { return 0; };

        // set quadrature orders mass part to be fixed!
        size_t size = mass_.size();
        for( size_t i=0; i<size; ++i )
          mass_[ i ].setQuadratureOrderFunctions( interiorOrder, surfaceOrder );
      }

      void setCommunicate( const bool communicate )
      {
        communicate_ = communicate;
        if( ! communicate_ && Dune::Fem::Parameter::verbose() )
        {
          std::cout << "MassLumpingOperator::setCommunicate: communicate was disabled!" << std::endl;
        }
      }

      void setQuadratureOrders(unsigned int interior, unsigned int surface)
      {
        // only set quadrature orders for Galerkin part, lumping quadratures are fixed
        size_t size = localOp_.size();
        for( size_t i=0; i<size; ++i )
          localOp_[ i ].setQuadratureOrders(interior,surface);
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        evaluate( u, w );
      }

      const GridPartType &gridPart () const { return op().gridPart(); }

      typedef Integrands ModelType;
      typedef Integrands DirichletModelType;
      ModelType &model() const { return localOperator().model(); }

      [[deprecated("Use localOperator instead!")]]
      const LocalGalerkinOperatorImplType& impl() const { return localOperator(); }

      //! return local operator holding instance of integrands
      const LocalGalerkinOperatorImplType& localOperator() const { return (*localOp_); }
      const MassOperatorImplType& mass() const { return (*mass_); }

      std::size_t gridSizeInterior () const { return gridSizeInterior_; }

    protected:
      const GalerkinOperatorImplType& op() const { return (*opImpl_); }

      template< class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        iterators_.update();
        w.clear();

        std::shared_mutex mutex;

        auto doEval = [this, &u, &w, &mutex] ()
        {
          // version with locking
          std::tuple< const LocalGalerkinOperatorImplType&,
                      const MassOperatorImplType& > integrands( localOperator(), mass() );


          // add galerkin and mass lumped part
          this->op().evaluate( u, w, this->iterators_, integrands, mutex );
        };

        bool singleThreadModeError = false ;

        try {
          // execute in parallel
          MPIManager :: run ( doEval );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = Impl::accumulateGridSize( opImpl_ );
        }
        catch ( const SingleThreadModeError& e )
        {
          singleThreadModeError = true;
        }

        // if error occurred, redo the whole evaluation
        if( singleThreadModeError )
        {
          // reset w from previous entries
          w.clear();

          // re-run in single thread mode if previous attempt failed
          std::tuple< const LocalGalerkinOperatorImplType&,
                      const MassOperatorImplType& > integrands( localOperator(), mass() );
          // add galerkin part and mass part
          op().evaluate( u, w, iterators_, integrands );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = op().gridSizeInterior();

        }

        // synchronize data
        if( communicate_ )
          w.communicate();
      }

      // GalerkinOperator implementation (see galerkin.hh)
      mutable ThreadIteratorType iterators_;
      ThreadSafeValue< GalerkinOperatorImplType > opImpl_;
      ThreadSafeValue< LocalGalerkinOperatorImplType > localOp_;
      ThreadSafeValue< MassOperatorImplType     > mass_;

      mutable std::size_t gridSizeInterior_;
      bool communicate_;
    };



    // DifferentiableMassLumpingOperator
    // ---------------------------------

    template< class Integrands, class MassIntegrands, class JacobianOperator >
    class MassLumpingDifferentiableOperator
      : public MassLumpingOperator< Integrands, MassIntegrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType >,
        public DifferentiableOperator< JacobianOperator >
    {
     typedef MassLumpingOperator< Integrands, MassIntegrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType > BaseType;

    public:
      typedef typename BaseType :: LocalGalerkinOperatorImplType  LocalGalerkinOperatorImplType;
      typedef typename BaseType :: MassOperatorImplType           MassOperatorImplType;

      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;

      typedef DiagonalAndNeighborStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType >  DiagonalAndNeighborStencilType;
      typedef DiagonalStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType >             DiagonalStencilType;

      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit MassLumpingDifferentiableOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                   const RangeDiscreteFunctionSpaceType &rSpace,
                                                   Args &&... args )
        : BaseType( rSpace.gridPart(), std::forward< Args >( args )... ),
          dSpace_(dSpace),
          rSpace_(rSpace),
          domainSpaceSequence_(dSpace.sequence()),
          rangeSpaceSequence_(rSpace.sequence()),
          stencilDAN_(), stencilD_()
      {
        if( hasSkeleton() )
          stencilDAN_.reset( new DiagonalAndNeighborStencilType( dSpace_, rSpace_ ) );
        else
          stencilD_.reset( new DiagonalStencilType( dSpace_, rSpace_ ) );
      }

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        // assemble Jacobian, same as GalerkinOperator
        assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        assemble( u, jOp );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return rSpace_;
      }
      // needed for DGHelmholtz operator
      const RangeDiscreteFunctionSpaceType& space() const
      {
        return rangeSpace();
      }

      using BaseType::gridPart;
      using BaseType::localOperator;

    protected:
      bool hasSkeleton() const
      {
        std::tuple< const LocalGalerkinOperatorImplType&, const MassOperatorImplType&  > integrands( localOperator(), mass() );
        return op().hasSkeleton( integrands );
      }

      using BaseType::op;
      using BaseType::mass;
      using BaseType::iterators_;
      using BaseType::gridSizeInterior_;

      void prepare( JacobianOperatorType& jOp ) const
      {
        if ( domainSpaceSequence_ != domainSpace().sequence()
             || rangeSpaceSequence_ != rangeSpace().sequence() )
        {
          domainSpaceSequence_ = domainSpace().sequence();
          rangeSpaceSequence_ = rangeSpace().sequence();
          if( hasSkeleton() )
          {
            assert( stencilDAN_ );
            stencilDAN_->update();
          }
          else
          {
            assert( stencilD_ );
            stencilD_->update();
          }
        }

        if( hasSkeleton() )
          jOp.reserve( *stencilDAN_ );
        else
          jOp.reserve( *stencilD_ );
        // set all entries to zero
        jOp.clear();
      }

      template < class GridFunction >
      void assemble( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        prepare( jOp );
        iterators_.update();

        bool singleThreadModeError = false;
        std::shared_mutex mutex;

        auto doAssemble = [this, &u, &jOp, &mutex] ()
        {
          // assemble Jacobian, same as GalerkinOperator
          std::tuple< const LocalGalerkinOperatorImplType&, const MassOperatorImplType&  > integrands( localOperator(), mass() );
          this->op().assemble( u, jOp, this->iterators_, integrands, mutex );
        };

        try {
          // execute in parallel
          MPIManager :: run ( doAssemble );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = Impl::accumulateGridSize( this->opImpl_ );
        }
        catch ( const SingleThreadModeError& e )
        {
          singleThreadModeError = true;
        }

        if (singleThreadModeError)
        {
          // redo matrix assembly since it failed
          jOp.clear();
          // add galerkin terms
          std::tuple< const LocalGalerkinOperatorImplType&, const MassOperatorImplType&  > integrands( localOperator(), mass() );
          op().assemble( u, jOp, iterators_, integrands );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = op().gridSizeInterior();
        }

        // note: assembly done without local contributions so need
        // to call flush assembly
        jOp.flushAssembly();
      }

      const DomainDiscreteFunctionSpaceType &dSpace_;
      const RangeDiscreteFunctionSpaceType &rSpace_;
      mutable int domainSpaceSequence_, rangeSpaceSequence_;

      mutable std::unique_ptr< DiagonalAndNeighborStencilType > stencilDAN_;
      mutable std::unique_ptr< DiagonalStencilType >            stencilD_;
    };



    // AutomaticDifferenceMassLumpingOperator
    // --------------------------------------

    template< class Integrands, class DomainFunction, class RangeFunction >
    class MassLumpingAutomaticDifferenceGalerkinOperator
      : public MassLumpingOperator< Integrands, DomainFunction, RangeFunction >,
        public AutomaticDifferenceOperator< DomainFunction, RangeFunction >
    {
      typedef MassLumpingOperator< Integrands, DomainFunction, RangeFunction > BaseType;
      typedef AutomaticDifferenceOperator< DomainFunction, RangeFunction > AutomaticDifferenceOperatorType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit MassLumpingAutomaticDifferenceGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : BaseType( gridPart, std::forward< Args >( args )... ), AutomaticDifferenceOperatorType()
      {}
    };



    namespace Impl
    {

      // MassLumpingSchemeImpl
      // ---------------------
      template< class Integrands, class MassIntegrands,
                class LinearOperator, bool addDirichletBC,
                template <class,class,class> class DifferentiableGalerkinOperatorImpl >
      struct MassLumpingSchemeTraits
      {
        template <class O, bool addDBC>
        struct DirichletBlockSelector { using type = void; };
        template <class O>
        struct DirichletBlockSelector<O,true> { using type = typename O::DirichletBlockVector; };

        using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
           DirichletWrapperOperator< DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator >>,
           DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator > >;

        using DirichletBlockVector = typename DirichletBlockSelector<
                 DirichletWrapperOperator<
                    DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator >>,
                 addDirichletBC>::type;

        typedef DifferentiableOperatorType type;
      };


      // MassLumpingSchemeImpl
      // ---------------------
      template< class Integrands, class MassIntegrands,
                class LinearOperator, class LinearInverseOperator, bool addDirichletBC,
                template <class,class,class> class DifferentiableGalerkinOperatorImpl = MassLumpingDifferentiableOperator >
      struct MassLumpingSchemeImpl : public FemScheme<
                     typename MassLumpingSchemeTraits< Integrands, MassIntegrands, LinearOperator,
                                                       addDirichletBC, DifferentiableGalerkinOperatorImpl
                                                     > :: type, // Operator
                                      LinearInverseOperator > // LinearInverseOperator
      {
        typedef FemScheme< typename MassLumpingSchemeTraits< Integrands, MassIntegrands, LinearOperator,
                                   addDirichletBC, DifferentiableGalerkinOperatorImpl>::type, // Operator
                                   LinearInverseOperator > // LinearInverseOperator
                        BaseType;

        // type of discrete function space used with this scheme
        typedef typename BaseType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;

        // needed for registerSchemeConstructor (Python export)
        typedef MassIntegrands  MassModelType;

        MassLumpingSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
                                const Integrands &integrands,
                                const MassIntegrands& massIntegrands,
                                const ParameterReader& parameter = Parameter::container() )
          : BaseType( dfSpace,
                      parameter,
                      std::move(integrands),
                      std::move(massIntegrands) )
        {}
      };

    } // end namespace Impl

    // MassLumpingScheme
    // -----------------

    template< class Integrands, class MassIntegrands, class LinearOperator, class InverseOperator, bool addDirichletBC >
    using MassLumpingScheme = Impl::MassLumpingSchemeImpl< Integrands, MassIntegrands, LinearOperator, InverseOperator, addDirichletBC,
                                                           MassLumpingDifferentiableOperator >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_MOLGALERKIN_HH
