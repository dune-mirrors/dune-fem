#ifndef DUNE_ASH_SOLVER_ISTL_HH
#define DUNE_ASH_SOLVER_ISTL_HH

#include <cassert>

#include <limits>
#include <map>
#include <memory>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/remoteindices.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/partitionset.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/mapper/parallel.hh>
#include <dune/fem/solver/parameter.hh>


namespace Dune
{

  namespace Amg
  {

    template< class M, class X, class Y >
    struct ConstructionTraits< SeqILDL< M, X, Y > >
    {
      typedef DefaultConstructionArgs< SeqILDL< M, X, Y > > Arguments;

      static SeqILDL< M, X, Y > *construct ( Arguments &args )
      {
        return new SeqILDL< M, X, Y >( args.getMatrix(), args.getArgs().relaxationFactor );
      }

      static void deconstruct ( SeqILDL< M, X, Y > *p ) { delete p; }
    };

  } // namespace Amg



  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class DiscreteFunctionSpace >
    class HierarchicalDiscreteFunction;

    template< class DomainFunction, class RangeFunction >
    class ISTLLinearOperator;

  } // namespace Fem



  namespace Fem
  {

    namespace ISTL
    {

      // VectorType
      // ----------

      template< class DiscreteFunction >
      using VectorType = std::decay_t< decltype( std::declval< const DiscreteFunction & >().dofVector().array() ) >;



      // MatrixType
      // ----------

      template< class LinearOperator >
      struct __MatrixType
      {
        typedef std::decay_t< decltype( std::declval< const LinearOperator & >().matrix() ) > Type;
      };

      template< class DomainFunction, class RangeFunction >
      struct __MatrixType< Dune::Fem::ISTLLinearOperator< DomainFunction, RangeFunction > >
      {
        typedef Dune::BCRSMatrix< typename Dune::Fem::ISTLLinearOperator< DomainFunction, RangeFunction >::LittleBlockType > Type;
      };


      template< class LinearOperator >
      using MatrixType = typename __MatrixType< LinearOperator >::Type;



      // Symmetry
      // --------

      enum Symmetry : bool { unsymmetric = false, symmetric = true };



      // FemCommunication
      // ----------------

      template< class DiscreteFunctionSpace >
      class FemCommunication
      {
        typedef FemCommunication< DiscreteFunctionSpace > ThisType;

      public:
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        typedef int GlobalLookupIndexSet;

        explicit FemCommunication ( const DiscreteFunctionSpaceType &dfSpace,  Dune::SolverCategory::Category solverCategory = Dune::SolverCategory::sequential )
          : dfSpace_( dfSpace ), solverCategory_( solverCategory )
        {}

        const typename DiscreteFunctionSpace::GridPartType::CommunicationType &communicator () const { return dfSpace_.gridPart().comm(); }

        template< class T >
        void copyOwnerToAll ( const T &x, T &y ) const
        {
          y = x;
          typedef Dune::Fem::HierarchicalDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
          DiscreteFunctionType z( "", dfSpace_, reinterpret_cast< typename DiscreteFunctionType::DofVectorType & >( y ) );
          z.communicate();
        }

        template< class T >
        void project ( T &x ) const
        {
          typedef typename T::field_type field_type;

          // clear auxiliary DoFs
          const auto &auxiliaryDofs = dfSpace_.auxiliaryDofs();
          for( int i : auxiliaryDofs )
            x[ i ] = field_type( 0 );
        }

        template< class T, class F >
        void dot ( const T &x, const T &y, F &scp ) const
        {
          const auto &auxiliaryDofs = dfSpace_.auxiliaryDofs();

          const int numAuxiliarys = auxiliaryDofs.size();
          for( int auxiliary = 0, i = 0; auxiliary < numAuxiliarys; ++auxiliary, ++i )
          {
            const int nextAuxiliary = auxiliaryDofs[ auxiliary ];
            for( ; i < nextAuxiliary; ++i )
              scp += x[ i ] * y[ i ];
          }

          scp = communicator().sum( scp );
        }

        template< class T >
        typename Dune::FieldTraits< typename T::field_type >::real_type norm ( const T &x ) const
        {
          using std::sqrt;
          typename Dune::FieldTraits< typename T::field_type >::real_type norm2( 0 );
          dot( x, x, norm2 );
          return sqrt( norm2 );
        }

        Dune::SolverCategory::Category getSolverCategory () const { return solverCategory_; }

      private:
        const DiscreteFunctionSpaceType &dfSpace_;
        Dune::SolverCategory::Category solverCategory_;
      };



      // BuildRemoteIndicesDataHandle
      // ----------------------------

      template< class Mapper, class GlobalLookup >
      struct BuildRemoteIndicesDataHandle
        : public Dune::CommDataHandleIF< BuildRemoteIndicesDataHandle< Mapper, GlobalLookup >, int >
      {
        typedef typename GlobalLookup::GlobalIndex GlobalIndexType;
        typedef typename GlobalLookup::LocalIndex::Attribute AttributeType;

        BuildRemoteIndicesDataHandle ( int rank, const Mapper &mapper, const GlobalLookup &globalLookup )
          : rank_( rank ), mapper_( mapper ), globalLookup_( globalLookup )
        {}

        bool contains ( int dim, int codim ) const { return mapper_.contains( codim ); }
        bool fixedSize( int dim, int codim ) const { return true; }

        template< class Buffer, class Entity >
        void gather ( Buffer &buffer, const Entity &entity ) const
        {
          buffer.write( rank_ );
          int attribute = -1;
          mapper_.mapEachEntityDof( entity, [ this, &attribute ] ( int, auto index ) {
              auto *pair = globalLookup_.pair( index );
              assert( pair && ((attribute == -1) || (attribute == pair->local().attribute())) );
              attribute = pair->local().attribute();
            } );
          buffer.write( static_cast< int >( attribute ) );
        }

        template< class Buffer, class Entity >
        void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
        {
          int rank, attribute;
          buffer.read( rank );
          buffer.read( attribute );
          assert( (attribute != -1) || (mapper_.numEntityDofs( entity ) == 0) );
          mapper_.mapEachEntityDof( entity, [ this, rank, attribute ] ( int, auto index ) {
              auto *pair = globalLookup_.pair( index );
              assert( pair );
              remotes[ rank ].emplace_back( static_cast< AttributeType >( attribute ), pair );
            } );
        }

        template< class Entity >
        std::size_t size ( const Entity &entity ) const
        {
          return 2;
        }

        std::map< int, std::vector< Dune::RemoteIndex< GlobalIndexType, AttributeType > > > remotes;

      private:
        int rank_;
        const Mapper &mapper_;
        const GlobalLookup &globalLookup_;
      };



      // buildCommunication
      // ------------------

      template< class DiscreteFunction >
      void buildCommunication ( const typename DiscreteFunction::DiscreteFunctionSpaceType &dfSpace,
                                Dune::SolverCategory::Category solverCategory,
                                std::shared_ptr< FemCommunication< DiscreteFunction > > &communication )
      {
        communication.reset( new FemCommunication< DiscreteFunction >( dfSpace, solverCategory ) );
      }

      template< class DiscreteFunctionSpace, class GlobalId, class LocalId >
      void buildCommunication ( const DiscreteFunctionSpace &dfSpace,
                                Dune::SolverCategory::Category solverCategory,
                                std::shared_ptr< Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId > > &communication )
      {
        typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
        typedef typename DiscreteFunctionSpace::BlockMapperType LocalMapperType;

        typedef typename Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId >::GlobalLookupIndexSet GlobalLookupType;

        typedef typename GlobalLookupType::LocalIndex LocalIndexType;

        communication.reset( new Dune::OwnerOverlapCopyCommunication< GlobalId, LocalId >( solverCategory ) );

        const GridPartType &gridPart = dfSpace.gridPart();
        LocalMapperType &localMapper = dfSpace.blockMapper();

        // create global index mapping
        Dune::Fem::ParallelDofMapper< GridPartType, LocalMapperType, GlobalId > globalMapper( gridPart, localMapper );

        // construct local attributes
        std::vector< typename LocalIndexType::Attribute > attribute( localMapper.size(), Dune::OwnerOverlapCopyAttributeSet::owner );
        for( const auto &auxiliary : dfSpace.auxiliaryDofs() )
          attribute[ auxiliary ] = Dune::OwnerOverlapCopyAttributeSet::copy;

        // build parallel index set
        communication->indexSet().beginResize();
        for( LocalId i = 0, size = localMapper.size(); i < size; ++i )
          communication->indexSet().add( globalMapper.mapping()[ i ], LocalIndexType( i, attribute[ i ] ) );
        communication->indexSet().endResize();

        // build remote indices
        communication->buildGlobalLookup();
        BuildRemoteIndicesDataHandle< LocalMapperType, GlobalLookupType > buildRemoteIndicesDataHandle( gridPart.comm().rank(), localMapper, communication->globalLookup() );
        gridPart.communicate( buildRemoteIndicesDataHandle, Dune::All_All_Interface, Dune::ForwardCommunication );
        communication->freeGlobalLookup();

        communication->remoteIndices().setIndexSets( communication->indexSet(), communication->indexSet(), communication->communicator() );
        if( !buildRemoteIndicesDataHandle.remotes.empty() )
        {
          for( auto &remote : buildRemoteIndicesDataHandle.remotes )
          {
            std::sort( remote.second.begin(), remote.second.end(), [] ( const auto &a, const auto &b ) { return (a.localIndexPair().global() < b.localIndexPair().global()); } );
            auto modifier = communication->remoteIndices().template getModifier< false, true >( remote.first );
            for( const auto &remoteIndex : remote.second )
              modifier.insert( remoteIndex );
          }
        }
        else
          communication->remoteIndices().template getModifier< false, true >( 0 );
      }



      // NamedType
      // ---------

      template< class T >
      struct NamedType
      {
        typedef T Type;

        NamedType ( std::string name ) : name( name ) {}

        std::string name;
      };



      // getEnum
      // -------

      template< class... T, class F >
      void getEnum ( const Dune::Fem::ParameterReader &parameter, const std::string &key, const std::tuple< NamedType< T >... > &types, const std::string &defaultValue, F &&f )
      {
        const std::string &value = parameter.getValue( key, defaultValue );
        bool success = false;
        Dune::Hybrid::forEach( std::index_sequence_for< T... >(), [ &types, &f, &value, &success ] ( auto &&i ) -> void {
            const auto &type = std::get< std::decay_t< decltype( i ) >::value >( types );
            if( value != type.name )
              return;

            assert( !success );
            success = true;

            f( type );
          } );
        if( !success )
          DUNE_THROW( Dune::Fem::ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      template< class... T, class F >
      void getEnum ( const std::string &key, const std::tuple< NamedType< T >... > &types, const std::string &defaultValue, F &&f )
      {
        getEnum( Dune::Fem::Parameter::container(), key, types, defaultValue, std::forward< F >( f ) );
      }



      // makePreconditioner
      // ------------------

      template< class AssembledOperator, class Communication >
      inline std::shared_ptr< Dune::Preconditioner< typename AssembledOperator::domain_type, typename AssembledOperator::range_type > >
      makePreconditioner ( const std::shared_ptr< AssembledOperator > op, const Communication &comm, Symmetry symmetry )
      {
        typedef typename AssembledOperator::matrix_type matrix_type;
        typedef typename AssembledOperator::domain_type domain_type;
        typedef typename AssembledOperator::range_type range_type;
        typedef typename Dune::FieldTraits< typename AssembledOperator::field_type >::real_type real_type;

        Dune::Amg::DefaultSmootherArgs< real_type > smootherArgs;
        smootherArgs.relaxationFactor = Dune::Fem::Parameter::getValue( "istl.preconditioner.relax", real_type( 1 ) );

        const auto smootherTypes = std::make_tuple( NamedType< Dune::SeqJac < matrix_type, domain_type, range_type > >( "jacobi" ),
                                                    NamedType< Dune::SeqSSOR< matrix_type, domain_type, range_type > >( "ssor" ),
                                                    NamedType< Dune::SeqILU < matrix_type, domain_type, range_type > >( "ilu" ),
                                                    NamedType< Dune::SeqILDL< matrix_type, domain_type, range_type > >( "ildl" ) );

        std::shared_ptr< Dune::Preconditioner< domain_type, range_type > > preconditioner;
        const std::string preconditionerTypes[] = { "richardson", "smoother", "amg" };
        switch( Dune::Fem::Parameter::getEnum( "istl.preconditioner.type", preconditionerTypes, 1 ) )
        {
        case 0:
          {
            auto sp = std::make_shared< Dune::Richardson< domain_type, range_type > >( smootherArgs.relaxationFactor );
            auto *bp = new Dune::BlockPreconditioner< domain_type, range_type, Communication, std::decay_t< decltype( *sp ) > >( *sp, comm );
            preconditioner.reset( bp, [ op, sp ] ( decltype( bp ) p ) { delete p; } );
          }
          break;

        case 1:
          getEnum( "istl.preconditioner.smoother", smootherTypes, "jacobi", [ op, &comm, &smootherArgs, &preconditioner ] ( auto namedType ) {
              typedef Dune::BlockPreconditioner< domain_type, range_type, Communication, typename decltype( namedType )::Type > SmootherType;
              typedef Dune::Amg::ConstructionTraits< SmootherType > ConstructionTraits;

              typename ConstructionTraits::Arguments args;
              args.setArgs( smootherArgs );
              args.setComm( comm );
              args.setMatrix( op->getmat() );

              preconditioner.reset( ConstructionTraits::construct( args ), [ op ] ( SmootherType *p ) { ConstructionTraits::deconstruct( p ); } );
            } );
          break;

        case 2:
          getEnum( "istl.preconditioner.smoother", smootherTypes, "jacobi", [ op, &comm, symmetry, &smootherArgs, &preconditioner ] ( auto namedType ) {
              typedef Dune::BlockPreconditioner< domain_type, range_type, Communication, typename decltype( namedType )::Type > SmootherType;

              Dune::Amg::Parameters amgParams;

              // coarsening parameters
              amgParams.setMaxLevel( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.maxlevel", 100 ) );
              amgParams.setCoarsenTarget( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.coarsentarget", 1000 ) );
              amgParams.setMinCoarsenRate( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.mincoarsenrate", 1.2 ) );
              const std::string accumulationModeNames[] = { "none", "once", "successive" };
              amgParams.setAccumulate( static_cast< Dune::Amg::AccumulationMode >( Dune::Fem::Parameter::getEnum( "istl.preconditioner.amg.accumulate", accumulationModeNames, 2 ) ) );
              amgParams.setProlongationDampingFactor( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.prolongation.dampingfactor", 1.6 ) );

              // parameters
              amgParams.setDebugLevel( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.debuglevel", 0 ) );
              amgParams.setNoPreSmoothSteps( Dune::Fem::Parameter::getValue< std::size_t >( "istl.preconditioner.amg.presmoothsteps", 2 ) );
              amgParams.setNoPostSmoothSteps( Dune::Fem::Parameter::getValue< std::size_t >( "istl.preconditioner.amg.postsmoothsteps", 2 ) );
              const std::string cycleNames[] = { "v-cycle", "w-cycle" };
              amgParams.setGamma( 1 + Dune::Fem::Parameter::getEnum( "istl.preconditioner.amg.cycle", cycleNames, 0 ) );
              amgParams.setAdditive( Dune::Fem::Parameter::getValue( "istl.preconditioner.amg.additive", false ) );

              typedef Dune::Amg::AMG< AssembledOperator, domain_type, SmootherType, Communication > AMG;
              if( symmetry == symmetric )
              {
                Dune::Amg::CoarsenCriterion< Dune::Amg::SymmetricCriterion< matrix_type, Dune::Amg::RowSum > > criterion( amgParams );
                preconditioner.reset( new AMG( *op, criterion, smootherArgs, comm ), [ op ] ( AMG *p ) { delete p; } );
              }
              else
              {
                Dune::Amg::CoarsenCriterion< Dune::Amg::UnSymmetricCriterion< matrix_type, Dune::Amg::RowSum > > criterion( amgParams );
                preconditioner.reset( new AMG( *op, criterion, smootherArgs, comm ), [ op ] ( AMG *p ) { delete p; } );
              }
            } );
          break;

        default:
          DUNE_THROW( Dune::InvalidStateException, "Invalid ISTL preconditioner type selected." );
        }
        return preconditioner;
      }



      // InverseOperator
      // ---------------

      template< class LinearOperator, Symmetry symmetry = unsymmetric >
      class InverseOperator final
        : public Dune::Fem::Operator< typename LinearOperator::RangeFunctionType, typename LinearOperator::DomainFunctionType >
      {
        static_assert( std::is_same< typename LinearOperator::DomainFunctionType, typename LinearOperator::RangeFunctionType >::value, "Domain function and range function must have the same type." );

      public:
        typedef LinearOperator LinearOperatorType;
        typedef LinearOperatorType  OperatorType;

        typedef typename LinearOperatorType::DomainFunctionType DiscreteFunctionType;

        typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type RealType;

        typedef SolverParameter SolverParameterType;

      private:
        typedef VectorType< DiscreteFunctionType > vector_type;
        typedef MatrixType< LinearOperatorType > matrix_type;

      public:
        // typedef FemCommunication< DiscreteFunctionSpaceType > CommunicationType;
        typedef Dune::OwnerOverlapCopyCommunication< std::size_t, typename DiscreteFunctionSpaceType::BlockMapperType::GlobalKeyType > CommunicationType;

        typedef Dune::OverlappingSchwarzOperator< matrix_type, vector_type, vector_type, CommunicationType > AssembledLinearOperatorType;
        typedef Dune::Preconditioner< vector_type, vector_type > PreconditionerType;
        typedef Dune::ScalarProduct< vector_type > ScalarProductType;

        typedef std::function< std::shared_ptr< PreconditionerType > ( std::shared_ptr< AssembledLinearOperatorType >, const CommunicationType &, Symmetry ) > PreconditionerFactory;

        InverseOperator ( PreconditionerFactory preconditionerFactory, RealType redEps, RealType absLimit, int maxIterations )
          : preconditionerFactory_( std::move( preconditionerFactory ) ),
            redEps_( redEps ), absLimit_( absLimit ), maxIterations_( maxIterations )
        {}

        InverseOperator ( RealType redEps, RealType absLimit, int maxIterations )
          : InverseOperator( makePreconditioner< AssembledLinearOperatorType, CommunicationType >, redEps, absLimit, maxIterations )
        {}

        InverseOperator ( RealType redEps, RealType absLimit, int maxIterations, bool verbose, const Dune::Fem::SolverParameter& parameter )
          : InverseOperator( redEps, absLimit, maxIterations )
        {}

        InverseOperator ( PreconditionerFactory preconditionerFactory, RealType redEps, RealType absLimit )
          : InverseOperator( std::move( preconditionerFactory ), redEps, absLimit, std::numeric_limits< int >::max() )
        {}

        InverseOperator ( RealType redEps, RealType absLimit )
          : InverseOperator( redEps, absLimit, std::numeric_limits< int >::max() )
        {}

        InverseOperator ( LinearOperator &op, RealType redEps, RealType absLimit, int maxIterations )
          : InverseOperator( redEps, absLimit, maxIterations )
        {
          bind( op );
        }

        InverseOperator ( LinearOperator &op, RealType redEps, RealType absLimit )
          : InverseOperator( op, redEps, absLimit, std::numeric_limits< int >::max() )
        {}

        InverseOperator ( const SolverParameter& parameter = SolverParameter( Parameter::container() ) )
          : InverseOperator( parameter.tolerance(), parameter.tolerance(), parameter.maxIterations() )
        {}

        void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const override
        {
          Dune::InverseOperatorResult result;
          vector_type b( u.dofVector().array() );
          if( absLimit_ < std::numeric_limits< RealType >::max() )
          {
            vector_type tmp( b );
            linearOperator_->apply( w.dofVector().array(), tmp );
            tmp -= b;
            RealType res = scalarProduct_->norm( tmp );
            RealType red = (res >0)? absLimit_ / res : 1e-3;
            solver_->apply( w.dofVector().array(), b, red, result );
          }
          else
            solver_->apply( w.dofVector().array(), b, result );
          iterations_ = result.iterations;
        }

        void bind ( LinearOperator &op )
        {
          buildCommunication( op.domainSpace(), Dune::SolverCategory::overlapping, communication_ );

          linearOperator_.reset( new AssembledLinearOperatorType( op.matrix(), *communication_ ) );
          scalarProduct_.reset( new Dune::OverlappingSchwarzScalarProduct< vector_type, CommunicationType >( *communication_ ) );

          // create preconditioner
          preconditioner_ = preconditionerFactory_( linearOperator_, *communication_, symmetry );

          // choose absolute error or reduction error
          const std::string reductionType[] = { "absolute", "relative" };
          int errorType = Dune::Fem::Parameter::getEnum( "istl.solver.errormeasure", reductionType, 0 );
          if( errorType != 0 )
            absLimit_ = std::numeric_limits< RealType >::max();

          // create linear solver
          const std::string verbosityTable[] = { "off", "on", "full" };
          int verbosity = Dune::Fem::Parameter::getEnum( "istl.solver.verbosity", verbosityTable, 0 );
          if( op.domainSpace().gridPart().comm().rank() != 0 )
            verbosity = 0;

          const std::string solverTypes[] = { "cg", "gcg", "minres", "bicgstab", "gmres" };
          switch( Dune::Fem::Parameter::getEnum( "istl.solver.type", solverTypes, (symmetry == symmetric ? 0 : 3) ) )
          {
          case 0:
            solver_.reset( new Dune::CGSolver< vector_type >( *linearOperator_, *scalarProduct_, *preconditioner_, redEps_, maxIterations_, verbosity ) );
            break;

          case 1:
            {
              const int restart = Dune::Fem::Parameter::getValue( "istl.solver.restart", 20 );
              solver_.reset( new Dune::GeneralizedPCGSolver< vector_type >( *linearOperator_, *scalarProduct_, *preconditioner_, redEps_, maxIterations_, verbosity, restart ) );
            }
            break;

          case 2:
            solver_.reset( new Dune::MINRESSolver< vector_type >( *linearOperator_, *scalarProduct_, *preconditioner_, redEps_, maxIterations_, verbosity ) );
            break;

          case 3:
            solver_.reset( new Dune::BiCGSTABSolver< vector_type >( *linearOperator_, *scalarProduct_, *preconditioner_, redEps_, maxIterations_, verbosity ) );
            break;

          case 4:
            {
              const int restart = Dune::Fem::Parameter::getValue( "istl.solver.restart", 20 );
              solver_.reset( new Dune::RestartedGMResSolver< vector_type >( *linearOperator_, *scalarProduct_, *preconditioner_, redEps_, restart, maxIterations_, verbosity ) );
            }
            break;
          }
        }

        void unbind ()
        {
          solver_.reset();
          preconditioner_.reset();
          scalarProduct_.reset();
          linearOperator_.reset();
          communication_.reset();
        }

        int iterations () const { return iterations_; }

        void setMaxIterations ( int maxIterations ) { maxIterations_ = maxIterations; }

      private:
        PreconditionerFactory preconditionerFactory_;
        RealType redEps_, absLimit_;
        int maxIterations_;

        std::shared_ptr< CommunicationType > communication_;
        std::shared_ptr< AssembledLinearOperatorType > linearOperator_;
        std::shared_ptr< ScalarProductType > scalarProduct_;
        std::shared_ptr< PreconditionerType > preconditioner_;
        std::shared_ptr< Dune::InverseOperator< vector_type, vector_type > > solver_;

        mutable int iterations_;
      };

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_ISTL

#endif // #ifndef DUNE_ASH_SOLVER_ISTL_HH
