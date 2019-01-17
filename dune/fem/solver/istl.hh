#ifndef DUNE_FEM_SOLVER_ISTL_HH
#define DUNE_FEM_SOLVER_ISTL_HH

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
#include <dune/common/typelist.hh>
#include <dune/common/typeutilities.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/partitionset.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/mapper/parallel.hh>

#include <dune/fem/solver/communication/fem.hh>
#include <dune/fem/solver/communication/hierarchical.hh>
#include <dune/fem/solver/communication/owneroverlapcopy.hh>


namespace Dune
{

  namespace Amg
  {

    // ConstructionTraits for SeqGS
    // ----------------------------

    template< class M, class X, class Y, int l >
    struct ConstructionTraits< SeqGS< M, X, Y, l > >
    {
      typedef DefaultConstructionArgs< SeqGS< M, X, Y, l > > Arguments;

      static SeqGS< M, X, Y, l > *construct ( Arguments &args )
      {
        return new SeqGS< M, X, Y, l >( args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor );
      }

      static void deconstruct ( SeqGS< M, X, Y, l > *p ) { delete p; }
    };



    // ConstructionTraits for SeqILDL
    // ------------------------------

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



      // IsBCRSMatrix
      // ------------

      template< class M >
      struct IsBCRSMatrix
        : std::false_type
      {};

      template< class B, class A >
      struct IsBCRSMatrix< Dune::BCRSMatrix< B, A > >
        : std::true_type
      {};



      // Symmetry
      // --------

      enum Symmetry : bool { unsymmetric = false, symmetric = true };



      // getEnum
      // -------

      template< class... T, class F >
      void getEnum ( const Dune::Fem::ParameterReader &parameter, const std::string &key, const std::tuple< std::pair< std::string, T >... > &choices, const std::string &defaultValue, F &&f )
      {
        const std::string &value = parameter.getValue( key, defaultValue );
        bool success = false;
        Dune::Hybrid::forEach( choices, [ &choices, &f, &value, &success ] ( const auto &choice ) {
            if( value != choice.first )
              return;

            assert( !success );
            success = true;

            f( choice.second );
          } );
        if( success )
          return;

        std::string values;
        Dune::Hybrid::forEach( choices, [ &choices, &values ] ( const auto &choice ) { values += (values.empty() ? "'" : ", '") + choice.first + "'"; } );
        DUNE_THROW( Dune::Fem::ParameterInvalid, "Parameter '" << key << "' invalid (choices are: " << values << ")." );
      }

      template< class... T, class F >
      void getEnum ( const std::string &key, const std::tuple< std::pair< std::string, T >... > &choices, const std::string &defaultValue, F &&f )
      {
        getEnum( Dune::Fem::Parameter::container(), key, choices, defaultValue, std::forward< F >( f ) );
      }



      // UniformLeafLevel
      // ----------------

      template< class T, class = void >
      struct UniformLeafLevelType;

      template< class T >
      using UniformLeafLevel = typename UniformLeafLevelType< T >::Type;

      template< class V1, class... V >
      struct UniformLeafLevelType< Dune::MultiTypeBlockVector< V1, V... >, std::enable_if_t< Std::are_all_same< UniformLeafLevel< V1 >, UniformLeafLevel< V >... >::value > >
      {
        typedef std::integral_constant< int, UniformLeafLevel< V1 >::value + 1 > Type;
      };

      template< class B, class A >
      struct UniformLeafLevelType< Dune::BlockVector< B, A > >
      {
        typedef std::integral_constant< int, UniformLeafLevel< B >::value + 1 > Type;
      };

      template< class K, int n >
      struct UniformLeafLevelType< Dune::FieldVector< K, n > >
      {
        typedef std::integral_constant< int, 0 > Type;
      };

      template< class R1, class... R >
      struct UniformLeafLevelType< Dune::MultiTypeBlockMatrix< R1, R... >, std::enable_if_t< Std::are_all_same< UniformLeafLevel< R1 >, UniformLeafLevel< R >... >::value > >
      {
        // Note: The rows of MultiTypeBlockMatrix are of type MultiTypeBlockVector, which already increases the level
        typedef std::integral_constant< int, UniformLeafLevel< R1 >::value > Type;
      };

      template< class B, class A >
      struct UniformLeafLevelType< Dune::BCRSMatrix< B, A > >
      {
        typedef std::integral_constant< int, UniformLeafLevel< B >::value + 1 > Type;
      };

      template< class K, int m, int n >
      struct UniformLeafLevelType< Dune::FieldMatrix< K, m, n > >
      {
        typedef std::integral_constant< int, 0 > Type;
      };




      // namedSmootherTypes
      // ------------------

      template< class M, class X, class Y, std::enable_if_t< IsBCRSMatrix< M >::value && (UniformLeafLevel< M >::value == 1), int > = 0 >
      inline static decltype( auto ) namedSmootherTypes ( PriorityTag< 2 > )
      {
        return std::make_tuple( std::make_pair( std::string( "jacobi" ), Dune::MetaType< Dune::SeqJac< M, X, Y, 1 > >() ),
                                std::make_pair( std::string( "gauss-seidel" ), Dune::MetaType< Dune::SeqGS< M, X, Y, 1 > >() ),
                                std::make_pair( std::string( "sor" ), Dune::MetaType< Dune::SeqSOR< M, X, Y, 1 > >() ),
                                std::make_pair( std::string( "ssor" ), Dune::MetaType< Dune::SeqSSOR< M, X, Y, 1 > >() ),
                                std::make_pair( std::string( "ilu0" ), Dune::MetaType< Dune::SeqILU0< M, X, Y > >() ),
                                std::make_pair( std::string( "ildl" ), Dune::MetaType< Dune::SeqILDL< M, X, Y > >() ) );
      }

      template< class M, class X, class Y, std::enable_if_t< (UniformLeafLevel< M >::value >= 0), int > = 0 >
      inline static decltype( auto ) namedSmootherTypes ( PriorityTag< 1 > )
      {
        return std::make_tuple( std::make_pair( std::string( "jacobi" ), Dune::MetaType< Dune::SeqJac< M, X, Y, UniformLeafLevel< M >::value > >() ),
                                std::make_pair( std::string( "gauss-seidel" ), Dune::MetaType< Dune::SeqGS< M, X, Y, UniformLeafLevel< M >::value > >() ),
                                std::make_pair( std::string( "sor" ), Dune::MetaType< Dune::SeqSOR< M, X, Y, UniformLeafLevel< M >::value > >() ),
                                std::make_pair( std::string( "ssor" ), Dune::MetaType< Dune::SeqSSOR< M, X, Y, UniformLeafLevel< M >::value > >() ) );
      }

      template< class M, class X, class Y >
      inline static decltype( auto ) namedSmootherTypes ( PriorityTag< 0 > )
      {
        return std::make_tuple();
      }

      template< class M, class X, class Y >
      inline static decltype( auto ) namedSmootherTypes ()
      {
        return namedSmootherTypes< M, X, Y >( PriorityTag< 42 >() );
      }



      // makeAMGPreconditioner
      // ---------------------

      template< class AssembledOperator, class Communication, std::enable_if_t< SupportsAMG< Communication >::value, int > = 0 >
      inline std::shared_ptr< Dune::Preconditioner< typename AssembledOperator::domain_type, typename AssembledOperator::range_type > >
      makeAMGPreconditioner ( const std::shared_ptr< AssembledOperator > op, const Communication &comm, Symmetry symmetry )
      {
        typedef typename AssembledOperator::matrix_type matrix_type;
        typedef typename AssembledOperator::domain_type domain_type;
        typedef typename AssembledOperator::range_type range_type;
        typedef typename Dune::FieldTraits< typename AssembledOperator::field_type >::real_type real_type;

        std::shared_ptr< Dune::Preconditioner< domain_type, range_type > > preconditioner;
        const auto smootherTypes = namedSmootherTypes< matrix_type, domain_type, range_type >();
        getEnum( "istl.preconditioner.smoother", smootherTypes, "jacobi", [ op, &comm, symmetry, &preconditioner ] ( auto type ) {
            typedef Dune::BlockPreconditioner< domain_type, range_type, Communication, typename decltype( type )::type > SmootherType;

            Dune::Amg::DefaultSmootherArgs< real_type > smootherArgs;
            smootherArgs.relaxationFactor = Dune::Fem::Parameter::getValue( "istl.preconditioner.relax", real_type( 1 ) );

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
        return preconditioner;
      }

      template< class AssembledOperator, class Communication, std::enable_if_t< !SupportsAMG< Communication >::value, int > = 0 >
      inline std::shared_ptr< Dune::Preconditioner< typename AssembledOperator::domain_type, typename AssembledOperator::range_type > >
      makeAMGPreconditioner ( const std::shared_ptr< AssembledOperator > op, const Communication &comm, Symmetry symmetry )
      {
        DUNE_THROW( Dune::InvalidStateException, "Communication does not support AMG." );
      }



      // makeSequentialPreconditioner
      // ----------------------------

      template< class Op >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::Richardson< typename Op::domain_type, typename Op::range_type > > &preconditioner )
      {
        typedef Dune::Richardson< typename Op::domain_type, typename Op::range_type > Preconditioner;
        preconditioner.reset( new Preconditioner( relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op, int l >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqJac< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > > &preconditioner )
      {
        typedef Dune::SeqJac< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), 1, relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op, int l >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqGS< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > > &preconditioner )
      {
        typedef Dune::SeqGS< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), 1, relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op, int l >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqSOR< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > > &preconditioner )
      {
        typedef Dune::SeqSOR< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), 1, relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op, int l >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqSSOR< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > > &preconditioner )
      {
        typedef Dune::SeqSSOR< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), 1, relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqILU0< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type > > &preconditioner )
      {
        typedef Dune::SeqILU0< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }

      template< class Op >
      inline void
      makeSequentialPreconditioner ( std::shared_ptr< Op > op, typename Dune::FieldTraits< typename Op::field_type >::real_type relaxationFactor,
                                     std::shared_ptr< Dune::SeqILDL< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type > > &preconditioner )
      {
        typedef Dune::SeqILDL< typename Op::matrix_type, typename Op::domain_type, typename Op::range_type > Preconditioner;
        preconditioner.reset( new Preconditioner( op->getmat(), relaxationFactor ), [ op ] ( Preconditioner *p ) { delete p; } );
      }



      // makeBlockPreconditioner
      // -----------------------

      template< class SeqPreconditioner, class Communication >
      inline std::shared_ptr< Dune::Preconditioner< typename SeqPreconditioner::domain_type, typename SeqPreconditioner::range_type > >
      makeBlockPreconditioner( std::shared_ptr< SeqPreconditioner > seqPreconditioner, const Communication &comm )
      {
        typedef typename SeqPreconditioner::domain_type domain_type;
        typedef typename SeqPreconditioner::range_type range_type;
        typedef Dune::Preconditioner< domain_type, range_type > Preconditioner;
        typedef Dune::BlockPreconditioner< domain_type, range_type, Communication, SeqPreconditioner > BlockPreconditioner;
        return std::shared_ptr< Preconditioner >( new BlockPreconditioner( *seqPreconditioner, comm ), [ seqPreconditioner ] ( Preconditioner *p ) { delete p; } );
      }



      // makePreconditioner
      // ------------------

      template< class Op, class Communication >
      inline std::shared_ptr< Dune::Preconditioner< typename Op::domain_type, typename Op::range_type > >
      makePreconditioner ( std::shared_ptr< Op > op, const Communication &comm, Symmetry symmetry )
      {
        typedef typename Op::matrix_type matrix_type;
        typedef typename Op::domain_type domain_type;
        typedef typename Op::range_type range_type;
        typedef typename Dune::FieldTraits< typename Op::field_type >::real_type real_type;

        const real_type relaxationFactor = Dune::Fem::Parameter::getValue( "istl.preconditioner.relax", real_type( 1 ) );

        const auto smootherTypes = namedSmootherTypes< matrix_type, domain_type, range_type >();

        std::shared_ptr< Dune::Preconditioner< domain_type, range_type > > preconditioner;
        const std::string preconditionerTypes[] = { "richardson", "smoother", "amg" };
        switch( Dune::Fem::Parameter::getEnum( "istl.preconditioner.type", preconditionerTypes, 1 ) )
        {
        case 0:
          {
            std::shared_ptr< Dune::Richardson< domain_type, range_type > > seqPreconditioner;
            makeSequentialPreconditioner( op, relaxationFactor, seqPreconditioner );
            preconditioner = makeBlockPreconditioner( seqPreconditioner, comm );
          }
          break;

        case 1:
          getEnum( "istl.preconditioner.smoother", smootherTypes, "jacobi", [ op, &comm, relaxationFactor, &preconditioner ] ( auto type ) {
              std::shared_ptr< typename decltype( type )::type > seqPreconditioner;
              makeSequentialPreconditioner( op, relaxationFactor, seqPreconditioner );
              preconditioner = makeBlockPreconditioner( seqPreconditioner, comm );
            } );
          break;

        case 2:
          preconditioner = makeAMGPreconditioner( op, comm, symmetry );
          break;

        default:
          DUNE_THROW( Dune::InvalidStateException, "Invalid ISTL preconditioner type selected." );
        }
        return preconditioner;
      }



      // InverseOperator
      // ---------------

      template< class LinearOperator, Symmetry symmetry = symmetric, template< class > class Communication = OwnerOverlapCopyCommunication >
      class InverseOperator final
        : public Dune::Fem::Operator< typename LinearOperator::RangeFunctionType, typename LinearOperator::DomainFunctionType >
      {
        static_assert( std::is_same< typename LinearOperator::DomainFunctionType, typename LinearOperator::RangeFunctionType >::value, "Domain function and range function must have the same type." );

      public:
        typedef LinearOperator LinearOperatorType;

        typedef typename LinearOperatorType::DomainFunctionType DiscreteFunctionType;

        typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef typename Dune::FieldTraits< typename DiscreteFunctionType::RangeFieldType >::real_type RealType;

      private:
        typedef VectorType< DiscreteFunctionType > vector_type;
        typedef MatrixType< LinearOperatorType > matrix_type;

      public:
        typedef Communication< DiscreteFunctionSpaceType > CommunicationType;

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

        InverseOperator ( RealType redEps, RealType absLimit, int maxIterations, bool verbose, const Dune::Fem::ParameterReader &parameter )
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



      // OperatorBasedPreconditionerFactory
      // ----------------------------------

      template< class Operator, Symmetry symmetry = symmetric, template< class > class Communication = OwnerOverlapCopyCommunication >
      class OperatorPreconditionerFactory
      {
        static_assert( std::is_same< typename Operator::DomainFunctionType, typename Operator::RangeFunctionType >::value, "Domain function and range function must have the same type." );

        typedef typename Operator::JacobianOperatorType LinearOperatorType;

        typedef typename LinearOperatorType::DomainFunctionType DiscreteFunctionType;

        typedef VectorType< DiscreteFunctionType > vector_type;
        typedef MatrixType< LinearOperatorType > matrix_type;

      public:
        typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef Communication< DiscreteFunctionSpaceType > CommunicationType;

        typedef Dune::Preconditioner< vector_type, vector_type > PreconditionerType;

        template< class... Args >
        OperatorPreconditionerFactory ( const DiscreteFunctionSpaceType &space, Args &&... args )
          : op_( std::forward< Args >( args )... ),
            zero_( "zero", space ),
            jacobian_( "preconditioner", space, space )
        {}

        template< class AssembledOperator >
        std::shared_ptr< PreconditionerType > operator() ( std::shared_ptr< AssembledOperator >, const CommunicationType &comm, Symmetry ) const
        {
          typedef Dune::OverlappingSchwarzOperator< matrix_type, vector_type, vector_type, CommunicationType > AssembledLinearOperatorType;
          op_.jacobian( zero_, jacobian_ );
          return makePreconditioner( std::make_shared< AssembledLinearOperatorType >( jacobian_.matrix(), comm ), comm, symmetry );
        }

        const Operator &op () const { return op_; }
        Operator &op () { return op_; }

      private:
        Operator op_;
        DiscreteFunctionType zero_;
        mutable LinearOperatorType jacobian_;
      };

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_ISTL_HH
