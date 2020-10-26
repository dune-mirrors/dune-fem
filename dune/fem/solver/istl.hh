#ifndef DUNE_FEM_SOLVER_ISTL_HH
#define DUNE_FEM_SOLVER_ISTL_HH

#if HAVE_DUNE_ISTL

#include <cassert>

// stl includes
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typelist.hh>
#include <dune/common/typeutilities.hh>

// dune-geometry includes
#include <dune/geometry/dimension.hh>

// dune-istl includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>

// dune-fem includes
#include <dune/fem/common/utility.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/space/mapper/parallel.hh>

// local includes
#include <dune/fem/solver/communication/fem.hh>
#include <dune/fem/solver/communication/hierarchical.hh>
#include <dune/fem/solver/communication/owneroverlapcopy.hh>



namespace Dune
{

  namespace Amg
  {

    // ConstructionTraits for SeqILDL
    // ------------------------------

    template<class M, class X, class Y>
    struct ConstructionTraits<SeqILDL<M, X, Y>>
    {
      using Arguments = DefaultConstructionArgs<SeqILDL<M, X, Y>>;

      static inline auto construct (Arguments& args) -> std::shared_ptr<SeqILDL<M, X, Y>>
      {
        return std::make_shared<SeqILDL<M, X, Y>>(args.getMatrix(), args.getArgs().relaxationFactor);
      }
    };

  } // namespace Amg



  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template<class DiscreteFunctionSpace>
    class HierarchicalDiscreteFunction;

    template<class DomainFunction, class RangeFunction>
    struct ISTLLinearOperator;



    namespace ISTL
    {

      // VectorType
      // ----------

      template< class DiscreteFunction >
      using VectorType = std::decay_t<decltype(std::declval<const DiscreteFunction&>().dofVector().array())>;



      // MatrixType
      // ----------

      template<class LinearOperator>
      struct __MatrixType
      {
        using Type = std::decay_t<decltype(std::declval<const LinearOperator&>().exportMatrix())>;
      };

      template<class DomainFunction, class RangeFunction>
      struct __MatrixType<Dune::Fem::ISTLLinearOperator<DomainFunction, RangeFunction>>
      {
        using Type = Dune::BCRSMatrix<typename Dune::Fem::ISTLLinearOperator<DomainFunction, RangeFunction>::LittleBlockType>;
      };


      template<class LinearOperator>
      using MatrixType = typename __MatrixType<LinearOperator>::Type;


      // __FillPrecondType
      // -----------------

      template<template<class, class, class, int...> class Prec, class Op, int ... l>
      using __FillPrecondType = Prec<typename Op::matrix_type, typename Op::domain_type, typename Op::range_type, l ...>;


      // IsBCRSMatrix
      // ------------

      template<class M>
      struct IsBCRSMatrix : std::false_type {};

      template<class B, class A>
      struct IsBCRSMatrix<Dune::BCRSMatrix<B, A>> : std::true_type {};



      // Symmetry
      // --------

      enum Symmetry : bool { unsymmetric = false, symmetric = true };



      // UniformLeafLevel
      // ----------------

      template<class T, class = void>
      struct UniformLeafLevelType;

      template<class T>
      using UniformLeafLevel = typename UniformLeafLevelType<T>::Type;

      template<class V1, class... V>
      struct UniformLeafLevelType<Dune::MultiTypeBlockVector<V1, V...>, std::enable_if_t<Std::are_all_same<UniformLeafLevel<V1>, UniformLeafLevel<V>...>::value>>
      {
        using Type = std::integral_constant<int, UniformLeafLevel<V1>::value + 1>;
      };

      template<class B, class A>
      struct UniformLeafLevelType<Dune::BlockVector<B, A>>
      {
        using Type = std::integral_constant<int, UniformLeafLevel<B>::value + 1>;
      };

      template<class K, int n>
      struct UniformLeafLevelType<Dune::FieldVector<K, n>>
      {
        using Type = std::integral_constant<int, 0>;
      };

      template<class R1, class... R>
      struct UniformLeafLevelType<Dune::MultiTypeBlockMatrix< R1, R...>, std::enable_if_t<Std::are_all_same<UniformLeafLevel<R1>, UniformLeafLevel<R>...>::value>>
      {
        // Note: The rows of MultiTypeBlockMatrix are of type MultiTypeBlockVector, which already increases the level
        using Type = std::integral_constant<int, UniformLeafLevel<R1>::value>;
      };

      template<class B, class A>
      struct UniformLeafLevelType<Dune::BCRSMatrix<B, A>>
      {
        using Type = std::integral_constant<int, UniformLeafLevel<B>::value + 1>;
      };

      template<class K, int m, int n>
      struct UniformLeafLevelType<Dune::FieldMatrix<K, m, n>>
      {
        using Type = std::integral_constant<int, 0>;
      };



      // namedSmootherTypes
      // ------------------

      template<class M, class X, class Y,
               std::enable_if_t<IsBCRSMatrix<M>::value && (UniformLeafLevel<M>::value == 1), int> = 0>
      inline static decltype(auto) namedSmootherTypes (PriorityTag<2>)
      {
        return std::make_tuple(std::make_pair(std::string("jacobi"), Dune::MetaType<Dune::SeqJac< M, X, Y, 1>>{}),
                               std::make_pair(std::string("gauss-seidel"), Dune::MetaType<Dune::SeqGS< M, X, Y, 1>>{}),
                               std::make_pair(std::string("sor"), Dune::MetaType<Dune::SeqSOR<M, X, Y, 1>>{}),
                               std::make_pair(std::string("ssor"), Dune::MetaType<Dune::SeqSSOR<M, X, Y, 1>>{}),
                               std::make_pair(std::string("ilu"), Dune::MetaType<Dune::SeqILU< M, X, Y>>{}),
                               std::make_pair(std::string("ildl"), Dune::MetaType<Dune::SeqILDL<M, X, Y>>{}));
      }

      template<class M, class X, class Y,
               std::enable_if_t<(UniformLeafLevel<M>::value >= 0), int> = 0>
      inline static decltype(auto) namedSmootherTypes (PriorityTag<1>)
      {
        return std::make_tuple(std::make_pair(std::string("jacobi"), Dune::MetaType<Dune::SeqJac<M, X, Y, UniformLeafLevel<M>::value>>{}),
                               std::make_pair(std::string("gauss-seidel"), Dune::MetaType<Dune::SeqGS<M, X, Y, UniformLeafLevel<M>::value>>{}),
                               std::make_pair(std::string("sor"), Dune::MetaType<Dune::SeqSOR<M, X, Y, UniformLeafLevel<M>::value>>{}),
                               std::make_pair(std::string("ssor"), Dune::MetaType<Dune::SeqSSOR<M, X, Y, UniformLeafLevel<M>::value>>{}));
      }

      template<class M, class X, class Y>
      inline static decltype(auto) namedSmootherTypes (PriorityTag<0>)
      {
        return std::make_tuple();
      }

      template<class M, class X, class Y>
      inline static decltype(auto) namedSmootherTypes ()
      {
        return namedSmootherTypes<M, X, Y>(PriorityTag<42>{});
      }


      inline static decltype(auto) namedAMGNormTypes ()
      {
        return std::make_tuple(std::make_pair(std::string("firstdiagonal"), Dune::MetaType<Dune::Amg::FirstDiagonal>{}),
                               std::make_pair(std::string("rowsum"), Dune::MetaType<Dune::Amg::RowSum>{}),
                               std::make_pair(std::string("frobenius"), Dune::MetaType<Dune::Amg::FrobeniusNorm>{}),
                               std::make_pair(std::string("one"), Dune::MetaType<Dune::Amg::AlwaysOneNorm>{}));
      }



      // SolverParameter
      // ---------------

      class SolverParameter : public LocalParameter<Fem::SolverParameter, SolverParameter>
      {
        using BaseType = LocalParameter<Fem::SolverParameter, SolverParameter>;

      public:
        using BaseType::parameter;
        using BaseType::keyPrefix;

      private:
        template<class... T, class F>
        void getEnum (const std::string& key,
                      const std::tuple<std::pair<std::string, T>...>& choices,
                      const std::string& defaultValue,
                      F&& f) const
        {
          const std::string& value = parameter().getValue(key, defaultValue);
          bool success = false;
          Hybrid::forEach(choices, [&f, &value, &success] (const auto& choice) {
            if (value != choice.first)
              return;

            assert(!success);
            success = true;

            f(choice.second);
          });
          if (success)
            return;

          std::string values;
          Hybrid::forEach(choices, [&values] (const auto& choice) { values += (values.empty() ? "'" : ", '") + choice.first + "'"; });
          DUNE_THROW(ParameterInvalid, "Parameter '" << key << "' invalid (choices are: " << values << ").");
        }

        // hide functions
        using BaseType::gmresRestart;
        using BaseType::relaxation;
        using BaseType::preconditionerIteration;
        using BaseType::preconditionerLevel;
        using BaseType::preconditionMethod;
        using BaseType::preconditionMethodTable;
        using BaseType::verbose;
        using BaseType::setVerbose;

      public:
        SolverParameter (const ParameterReader& parameter = Parameter::container())
          : BaseType("istl.", parameter)
        {}

        SolverParameter (const std::string& keyPrefix, const ParameterReader& parameter = Parameter::container())
         : BaseType(keyPrefix, parameter)
        {}

        SolverParameter (const Fem::SolverParameter& parameter)
          : SolverParameter(parameter.keyPrefix(), parameter.parameter())
        {}

        virtual int errorMeasure () const override
        {
          if (errorMeasure_ < 0)
            errorMeasure_ = BaseType::errorMeasure();
          return errorMeasure_;
        }

        virtual int verbosity () const
        {
          const std::string verbosityTable[] = { "off", "on", "full" };
          return parameter().getEnum(keyPrefix() + "verbosity", verbosityTable, 0);
        }

        virtual int restart () const { return parameter().getValue(keyPrefix() + "restart", 20); }
        virtual int fcgMmax () const { return parameter().getValue(keyPrefix() + "fcg.mmax", 10); }

        virtual double precRelaxation () const { return parameter().getValue(keyPrefix() + "preconditioner.relax", 1.0); }
        virtual int precIterations () const { return parameter().getValue(keyPrefix() + "preconditioner.iterations", 1); }

        virtual int precType () const
        {
          const std::string types[] = { "richardson", "smoother", "amg" };
          return parameter().getEnum(keyPrefix() + "preconditioner.type", types, 1);
        }

        template<class... T, class F>
        void precSmoother (const std::tuple<std::pair<std::string, T>...>& choices, const std::string& defaultChoice, F&& f) const
        {
          getEnum(keyPrefix() + "preconditioner.smoother", choices, defaultChoice, std::forward<F>(f));
        }

        virtual int iluFillin () const { return parameter().getValue(keyPrefix() + "ilu.fillin", 0); }
        virtual bool iluReorder () const { return parameter().getValue(keyPrefix() + "ilu.reorder", false); }

        virtual int amgMaxLevel () const { return parameter().getValue(keyPrefix() + "amg.maxlevel", 100); }
        virtual int amgCoarsenTarget () const { return parameter().getValue(keyPrefix() + "amg.coarsentarget", 1000); }
        virtual double amgMinCoarsenRate () const { return parameter().getValue(keyPrefix() + "amg.mincoarsenrate", 1.2); }
        virtual double amgProlongDamping () const { return parameter().getValue(keyPrefix() + "amg.prolongation.dampingfactor", 1.6); }
        virtual int amgDebugLevel () const { return parameter().getValue(keyPrefix() + "amg.debuglevel", 0); }
        virtual int amgPreSmoothSteps () const { return parameter().getValue(keyPrefix() + "amg.presmoothsteps", 2); }
        virtual int amgPostSmoothSteps () const { return parameter().getValue(keyPrefix() + "amg.postsmoothsteps", 2); }
        virtual bool amgAdditive () const { return parameter().getValue(keyPrefix() + "amg.additive", false); }
        virtual double amgAlpha () const { return parameter().getValue(keyPrefix() + "amg.alpha", 1.0/3.0); }
        virtual double amgBeta () const { return parameter().getValue(keyPrefix() + "amg.beta", 1.0e-5); }

        virtual int amgCycle () const
        {
          const std::string cycles[] = { "v-cycle", "w-cycle" };
          return 1 + parameter().getEnum(keyPrefix() + "amg.cycle", cycles, 0);
        }

        virtual int amgAccumulate () const
        {
          const std::string modes[] = { "none", "once", "successive" };
          return parameter().getEnum(keyPrefix() + "amg.accumulate", modes, 2);
        }

        virtual std::size_t amgAggregationDimension () const { return parameter().getValue<std::size_t>(keyPrefix() + "amg.aggregation.dimension", 2); }
        virtual std::size_t amgAggregationDistance () const { return parameter().getValue<std::size_t>(keyPrefix() + "amg.aggregation.distance", 2); }
        virtual std::size_t amgAggregationMinSize () const { return parameter().getValue<std::size_t>(keyPrefix() + "amg.aggregation.min", 4); }
        virtual std::size_t amgAggregationMaxSize () const { return parameter().getValue<std::size_t>(keyPrefix() + "amg.aggregation.max", 6); }
        virtual std::size_t amgAggregationConnectivity () const { return parameter().getValue<std::size_t>(keyPrefix() + "amg.aggregation.connectivity", 15); }
        virtual bool amgAggregationSkipIsolated () const { return parameter().getValue(keyPrefix() + "amg.aggregation.skipisolated", false); }

        virtual int amgAggregationMode () const
        {
          const std::string modes[] = { "isotropic", "anisotropic", "parameter" };
          return parameter().getEnum(keyPrefix() + "amg.aggregation.mode", modes, 0);
        }

        template<class... T, class F>
        void amgNorm (const std::tuple<std::pair<std::string, T>...>& choices, const std::string& defaultChoice, F&& f) const
        {
          getEnum(keyPrefix() + "amg.norm", choices, defaultChoice, std::forward<F>(f));
        }

      private:
        mutable int errorMeasure_ = -1;
      };



      // makeAMGPreconditioner
      // ---------------------

      template<class AssembledOperator, class Communication,
               std::enable_if_t<SupportsAMG<Communication>::value, int> = 0>
      inline auto makeAMGPreconditioner (const std::shared_ptr<AssembledOperator> op, const Communication& comm, Symmetry symmetry,
                                         const SolverParameter& parameter = {})
        -> std::shared_ptr<Dune::Preconditioner<typename AssembledOperator::domain_type, typename AssembledOperator::range_type>>
      {
        using matrix_type = typename AssembledOperator::matrix_type;
        using domain_type = typename AssembledOperator::domain_type;
        using range_type  = typename AssembledOperator::range_type;
        using real_type   = real_t<typename AssembledOperator::field_type>;

        std::shared_ptr<Dune::Preconditioner<domain_type, range_type>> preconditioner;
        const auto smoothers = namedSmootherTypes<matrix_type, domain_type, range_type>();
        parameter.precSmoother(smoothers, "jacobi", [op, &comm, symmetry, &parameter, &preconditioner] (auto type) {
            using SmootherType  = Dune::BlockPreconditioner<domain_type, range_type, Communication, typename decltype(type)::type>;
            using AMG           = Dune::Amg::AMG<AssembledOperator, domain_type, SmootherType, Communication>;

            Dune::Amg::DefaultSmootherArgs <real_type> smootherArgs;
            smootherArgs.relaxationFactor = parameter.precRelaxation();
            smootherArgs.iterations = parameter.precIterations();

            Dune::Amg::Parameters amgParams(
              parameter.amgMaxLevel(),
              parameter.amgCoarsenTarget(),
              parameter.amgMinCoarsenRate(),
              parameter.amgProlongDamping(),
              static_cast<Dune::Amg::AccumulationMode>(parameter.amgAccumulate())
            );

            // parameters
            switch(parameter.amgAggregationMode())
            {
              case 0:
                amgParams.setDefaultValuesIsotropic(parameter.amgAggregationDimension(), parameter.amgAggregationDistance()); // ...
                break;

              case 1:
                amgParams.setDefaultValuesAnisotropic(parameter.amgAggregationDimension(), parameter.amgAggregationDistance()); // ...
                break;

              case 2:
                amgParams.setMaxDistance(parameter.amgAggregationDistance());
                amgParams.setMinAggregateSize(parameter.amgAggregationMinSize());
                amgParams.setMaxAggregateSize(parameter.amgAggregationMaxSize());
                break;
            }

            amgParams.setMaxConnectivity(parameter.amgAggregationConnectivity());
            amgParams.setSkipIsolated(parameter.amgAggregationSkipIsolated());

            amgParams.setDebugLevel(parameter.amgDebugLevel());
            amgParams.setNoPreSmoothSteps(parameter.amgPreSmoothSteps());
            amgParams.setNoPostSmoothSteps(parameter.amgPostSmoothSteps());
            amgParams.setAlpha(parameter.amgAlpha());
            amgParams.setBeta(parameter.amgBeta());
            amgParams.setGamma(parameter.amgCycle());
            amgParams.setAdditive(parameter.amgAdditive());

            parameter.amgNorm(namedAMGNormTypes(), "rowsum", [op, &comm, symmetry, &amgParams, &smootherArgs, &preconditioner](auto type){
                if (symmetry == symmetric)
                {
                  Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion< matrix_type, typename decltype(type)::type>> criterion(amgParams);
                  preconditioner.reset(new AMG(*op, criterion, smootherArgs, comm), [op] (AMG* p) { delete p; });
                }
                else
                {
                  Dune::Amg::CoarsenCriterion<Dune::Amg::UnSymmetricCriterion<matrix_type, typename decltype(type)::type>> criterion( amgParams );
                  preconditioner.reset(new AMG(*op, criterion, smootherArgs, comm), [op] (AMG* p) { delete p; });
                }
            });
          });
        return preconditioner;
      }

      template<class AssembledOperator, class Communication,
               std::enable_if_t<!SupportsAMG<Communication>::value, int> = 0>
      inline auto makeAMGPreconditioner (const std::shared_ptr<AssembledOperator> op, const Communication& comm, Symmetry symmetry,
                                         const SolverParameter& parameter = {})
        -> std::shared_ptr<Dune::Preconditioner<typename AssembledOperator::domain_type, typename AssembledOperator::range_type>>
      {
        DUNE_THROW(Dune::InvalidStateException, "Communication does not support AMG.");
      }



      // makeSequentialPreconditioner
      // ----------------------------

      template<class Op>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<Dune::Richardson<typename Op::domain_type, typename Op::range_type>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = Dune::Richardson<typename Op::domain_type, typename Op::range_type>;
        preconditioner.reset(new Preconditioner(parameter.precRelaxation()),
                             [op] ( Preconditioner* p) { delete p; });
      }

      template<class Op, int l>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<__FillPrecondType<Dune::SeqJac, Op, l>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = __FillPrecondType<Dune::SeqJac, Op, l>;
        preconditioner.reset(new Preconditioner(op->getmat(), parameter.precIterations(), parameter.precRelaxation()),
                             [op] (Preconditioner* p) { delete p; });
      }

      template<class Op, int l>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<__FillPrecondType<Dune::SeqSOR, Op, l>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = __FillPrecondType<Dune::SeqSOR, Op, l>;
        preconditioner.reset(new Preconditioner(op->getmat(), parameter.precIterations(), parameter.precRelaxation()),
                             [op] (Preconditioner* p) { delete p; });
      }

      template<class Op, int l>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<__FillPrecondType<Dune::SeqSSOR, Op, l>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = __FillPrecondType<Dune::SeqSSOR, Op, l>;
        preconditioner.reset(new Preconditioner(op->getmat(), parameter.precIterations(), parameter.precRelaxation()),
                             [op] (Preconditioner* p) { delete p; });
      }

      template<class Op>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<__FillPrecondType<Dune::SeqILU, Op>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = __FillPrecondType<Dune::SeqILU, Op>;
        preconditioner.reset(new Preconditioner(op->getmat(), parameter.iluFillin(), parameter.precRelaxation(), parameter.iluReorder()),
                             [op] (Preconditioner* p) { delete p; });
      }

      template<class Op>
      inline void makeSequentialPreconditioner (std::shared_ptr<Op> op,
                                                std::shared_ptr<__FillPrecondType<Dune::SeqILDL, Op>>& preconditioner,
                                                const SolverParameter& parameter = {})
      {
        using Preconditioner = __FillPrecondType<Dune::SeqILDL, Op>;
        preconditioner.reset(new Preconditioner(op->getmat(), parameter.precRelaxation()),
                             [op] (Preconditioner* p) { delete p; });
      }



      // makeBlockPreconditioner
      // -----------------------

      template<class SeqPreconditioner, class Communication>
      inline auto makeBlockPreconditioner(std::shared_ptr<SeqPreconditioner> seqPreconditioner, const Communication& comm)
        -> std::shared_ptr<Dune::Preconditioner<typename SeqPreconditioner::domain_type, typename SeqPreconditioner::range_type>>
      {
        using domain_type = typename SeqPreconditioner::domain_type;
        using range_type  = typename SeqPreconditioner::range_type;

        using BlockPreconditioner = Dune::BlockPreconditioner<domain_type, range_type, Communication, SeqPreconditioner>;

        return std::make_shared<BlockPreconditioner>(seqPreconditioner, comm);
      }



      // makePreconditioner
      // ------------------

      template<class Op, class Communication>
      inline auto makePreconditioner (std::shared_ptr<Op> op, const Communication& comm, Symmetry symmetry,
                                      const SolverParameter& parameter = {})
        -> std::shared_ptr<Dune::Preconditioner<typename Op::domain_type, typename Op::range_type>>
      {
        using domain_type = typename Op::domain_type;
        using range_type  = typename Op::range_type;

        const auto smoothers = namedSmootherTypes<typename Op::matrix_type, domain_type, range_type>();
        std::shared_ptr<Dune::Preconditioner<domain_type, range_type>> preconditioner;

        switch (parameter.precType())
        {
        case 0:
          {
            std::shared_ptr<Dune::Richardson<domain_type, range_type>> seqPreconditioner;
            makeSequentialPreconditioner(op, seqPreconditioner, parameter);
            preconditioner = makeBlockPreconditioner(seqPreconditioner, comm);
          }
          break;

        case 1:
          parameter.precSmoother(smoothers, "jacobi", [op, &comm, &parameter, &preconditioner] (auto type) {
              std::shared_ptr<typename decltype(type)::type> seqPreconditioner;
              makeSequentialPreconditioner(op, seqPreconditioner, parameter);
              preconditioner = makeBlockPreconditioner(seqPreconditioner, comm);
            } );
          break;

        case 2:
          preconditioner = makeAMGPreconditioner(op, comm, symmetry, parameter);
          break;

        default:
          DUNE_THROW(Dune::InvalidStateException, "Invalid ISTL preconditioner type selected.");
        }
        return preconditioner;
      }



      // InverseOperator
      // ---------------

      template<class LinearOperator, Symmetry symmetry = symmetric, template<class> class Communication = OwnerOverlapCopyCommunication>
      class InverseOperator final
        : public Dune::Fem::Operator<typename LinearOperator::RangeFunctionType, typename LinearOperator::DomainFunctionType>
      {
        static_assert(std::is_same< typename LinearOperator::DomainFunctionType, typename LinearOperator::RangeFunctionType>::value,
                      "Domain function and range function must have the same type.");

      public:
        using LinearOperatorType = LinearOperator;
        using OperatorType = LinearOperatorType;

        using DiscreteFunctionType      = typename LinearOperatorType::DomainFunctionType;
        using DiscreteFunctionSpaceType = typename DiscreteFunctionType::DiscreteFunctionSpaceType;

        using RealType = real_t<typename DiscreteFunctionType::RangeFieldType>;

        using SolverParameterType = SolverParameter;

      private:
        using vector_type = VectorType<DiscreteFunctionType>;
        using matrix_type = MatrixType<LinearOperatorType>;

      public:
        using CommunicationType = Communication<DiscreteFunctionSpaceType>;

        using AssembledLinearOperatorType = Dune::OverlappingSchwarzOperator<matrix_type, vector_type, vector_type, CommunicationType>;
        using PreconditionerType          = Dune::Preconditioner<vector_type, vector_type>;
        using ScalarProductType           = Dune::ScalarProduct<vector_type>;

        using PreconditionerFactory = std::function<std::shared_ptr<PreconditionerType>(std::shared_ptr<AssembledLinearOperatorType>, const CommunicationType&, Symmetry, const SolverParameterType&)>;

        InverseOperator (const SolverParameterType& parameter = {})
          : InverseOperator(makePreconditioner<AssembledLinearOperatorType, CommunicationType>, parameter)
        {}

        InverseOperator (PreconditionerFactory factory, const SolverParameterType& parameter = {})
          : preconditionerFactory_{std::move(factory)},
            parameter_{std::make_shared<SolverParameterType>(parameter)}
        {}

        InverseOperator (const LinearOperatorType& op, const SolverParameterType& parameter = {})
          : InverseOperator(parameter)
        {
          bind(op);
        }

        InverseOperator (const LinearOperatorType& op, PreconditionerFactory factory, const SolverParameterType& parameter = {})
          : InverseOperator(std::move(factory), parameter)
        {
          bind(op);
        }

        void operator() (const DiscreteFunctionType& u, DiscreteFunctionType& w) const override
        {
          result_.clear();
          vector_type b(u.dofVector().array());
          if (parameter_->errorMeasure() == 0)
          {
            vector_type residuum(b);
            linearOperator_->applyscaleadd(-1.0, w.dofVector().array(), residuum);
            RealType res = scalarProduct_->norm(residuum);
            RealType reduction = (res > 0) ? parameter_->tolerance() / res : 1e-3;
            solver_->apply(w.dofVector().array(), b, reduction, result_);
          }
          else
            solver_->apply(w.dofVector().array(), b, result_);
        }

        void bind (const LinearOperatorType& op)
        {
          buildCommunication(op.domainSpace(), Dune::SolverCategory::overlapping, communication_);

          linearOperator_.reset(new AssembledLinearOperatorType(op.exportMatrix(), *communication_));
          scalarProduct_.reset(new Dune::OverlappingSchwarzScalarProduct<vector_type, CommunicationType>(*communication_));

          // create preconditioner
          preconditioner_ = preconditionerFactory_(linearOperator_, *communication_, symmetry, parameter());

          // create linear solver
          auto reduction      = parameter().tolerance();
          auto maxIterations  = parameter().maxIterations();
          auto verbosity      = op.domainSpace().gridPart().comm().rank() == 0 ? parameter().verbosity() : 0;

          switch (parameter().solverMethod({SolverParameterType::cg,
                                            SolverParameterType::bicgstab,
                                            SolverParameterType::gmres,
                                            SolverParameterType::minres,
                                            SolverParameterType::gradient,
                                            SolverParameterType::loop},
                                           {"pcg", "fcg", "fgmres"},
                                           (symmetry == symmetric ? 0 : 1)))
          {
          case -2:
            solver_.reset(new Dune::RestartedFlexibleGMResSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, parameter().restart(), maxIterations, verbosity));
            break;

          case -1:
            solver_.reset(new Dune::RestartedFCGSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity, parameter().fcgMmax()));
            break;

          case 0:
            solver_.reset(new Dune::GeneralizedPCGSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity, parameter().restart()));
            break;

          case 1:
            solver_.reset(new Dune::CGSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity));
            break;

          case 2:
            solver_.reset(new Dune::BiCGSTABSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity));
            break;

          case 3:
            solver_.reset(new Dune::RestartedGMResSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, parameter().restart(), maxIterations, verbosity));
            break;

          case 4:
            solver_.reset(new Dune::MINRESSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity));
            break;

          case 5:
            solver_.reset(new Dune::GradientSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity));
            break;

          case 6:
            solver_.reset(new Dune::LoopSolver<vector_type>(linearOperator_, scalarProduct_, preconditioner_, reduction, maxIterations, verbosity));
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

        auto parameter () -> SolverParameterType& { return *parameter_; }
        void setParamters (const SolverParameterType& parameter) { parameter_ = std::make_shared<SolverParameterType>(parameter); }

        int iterations () const { return result_.iterations; }
        bool converged () const { return result_.converged; }

        void setMaxIterations (int maxIterations) { parameter().setMaxIterations(maxIterations); }

      private:
        PreconditionerFactory preconditionerFactory_;
        std::shared_ptr<SolverParameterType> parameter_;

        std::shared_ptr<CommunicationType> communication_;
        std::shared_ptr<AssembledLinearOperatorType> linearOperator_;
        std::shared_ptr<ScalarProductType> scalarProduct_;
        std::shared_ptr<PreconditionerType> preconditioner_;
        std::shared_ptr<Dune::InverseOperator<vector_type, vector_type>> solver_;

        mutable Dune::InverseOperatorResult result_;
      };



      // OperatorBasedPreconditionerFactory
      // ----------------------------------

      template<class Operator, Symmetry symmetry = symmetric, template<class> class Communication = OwnerOverlapCopyCommunication>
      class OperatorPreconditionerFactory
      {
        static_assert(std::is_same<typename Operator::DomainFunctionType, typename Operator::RangeFunctionType>::value,
                      "Domain function and range function must have the same type.");

        using LinearOperatorType    = typename Operator::JacobianOperatorType;
        using DiscreteFunctionType  = typename LinearOperatorType::DomainFunctionType;

        using vector_type = VectorType<DiscreteFunctionType>;
        using matrix_type = MatrixType<LinearOperatorType>;

      public:
        using DiscreteFunctionSpaceType = typename DiscreteFunctionType::DiscreteFunctionSpaceType;

        using CommunicationType = Communication<DiscreteFunctionSpaceType>;

        using PreconditionerType = Dune::Preconditioner<vector_type, vector_type>;

        template<class... Args>
        OperatorPreconditionerFactory (const DiscreteFunctionSpaceType& space, Args&& ...args)
          : op_(std::forward<Args>(args)...),
            zero_("zero", space),
            jacobian_("preconditioner", space, space)
        {}

        template<class AssembledOperator>
        auto operator() (std::shared_ptr<AssembledOperator>, const CommunicationType& comm, Symmetry,
                         const SolverParameter& parameter = {}) const
          -> std::shared_ptr<PreconditionerType>
        {
          using AssembledLinearOperatorType = Dune::OverlappingSchwarzOperator<matrix_type, vector_type, vector_type, CommunicationType>;
          op_.jacobian(zero_, jacobian_);
          return makePreconditioner(std::make_shared<AssembledLinearOperatorType>(jacobian_.matrix(), comm), comm, symmetry, parameter);
        }

        auto op () const -> const Operator& { return op_; }
        auto op () -> Operator& { return op_; }

      private:
        Operator op_;
        DiscreteFunctionType zero_;
        mutable LinearOperatorType jacobian_;
      };

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifdef HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_SOLVER_ISTL_HH
