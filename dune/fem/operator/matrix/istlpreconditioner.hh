#ifndef DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH
#define DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH

#include <memory>
#include <type_traits>

// standard diagonal preconditioner
#include <dune/fem/solver/diagonalpreconditioner.hh>

#if HAVE_DUNE_ISTL
#include <dune/common/version.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#endif

#include <dune/fem/solver/parameter.hh> // SolverParameter
#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

namespace Dune
{

  namespace Fem
  {

    struct ISTLPreconditionMethods
    {
      static std::vector< int > supportedPreconditionMethods() {
        return std::vector< int > ({ SolverParameter::none,        // no preconditioner
                                     SolverParameter::ssor,        // SSOR preconditioner
                                     SolverParameter::sor ,        // SOR preconditioner
                                     SolverParameter::ilu ,        // ILU preconditioner
                                     SolverParameter::gauss_seidel,// Gauss-Seidel preconditioner
                                     SolverParameter::jacobi,      // Jacobi preconditioner
                                     SolverParameter::amg_ilu,     // AMG with ILU-0 smoother (deprecated)
                                     SolverParameter::amg_jacobi,  // AMG with Jacobi smoother
                                     SolverParameter::ildl         // ILDL from istl
                                    });
      }
    };


    template <class MatrixImp>
    class LagrangeParallelMatrixAdapter;

    template <class MatrixImp>
    class LagrangeParallelMatrixAdapter;

    template <class MatrixImp>
    class DGParallelMatrixAdapter;

    struct ISTLSolverParameter : public LocalParameter< SolverParameter, ISTLSolverParameter >
    {
      typedef LocalParameter< SolverParameter, ISTLSolverParameter > BaseType;
    public:
      using BaseType :: parameter ;
      using BaseType :: keyPrefix ;

      ISTLSolverParameter( const ParameterReader &parameter = Parameter::container() )
        : BaseType( parameter )
      {}

      ISTLSolverParameter( const std::string &keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : BaseType( keyPrefix, parameter )
      {}

      ISTLSolverParameter( const SolverParameter& sp )
        : ISTLSolverParameter( sp.keyPrefix(), sp.parameter() )
      {}

      virtual double overflowFraction () const
      {
        return parameter().getValue< double >( keyPrefix() + "matrix.overflowfraction", 1.0 );
      }

      virtual bool fastILUStorage () const
      {
        return parameter().getValue< bool >( keyPrefix() + "preconditioning.fastilustorage", true );
      }
    };

#if HAVE_DUNE_ISTL
    /** \brief Iterative linear methods as preconditioners. Implemented are
     *         Jacobi, Gauss-Seidel, SOR and SSOR.
     *
     *         Note: On rows corresponding to border DoFs this scheme falls back to Jacobi, to avoid
     *               complicated coloring strategies. This will lead to a
     *               slightly increased number of linear iterations in parallel
     *               runs.
     */
    template< class MatrixObject, class X, class Y, int method, int l=1 >
    class ParallelIterative : public Preconditioner<X,Y>
    {
    public:
      typedef typename MatrixObject :: ColumnDiscreteFunctionType DiscreteFunctionType ;

    protected:
      typedef typename X::field_type field_type;
      typedef typename MatrixObject::MatrixType  MatrixType;
      typedef MatrixType M;
      typedef typename M::block_type MatrixBlockType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType ;
      typedef typename DiscreteFunctionSpaceType:: DomainFieldType DomainFieldType;

      static const bool isNumber = (M::block_type::rows == M::block_type::cols) && (M::block_type::rows == 1);

      typedef FieldVector< field_type, MatrixBlockType::rows*MatrixBlockType::cols > FlatMatrixBlock;
      // BlockVector from dune-istl
      typedef BlockVector< FlatMatrixBlock > DiagonalType;

      typedef FieldMatrixConverter< FlatMatrixBlock, MatrixBlockType > MatrixConverterType;

      // Implementation of SOR and Jacobi's method
      // Note: for SOR we have xnew == xold
      template <class rowiterator, bool forward>
      void iterate (const M& A, X& xnew, const X& xold, const Y& b, const field_type& w,
                    const DiagonalType& diagInv,
                    rowiterator i,
                    const rowiterator endi,
                    const std::integral_constant<bool, forward> fb ) const
      {
        typedef typename M::ConstColIterator coliterator;

        typedef typename Y::block_type bblock;
        typedef typename X::block_type xblock;

        bblock rhs;
        xblock v;

        // Initialize nested data structure if there are entries
        if( i != endi )
          v=xnew[0];

        const bool hasDiagInv = diagInv.size() > 0 ;

        const auto nextRow = [&i]()
        {
          // increment/decrement iterator
          if constexpr ( forward )
            ++i;
          else
            --i;
        };

        for (; i!=endi; nextRow() )
        {
          const coliterator endj=(*i).end();           // iterate over a_ij with j < i
          coliterator j=(*i).begin();
          const auto rowIdx = i.index();

          // treat missing rows as unit rows
          if( j == endj )
          {
            xnew[rowIdx] += w*b[rowIdx];
            continue ;
          }

          rhs = b[rowIdx];           // rhs = b_i
          if constexpr (isNumber)
          {
            // note that for SOR xnew == xold
            for (; j.index()<rowIdx; ++j)
              rhs -= (*j) * xold[j.index()];  //  rhs -= sum_{j<i} a_ij * x_j

            // not needed, since we store the diagonal separately
            coliterator diag=j;               // *diag = a_ii
            for (; j!=endj; ++j)
              rhs -= (*j) * xold[j.index()];  //  rhs -= sum_{j<i} a_ij * x_j

            // v = rhs / diag
            if( hasDiagInv )
              v = rhs * diagInv[ rowIdx ];
            else
              v = rhs / (*diag);

            xnew[ rowIdx ] += w*v;            //  x_i = w / a_ii * (b_i - sum_{j<i} a_ij * xnew_j - sum_{j>=i} a_ij * xold_j)
          }
          else
          {
            // note that for SOR xnew == xold
            for (; j.index()<rowIdx; ++j)
              (*j).mmv(xold[j.index()],rhs);  //  rhs -= sum_{j<i} a_ij * x_j
            coliterator diag=j;               // *diag = a_ii
            for (; j!=endj; ++j)
              (*j).mmv(xold[j.index()],rhs);  //  rhs -= sum_{j<i} a_ij * x_j

            // v = rhs / diag
            if( hasDiagInv )
            {
              MatrixConverterType m( diagInv[ rowIdx ] );
              m.solve( v, rhs );
            }
            else
              (*diag).solve(v, rhs);
            xnew[rowIdx].axpy(w,v);           //  x_i = w / a_ii * (b_i - sum_{j<i} a_ij * xnew_j - sum_{j>=i} a_ij * xold_j)
          }
        }
      }

      // Implementation of SOR and Jacobi's method
      // Note: for SOR we have xnew == xold
      void forwardIterate (const M& A, X& xnew, const X& xold, const Y& b, const field_type& w,
                           const DiagonalType& diagInv ) const
      {
        bool runParallel = threading_ && (&xnew != &xold); // Jacobi only

        if( runParallel )
        {
          auto doIterate = [this, &A, &xnew, &xold, &b, &w, &diagInv] ()
          {
            const auto slice = A.sliceBeginEnd( MPIManager::thread(), MPIManager::numThreads() );
            // standard forwards iteration
            iterate( A, xnew, xold, b, w, diagInv,
                     A.slicedBegin( slice.first ), A.slicedEnd( slice.second ), std::true_type() );
          };

          try {
            // execute in parallel
            MPIManager :: run ( doIterate );
          }
          catch ( const SingleThreadModeError& e )
          {
            runParallel = false;
          }
        }

        if( ! runParallel )
        {
          // serial version
          iterate( A, xnew, xold, b, w, diagInv, A.begin(), A.end(), std::true_type() );
        }
      }

      // Implementation of SOR and Jacobi's method
      // Note: for SOR we have xnew == xold
      void backwardIterate (const M& A, X& xnew, const X& xold, const Y& b, const field_type& w,
                            const DiagonalType& diagInv ) const
      {
        // backwards iteration (use r-iterators)
        iterate( A, xnew, xold, b, w, diagInv, A.beforeEnd(), A.beforeBegin(), std::false_type() );
      }

    protected:
      const MatrixObject& mObj_;
      const MatrixType& _A_;
      const int _n;
      const double _w;

      DiagonalType diagonalInv_;

      const bool threading_;

    public:
      ParallelIterative( const MatrixObject& mObj, int n=1, double relax=1.0  )
        : mObj_( mObj ),
          _A_( mObj.exportMatrix() ),
          _n( n ), _w(relax),
          threading_( mObj.threading()  )
      {
        Dune::CheckIfDiagonalPresent<MatrixType,l>::check(_A_);

        const auto& space = mObj_.domainSpace();
        if( space.continuous() )
        {
          // only communicate diagonal if matrix blocks are small
          if constexpr ( size_t(MatrixBlockType::rows) == size_t(DiscreteFunctionSpaceType::dimRange) )
          {
            // create BlockVectorWrapper for BlockVector
            ISTLBlockVector< FlatMatrixBlock > diagonalInv( &diagonalInv_ );

            diagonalInv.resize(space.blockMapper().size());
            assert( space.blockMapper().size() == _A_.N() );
            diagonalInv.clear();

            // extract diagonal elements from matrix object
            {
              const auto end = _A_.end();
              for( auto row = _A_.begin(); row != end; ++row )
              {
                // get diagonal entry of matrix if existent
                auto diag = (*row).find( row.index() );
                if( diag != (*row).end() )
                {
                  MatrixConverterType m( diagonalInv[ row.index() ] );
                  m = (*diag);
                }
              }
            }

            // make diagonal consistent (communicate at border dofs)
            // get default communication operation type
            typename DiscreteFunctionSpaceType :: template
              CommDataHandle< DiscreteFunctionType > :: OperationType operation;
            space.communicate( diagonalInv, operation );

            // In general: store 1/diag if diagonal is number
            //
            // note: We set near-zero entries to 1 to avoid NaNs. Such entries occur
            //       if DoFs are excluded from matrix setup
            if constexpr( isNumber )
            {
              const double eps = 16.*std::numeric_limits< double >::epsilon();
              for( auto& diag : diagonalInv )
              {
                diag = (std::abs( diag ) < eps ? 1. : 1. / diag );
              }
            }
          }
          else
          {
            DUNE_THROW(NotImplemented,"ParallelIterative: communicating diagonal does not work for large block matrices");
          }
        }
      }

      ParallelIterative( const ParallelIterative& org )
        : mObj_( org.mObj_ ),
          _A_( org._A_ ),
          _n( org._n ), _w( org._w ),
          diagonalInv_( org.diagonalInv_ )
      {
      }

      //! \copydoc Preconditioner
      void pre (X& x, Y& b) override {}

      //! \copydoc Preconditioner
      void apply (X& v, const Y& d) override
      {
        const auto& space = mObj_.domainSpace();
        const bool continuous = space.continuous();

        // for communication
        DiscreteFunctionType tmp("ParIt::apply", space, v );

        static constexpr bool jacobi = (method == SolverParameter::jacobi);

        std::unique_ptr< X > xTmp;
        if constexpr ( jacobi )
          xTmp.reset( new X(v) );

        X& x = (jacobi) ? (*xTmp) : v;

        for (int i=0; i<_n; ++i)
        {
          forwardIterate(_A_, v, x, d, _w, diagonalInv_ );

          if constexpr ( method == SolverParameter::ssor )
          {
            // seems not to be necessary
            // communicate border unknowns
            //tmp.communicate();

            // create symmetry by iterating backwards
            backwardIterate(_A_, v, x, d, _w, diagonalInv_ );
          }

          // communicate border unknowns
          tmp.communicate();

          // for Jacobi now update the result
          if constexpr ( jacobi )
          {
            x = v;
            // for Jacobi skip last communication since this
            // is already in consistent state at this point
            if( continuous && i == _n-1)
              continue ;
          }
        }
      }

      //! \copydoc Preconditioner
      void post (X& x) override {}

      //! \brief The category the preconditioner is part of.
      SolverCategory::Category category () const override { return SolverCategory::sequential; }
    };

    template< class MatrixObject, class X, class Y >
    using ParallelSOR = ParallelIterative< MatrixObject, X, Y, SolverParameter::sor>;

    template< class MatrixObject, class X, class Y >
    using ParallelSSOR = ParallelIterative< MatrixObject, X, Y, SolverParameter::ssor>;

    template< class MatrixObject, class X, class Y >
    using ParallelJacobi = ParallelIterative< MatrixObject, X, Y, SolverParameter::jacobi>;

    template< class X, class Y >
    class IdentityPreconditionerWrapper : public Preconditioner<X,Y>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::field_type field_type;

      //! default constructor
      IdentityPreconditionerWrapper(){}

      //! \copydoc Preconditioner
      void pre (X& x, Y& b) override {}

      //! \copydoc Preconditioner
      void apply (X& v, const Y& d) override
      {
        if constexpr ( std::is_same< X, Y> :: value )
        {
          v = d;
        }
      }

      //! \copydoc Preconditioner
      void post (X& x) override {}

      SolverCategory::Category category () const override { return SolverCategory::sequential; }
    };


    //! wrapper class to store perconditioner
    //! as the interface class does not have to category
    //! enum
    template<class MatrixImp>
    class PreconditionerWrapper
      : public Preconditioner<typename MatrixImp :: RowBlockVectorType,
                              typename MatrixImp :: ColBlockVectorType>
    {
      typedef MatrixImp MatrixType;
      typedef typename MatrixType :: RowBlockVectorType X;
      typedef typename MatrixType :: ColBlockVectorType Y;

      // use BCRSMatrix type because of specializations in dune-istl
      typedef typename MatrixType :: BaseType ISTLMatrixType ;

      typedef typename MatrixType :: CommunicationType
        CommunicationType;

      // matrix adapter for AMG
  //#if HAVE_MPI
  //    typedef Dune::OverlappingSchwarzOperator<
  //     ISTLMatrixType, X, Y, CommunicationType> OperatorType ;
  //#else
      typedef MatrixAdapter< ISTLMatrixType, X, Y > OperatorType;
  //#endif
      mutable std::shared_ptr< OperatorType > op_;

      // auto pointer to preconditioning object
      typedef Preconditioner<X,Y> PreconditionerInterfaceType;
      mutable std::shared_ptr<PreconditionerInterfaceType> preconder_;

      // flag whether we have preconditioning, and if yes if it is AMG
      const int preEx_;
      const bool verbose_;

      template <class XImp, class YImp>
      struct Apply
      {
        inline static void apply(PreconditionerInterfaceType& preconder,
                                 XImp& v, const YImp& d)
        {
        }

        inline static void copy(XImp& v, const YImp& d)
        {
        }
     };

      template <class XImp>
      struct Apply<XImp,XImp>
      {
        inline static void apply(PreconditionerInterfaceType& preconder,
                                 XImp& v, const XImp& d)
        {
          preconder.apply( v ,d );
        }

        inline static void copy(X& v, const Y& d)
        {
          v = d;
        }
      };
    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::field_type field_type;

      //! copy constructor
      PreconditionerWrapper (const PreconditionerWrapper& org, bool verbose)
        : op_( org.op_ )
        , preconder_(org.preconder_)
        , preEx_(org.preEx_)
        , verbose_( verbose )
      {
      }

      //! default constructor
      PreconditionerWrapper(bool verbose)
        : op_()
        , preconder_()
        , preEx_( 0 )
        , verbose_( verbose )
      {}

      //! create preconditioner of given type
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix,
                            int iter,
                            field_type relax,
                            bool verbose,
                            const PreconditionerType* p)
        : op_()
        , preconder_( new PreconditionerType( matrix, iter, relax ) )
        , preEx_( 1 )
        , verbose_( verbose )
      {
      }


      //! create preconditioner with given preconditioner object
      //! owner ship is taken over here
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix, bool verbose,
                            PreconditionerType* p)
        : op_()
        , preconder_( p )
        , preEx_( 1 )
        , verbose_( verbose )
      {
      }

      //! create preconditioner of given type
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix,
                            int iter,
                            field_type relax,
                            bool verbose,
                            const PreconditionerType* p ,
                            const CommunicationType& comm)
  //#if HAVE_MPI
  //      : op_( new OperatorType( matrix, comm ) )
  //#else
        : op_( new OperatorType( matrix ) )
  //#endif
        , preconder_( createAMGPreconditioner(comm, iter, relax, p) )
        , preEx_( 2 )
        , verbose_( verbose )
      {
      }

      //! \copydoc Preconditioner
      void pre (X& x, Y& b) override
      {
        // all the sequentiel implemented Preconditioners do nothing in pre and post
        // apply preconditioner
        if( preEx_ > 1 )
        {
          //X tmp (x);
          preconder_->pre(x,b);
          //assert( std::abs( x.two_norm() - tmp.two_norm() ) < 1e-15);
        }
      }

      //! \copydoc Preconditioner
      void apply (X& v, const Y& d) override
      {
        if( preEx_ )
        {
          // apply preconditioner
          Apply<X,Y>::apply( *preconder_ , v, d );
        }
        else
        {
          Apply<X,Y>::copy( v, d );
        }
      }

      //! \copydoc Preconditioner
      void post (X& x) override
      {
        // all the implemented Preconditioners do nothing in pre and post
        // apply preconditioner
        if( preEx_ > 1 )
        {
          preconder_->post(x);
        }
      }

      //! \brief The category the precondtioner is part of.
      SolverCategory::Category category () const override
      {
        return (preconder_ ? preconder_->category() : SolverCategory::sequential);
      }

    protected:
      template <class Smoother>
      PreconditionerInterfaceType*
      createAMGPreconditioner(const CommunicationType& comm,
                int iter, field_type relax, const Smoother* )
      {
        typedef typename Dune::FieldTraits< field_type>::real_type real_type;
        typedef typename std::conditional< std::is_convertible<field_type,real_type>::value,
                         Dune::Amg::FirstDiagonal, Dune::Amg::RowSum >::type Norm;
        typedef Dune::Amg::CoarsenCriterion<
                Dune::Amg::UnSymmetricCriterion<ISTLMatrixType, Norm > > Criterion;

        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;

        SmootherArgs smootherArgs;

        smootherArgs.iterations = iter;
        smootherArgs.relaxationFactor = relax ;

        int coarsenTarget=1200;
        Criterion criterion(15,coarsenTarget);
        criterion.setDefaultValuesIsotropic(2);
        criterion.setAlpha(.67);
        criterion.setBeta(1.0e-8);
        criterion.setMaxLevel(10);
        if( verbose_ && Parameter :: verbose() )
          criterion.setDebugLevel( 1 );
        else
          criterion.setDebugLevel( 0 );

        /*
        if( comm.size() > 1 )
        {
          typedef Dune::OwnerOverlapCopyCommunication<int> ParallelInformation;
          ParallelInformation pinfo(MPI_COMM_WORLD);
          typedef Dune::Amg::AMG<OperatorType, X, Smoother, ParallelInformation> AMG;
          return new AMG(*op_, criterion, smootherArgs, pinfo);
        }
        else
        */
        {
          // X == Y is needed for AMG
          typedef Dune::Amg::AMG<OperatorType, X, Smoother> AMG;
          return new AMG(*op_, criterion, smootherArgs);
          //return new AMG(*op_, criterion, smootherArgs, 1, 1, 1, false);
        }
      }
    };


    // ISTLPreconditionerFactory
    // -------------------------

    template < class MatrixObject >
    class ISTLMatrixAdapterFactory ;

    template < class DomainSpace, class RangeSpace, class DomainBlock, class RangeBlock,
               template <class, class, class, class> class MatrixObject >
    class ISTLMatrixAdapterFactory< MatrixObject< DomainSpace, RangeSpace, DomainBlock, RangeBlock > >
    {
      typedef DomainSpace       DomainSpaceType ;
      typedef RangeSpace        RangeSpaceType;

    public:
      typedef MatrixObject< DomainSpaceType, RangeSpaceType, DomainBlock, RangeBlock >  MatrixObjectType;
      typedef typename MatrixObjectType :: MatrixType            MatrixType;
      typedef ISTLParallelMatrixAdapterInterface< MatrixType >   MatrixAdapterInterfaceType;

      //! return matrix adapter object that works with ISTL linear solvers
      static std::unique_ptr< MatrixAdapterInterfaceType >
      matrixAdapter( const MatrixObjectType& matrixObj,
                     const ISTLSolverParameter& param)
      {
        std::unique_ptr< MatrixAdapterInterfaceType > ptr;
        if( matrixObj.domainSpace().continuous() )
        {
          typedef LagrangeParallelMatrixAdapter< MatrixType > MatrixAdapterImplementation;
          ptr.reset( matrixAdapterObject( matrixObj, (MatrixAdapterImplementation *) nullptr, param ) );
        }
        else
        {
          typedef DGParallelMatrixAdapter< MatrixType > MatrixAdapterImplementation;
          ptr.reset( matrixAdapterObject( matrixObj, (MatrixAdapterImplementation *) nullptr, param ) );
        }
        return ptr;
      }

    protected:
      template <class MatrixAdapterType>
      static MatrixAdapterType*
      matrixAdapterObject( const MatrixObjectType& matrixObj,
                           const MatrixAdapterType*,
                           const ISTLSolverParameter& param )
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        return new MatrixAdapterType( matrixObj.exportMatrix(),
                                      matrixObj.domainSpace(), matrixObj.rangeSpace(), PreConType(param.verbose()),
                                      matrixObj.threading() );
      }

    };

#ifndef DISABLE_ISTL_PRECONDITIONING
    //! Specialization for domain space == range space
    template < class Space, class DomainBlock, class RangeBlock,
               template <class, class, class, class> class MatrixObject >
    class ISTLMatrixAdapterFactory< MatrixObject< Space, Space, DomainBlock, RangeBlock > >
    {
    public:
      typedef Space       DomainSpaceType ;
      typedef Space       RangeSpaceType;

      typedef MatrixObject< DomainSpaceType, RangeSpaceType, DomainBlock, RangeBlock >  MatrixObjectType;
      typedef typename MatrixObjectType :: MatrixType              MatrixType;
      typedef typename MatrixType   :: BaseType                    ISTLMatrixType ;
      typedef Fem::PreconditionerWrapper< MatrixType >             PreconditionAdapterType;

    protected:
      template <class MatrixAdapterType, class PreconditionerType>
      static MatrixAdapterType*
      createMatrixAdapter(const MatrixAdapterType*,
                          const PreconditionerType* preconditioning,
                          MatrixType& matrix,
                          const DomainSpaceType& domainSpace,
                          const RangeSpaceType& rangeSpace,
                          const double relaxFactor,
                          std::size_t numIterations,
                          bool verbose)
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix, numIterations, relaxFactor, verbose, preconditioning );
        return new MatrixAdapterType(matrix, domainSpace, rangeSpace, preconAdapter );
      }

      template <class MatrixAdapterType, class PreconditionerType>
      static MatrixAdapterType*
      createAMGMatrixAdapter(const MatrixAdapterType*,
                             const PreconditionerType* preconditioning,
                             MatrixType& matrix,
                             const DomainSpaceType& domainSpace,
                             const RangeSpaceType& rangeSpace,
                             const double relaxFactor,
                             std::size_t numIterations,
                             bool solverVerbose)
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix, numIterations, relaxFactor, solverVerbose,
            preconditioning, domainSpace.gridPart().comm() );
        return new MatrixAdapterType(matrix, domainSpace, rangeSpace, preconAdapter );
      }

      template <class MatrixAdapterType>
      static MatrixAdapterType*
      matrixAdapterObject( const MatrixObjectType& matrixObj,
                           const MatrixAdapterType*,
                           const ISTLSolverParameter& param )
      {
        int preconditioning = param.preconditionMethod(
            ISTLPreconditionMethods::supportedPreconditionMethods() );
        const double relaxFactor   = param.relaxation();
        const size_t numIterations = param.preconditionerIteration();

        const DomainSpaceType& domainSpace = matrixObj.domainSpace();
        const RangeSpaceType&  rangeSpace  = matrixObj.rangeSpace();

        MatrixType& matrix = matrixObj.exportMatrix();
        const auto procs = domainSpace.gridPart().comm().size();

        typedef typename MatrixType :: BaseType ISTLMatrixType;
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;

        typedef typename MatrixObjectType :: RowBlockVectorType      RowBlockVectorType;
        typedef typename MatrixObjectType :: ColumnBlockVectorType   ColumnBlockVectorType;

        // no preconditioner
        if( preconditioning == SolverParameter::none )
        {
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, PreConType(param.verbose()) );
        }
        // SSOR
        else if( preconditioning == SolverParameter::ssor )
        {
          typedef ParallelSSOR<MatrixObjectType, RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrixObj, numIterations, relaxFactor ) );
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
        }
        // SOR and Gauss-Seidel
        else if(preconditioning == SolverParameter::sor || preconditioning == SolverParameter::gauss_seidel)
        {
          typedef ParallelSOR<MatrixObjectType, RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrixObj, numIterations,
                                    (preconditioning == SolverParameter::gauss_seidel) ? 1.0 : relaxFactor ) );
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
        }
        // ILU
        else if(preconditioning == SolverParameter::ilu)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILU not working in parallel computations");

          typedef SeqILU<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          // might need to set verbosity here?
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrix, numIterations, relaxFactor, param.fastILUStorage()) );
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
        }
        // Jacobi
        else if(preconditioning == SolverParameter::jacobi)
        {
          typedef ParallelJacobi<MatrixObjectType, RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrixObj, numIterations, relaxFactor ) );
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
        }
        // AMG ILU-0
        else if(preconditioning == SolverParameter::amg_ilu)
        {
          // use original SeqILU because of some AMG traits classes.
          typedef SeqILU<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createAMGMatrixAdapter( (MatrixAdapterType *)nullptr,
                                         (PreconditionerType*)nullptr,
                                         matrix, domainSpace, rangeSpace, relaxFactor, numIterations,
                                         param.verbose());
        }
        // AMG Jacobi
        else if(preconditioning == SolverParameter::amg_jacobi)
        {
          typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType,1> PreconditionerType;
          return createAMGMatrixAdapter( (MatrixAdapterType *)nullptr,
                                         (PreconditionerType*)nullptr,
                                         matrix, domainSpace, rangeSpace, relaxFactor, numIterations,
                                         param.verbose());
        }
        // ILDL
        else if(preconditioning == SolverParameter::ildl)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILDL not working in parallel computations");

          PreConType preconAdapter( matrix, param.verbose(), new SeqILDL<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType>( matrix , relaxFactor ) );
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
        }
        else
        {
          preConErrorMsg(preconditioning);
        }
        return new MatrixAdapterType(matrix, domainSpace, rangeSpace,  PreConType(param.verbose()) );
      }

      typedef ISTLParallelMatrixAdapterInterface< MatrixType >   MatrixAdapterInterfaceType;

      static void preConErrorMsg(int preCon)
      {
        // TODO: re-write
        std::cerr << "ERROR: Wrong precoditioning number (p = " << preCon;
        std::cerr <<") in ISTLMatrixObject! \n";
        std::cerr <<"Valid values are: \n";
        std::cerr <<"0 == no \n";
        std::cerr <<"1 == SSOR \n";
        std::cerr <<"2 == SOR \n";
        std::cerr <<"3 == ILU-0 \n";
        std::cerr <<"4 == ILU-n \n";
        std::cerr <<"5 == Gauss-Seidel \n";
        std::cerr <<"6 == Jacobi \n";

        assert( false );
        DUNE_THROW(NotImplemented,"Wrong precoditioning selected");
      }

    public:
      //! return matrix adapter object that works with ISTL linear solvers
      static std::unique_ptr< MatrixAdapterInterfaceType >
      matrixAdapter( const MatrixObjectType& matrixObj,
                     const ISTLSolverParameter& param)
      {
        const ISTLSolverParameter* parameter = dynamic_cast< const ISTLSolverParameter* > (&param);
        std::unique_ptr< ISTLSolverParameter > paramPtr;
        if( ! parameter )
        {
          paramPtr.reset( new ISTLSolverParameter( param ) );
          parameter = paramPtr.operator->();
        }

        std::unique_ptr< MatrixAdapterInterfaceType > ptr;
        if( matrixObj.domainSpace().continuous() )
        {
          typedef LagrangeParallelMatrixAdapter< MatrixType > MatrixAdapterImplementation;
          ptr.reset( matrixAdapterObject( matrixObj, (MatrixAdapterImplementation *) nullptr, *parameter ) );
        }
        else
        {
          typedef DGParallelMatrixAdapter< MatrixType > MatrixAdapterImplementation;
          ptr.reset( matrixAdapterObject( matrixObj, (MatrixAdapterImplementation *) nullptr, *parameter ) );
        }
        return ptr;
      }

    }; // end specialization of ISTLMatrixAdapterFactor for domainSpace == rangeSpace
#endif

#endif // end HAVE_DUNE_ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PRECONDITIONERWRAPPER_HH
