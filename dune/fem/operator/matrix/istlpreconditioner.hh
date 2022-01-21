#ifndef DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH
#define DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH

#include <memory>

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

namespace Dune
{

  namespace Fem
  {

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
    /** \brief Iterative linear methods as preconditioners.
     *         Note: For Lagrange type spaces this only works moderately, since
     *               the communication of the diagonal is not yet completely done,
     *               only in the scalar case.
     *         Note: On rows corresponding to border DoFs this scheme falls back to Jacobi, to avoid
     *               complicated coloring strategies.
     */
    template< class MatrixObject, class X, class Y, bool jacobi, int l=1 >
    class ParallelIterative : public Preconditioner<X,Y>
    {
    public:
      typedef typename MatrixObject :: ColumnDiscreteFunctionType DiscreteFunctionType ;

    protected:
      typedef typename X::field_type field_type;
      typedef typename MatrixObject::MatrixType  MatrixType;
      typedef MatrixType M;
      typedef typename M::block_type MatrixBlockType;

      static const bool isNumber = (M::block_type::rows == M::block_type::cols) == 1;

      // implementation of SOR and Jacobi's method
      void iterate (const M& A, X& x, const X& xold, const Y& b, const field_type& w,
                    const std::unique_ptr< DiscreteFunctionType >& diagInv) const
      {
        typedef typename M::ConstRowIterator rowiterator;
        typedef typename M::ConstColIterator coliterator;
        typedef typename Y::block_type bblock;
        typedef typename X::block_type xblock;

        bblock rhs;
        xblock v;

        // Initialize nested data structure if there are entries
        if(A.begin()!=A.end())
          v=x[0];

        const rowiterator endi=A.end();
        for (rowiterator i=A.begin(); i!=endi; ++i)
        {
          const coliterator endj=(*i).end();           // iterate over a_ij with j < i
          coliterator j=(*i).begin();
          const auto rowIdx = i.index();
          // treat missing rows as unit rows
          if( j == endj )
          {
            x[rowIdx] += w*b[rowIdx];
            continue ;
          }

          rhs = b[rowIdx];           // rhs = b_i
          if constexpr (isNumber)
          {
            for (; j.index()<rowIdx; ++j)
              rhs -= (*j) * x[j.index()];  //  rhs -= sum_{j<i} a_ij * xnew_j
            // not needed, since we store the diagonal separately
            coliterator diag=j;            // *diag = a_ii
            for (; j!=endj; ++j)
              rhs -= (*j) * x[j.index()];  //  rhs -= sum_{j<i} a_ij * xnew_j

            // v = rhs / diag
            if( diagInv )
              v = rhs * diagInv->dofVector()[ rowIdx ];
            else
              v = rhs / (*diag);

            x[ rowIdx ] += w*v;           //  x_i = w / a_ii * (b_i - sum_{j<i} a_ij * xnew_j - sum_{j>=i} a_ij * xold_j)
          }
          else
          {
            for (; j.index()<rowIdx; ++j)
              (*j).mmv(x[j.index()],rhs);  //  rhs -= sum_{j<i} a_ij * xnew_j
            coliterator diag=j;            // *diag = a_ii
            for (; j!=endj; ++j)
              (*j).mmv(x[j.index()],rhs);  //  rhs -= sum_{j<i} a_ij * xnew_j

            (*diag).solve(v, rhs);
            x[rowIdx].axpy(w,v);        //  x_i = w / a_ii * (b_i - sum_{j<i} a_ij * xnew_j - sum_{j>=i} a_ij * xold_j)
          }
        }

      }

    protected:
      const MatrixObject& mObj_;
      const MatrixType& _A_;
      const int _n;
      const double _w;
      std::unique_ptr< DiscreteFunctionType > diagonalInv_;

    public:
      ParallelIterative( const MatrixObject& mObj, int n=1, double relax=1.0  )
        : mObj_( mObj ),
          _A_( mObj.exportMatrix() ),
          _n( n ), _w(relax)
      {
        Dune::CheckIfDiagonalPresent<MatrixType,l>::check(_A_);

        if( mObj_.domainSpace().continuous() && isNumber )
        {
          diagonalInv_.reset( new DiscreteFunctionType("ParIt::diag", mObj_.domainSpace() ) );
          // extract diagonal elements form matrix object
          mObj_.extractDiagonal( *diagonalInv_ );

          // make consistent at border dofs
          diagonalInv_->communicate();

          // In general: store 1/diag
          //
          // note: We set near-zero entries to 1 to avoid NaNs. Such entries occur
          //       if DoFs are excluded from matrix setup
          const double eps = 16.*std::numeric_limits< double >::epsilon();
          const auto dend = diagonalInv_->dend();
          for( auto dit = diagonalInv_->dbegin(); dit != dend; ++dit )
            *dit = (std::abs( *dit ) < eps ? 1. : 1. / *dit);
        }
        else if (mObj_.domainSpace().continuous())
        {
          DUNE_THROW(NotImplemented,"ParallelIterative preconditioner needs communication of diagonal block with blockSize > 1");
        }
      }

      //! \copydoc Preconditioner
      void pre (X& x, Y& b) override {}

      //! \copydoc Preconditioner
      void apply (X& v, const Y& d) override
      {
        const auto& space = mObj_.domainSpace();
        // for communication
        DiscreteFunctionType tmp("ParIt::tmp", space, v );

        std::unique_ptr< X > xTmp;
        if constexpr ( jacobi )
          xTmp.reset( new X(v) );

        X& x = (jacobi) ? (*xTmp) : v;

        for (int i=0; i<_n; ++i)
        {
          iterate(_A_, x, v, d, _w, diagonalInv_ );
          if constexpr ( jacobi )
            v = x;
          // communicate border unknowns
          tmp.communicate();
        }
      }

      //! \copydoc Preconditioner
      void post (X& x) override {}

      //! \brief The category the preconditioner is part of.
      SolverCategory::Category category () const override { return SolverCategory::sequential; }
    };

    template< class MatrixObject, class X, class Y >
    using ParallelSOR = ParallelIterative< MatrixObject, X, Y, false>;

    template< class MatrixObject, class X, class Y >
    using ParallelJacobi = ParallelIterative< MatrixObject, X, Y, true>;

    template< class X, class Y >
    class IdentityPreconditionerWrapper : public Preconditioner<X,Y>
    {
      template <class XImp, class YImp>
      struct Apply
      {
        inline static void copy(XImp& v, const YImp& d)
        {
        }
      };

      template <class XImp>
      struct Apply<XImp,XImp>
      {
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

      //! default constructor
      IdentityPreconditionerWrapper(){}

      //! \copydoc Preconditioner
      void pre (X& x, Y& b) override {}

      //! \copydoc Preconditioner
      void apply (X& v, const Y& d) override
      {
        Apply< X, Y> :: copy( v, d );
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

      typedef typename MatrixType :: CollectiveCommunictionType
        CollectiveCommunictionType;

      // matrix adapter for AMG
  //#if HAVE_MPI
  //    typedef Dune::OverlappingSchwarzOperator<
  //     ISTLMatrixType, X, Y, CollectiveCommunictionType> OperatorType ;
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
                            const CollectiveCommunictionType& comm)
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
      createAMGPreconditioner(const CollectiveCommunictionType& comm,
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
                                      matrixObj.domainSpace(), matrixObj.rangeSpace(), PreConType(param.verbose()) );
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
                  { SolverParameter::none,        // no preconditioner
                    SolverParameter::ssor,        // SSOR preconditioner
                    SolverParameter::sor ,        // SOR preconditioner
                    SolverParameter::ilu ,        // ILU preconditioner (deprecated)
                    SolverParameter::gauss_seidel,// Gauss-Seidel preconditioner
                    SolverParameter::jacobi,      // Jacobi preconditioner
                    SolverParameter::amg_ilu,     // AMG with ILU-0 smoother (deprecated)
                    SolverParameter::amg_jacobi,  // AMG with Jacobi smoother
                    SolverParameter::ildl         // ILDL from istl
                  } );
        const double relaxFactor         = param.relaxation();
        const size_t numIterations       = param.preconditionerIteration();

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
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqSSOR not working in parallel computations");

          typedef SeqSSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations, param.verbose() );
        }
        // SOR
        else if(preconditioning == SolverParameter::sor )
        {
          typedef ParallelSOR<MatrixObjectType, RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrixObj, numIterations, relaxFactor ) );
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
        // Gauss-Seidel
        else if(preconditioning == SolverParameter::gauss_seidel)
        {
          typedef ParallelSOR<MatrixObjectType, RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix, param.verbose(), new PreconditionerType( matrixObj, numIterations, 1.0 ) );
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
