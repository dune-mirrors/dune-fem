#ifndef DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH
#define DUNE_FEM_ISTLPRECONDITIONERWRAPPER_HH

#include <memory>
// standard diagonal preconditioner
#include <dune/fem/solver/diagonalpreconditioner.hh>


#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#endif

#include <dune/fem/operator/matrix/spmatrix.hh> // MatrixParameter
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

    struct ISTLMatrixParameter
      : public MatrixParameter
    {

      ISTLMatrixParameter( const std::string keyPrefix = "istl." )
        : keyPrefix_( keyPrefix )
      {}

      virtual double overflowFraction () const
      {
        return Parameter::getValue< double >( keyPrefix_ + "matrix.overflowfraction", 1.0 );
      }

      virtual int numIterations () const
      {
        return Parameter::getValue< int >( keyPrefix_ + "preconditioning.iterations", 5 );
      }

      virtual double relaxation () const
      {
        return Parameter::getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
      }

      virtual int method () const
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu-0", "ilu-n", "gauss-seidel", "jacobi", "amg-ilu-0", "amg-ilu-n", "amg-jacobi" };
        return Parameter::getEnum(  keyPrefix_ + "preconditioning.method", preConTable, 0 );
      }

      virtual std::string preconditionName() const
      {
        static const std::string preConTable[]
          = { "None", "SSOR", "SOR", "ILU-0", "ILU-n", "Gauss-Seidel", "Jacobi", "AMG-ILU-0", "AMG-ILU-n", "AMG-Jacobi" };
        const int precond = method();
        std::stringstream tmp;
        tmp << preConTable[precond];

        if( precond != 3 )
          tmp << " n=" << numIterations();
        tmp << " relax=" << relaxation();
        return tmp.str();
      }

     private:
      std::string keyPrefix_;

    };

#if HAVE_DUNE_ISTL
    template<class M, class X, class Y, int l=1>
    class FemSeqILU0 : public SeqILU0<M,X,Y,l>
    {
      typedef SeqILU0<M,X,Y,l>  BaseType ;
    public:
      typedef typename X::field_type field_type;
      FemSeqILU0 (const M& A, int iter, field_type w)
        : BaseType( A, w )
      {}
      FemSeqILU0 (const M& A, field_type w)
        : BaseType( A, w )
      {}
    };


    template< class MatrixObject, class X, class Y >
    class FemDiagonalPreconditioner : public Preconditioner<X,Y>
    {
    public:
      typedef typename MatrixObject :: ColumnDiscreteFunctionType DiscreteFunctionType ;
      // select preconditioner for assembled operators
      typedef DiagonalPreconditionerBase< DiscreteFunctionType, MatrixObject, true > PreconditionerType ;

    protected:
      PreconditionerType diagonalPrecon_;

    public:
      typedef typename X::field_type field_type;
      FemDiagonalPreconditioner( const MatrixObject& mObj )
        : diagonalPrecon_( mObj )
      {}

      //! \copydoc Preconditioner
      virtual void pre (X& x, Y& b) {}

      //! \copydoc Preconditioner
      virtual void apply (X& v, const Y& d)
      {
        diagonalPrecon_.applyToISTLBlockVector( d, v );
      }

      //! \copydoc Preconditioner
      virtual void post (X& x) {}
    };


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

      enum {
        //! \brief The category the precondtioner is part of.
        category=SolverCategory::sequential };

      //! default constructor
      IdentityPreconditionerWrapper(){}

      //! \copydoc Preconditioner
      virtual void pre (X& x, Y& b) {}

      //! \copydoc Preconditioner
      virtual void apply (X& v, const Y& d)
      {
        Apply< X, Y> :: copy( v, d );
      }

      //! \copydoc Preconditioner
      virtual void post (X& x) {}
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

      enum {
        //! \brief The category the precondtioner is part of.
        category=SolverCategory::sequential };

      //! copy constructor
      PreconditionerWrapper (const PreconditionerWrapper& org)
        : op_( org.op_ )
        , preconder_(org.preconder_)
        , preEx_(org.preEx_)
      {
      }

      //! default constructor
      PreconditionerWrapper()
        : op_()
        , preconder_()
        , preEx_( 0 )
      {}

      //! create preconditioner of given type
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix,
                            int iter,
                            field_type relax,
                            const PreconditionerType* p )
        : op_()
        , preconder_( new PreconditionerType( matrix, iter, relax ) )
        , preEx_( 1 )
      {
      }


      //! create preconditioner with given preconditioner object
      //! owner ship is taken over here
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix,
                            PreconditionerType* p )
        : op_()
        , preconder_( p )
        , preEx_( 1 )
      {
      }

      //! create preconditioner of given type
      template <class PreconditionerType>
      PreconditionerWrapper(MatrixType & matrix,
                            int iter,
                            field_type relax,
                            const PreconditionerType* p ,
                            const CollectiveCommunictionType& comm )
  //#if HAVE_MPI
  //      : op_( new OperatorType( matrix, comm ) )
  //#else
        : op_( new OperatorType( matrix ) )
  //#endif
        , preconder_( createAMGPreconditioner(comm, iter, relax, p) )
        , preEx_( 2 )
      {
      }

      //! \copydoc Preconditioner
      virtual void pre (X& x, Y& b)
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
      virtual void apply (X& v, const Y& d)
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
      virtual void post (X& x)
      {
        // all the implemented Preconditioners do nothing in pre and post
        // apply preconditioner
        if( preEx_ > 1 )
        {
          preconder_->post(x);
        }
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

    template < class DomainSpace, class RangeSpace,
               template <class, class> class MatrixObject >
    class ISTLMatrixAdapterFactory< MatrixObject< DomainSpace, RangeSpace > >
    {
      typedef DomainSpace       DomainSpaceType ;
      typedef RangeSpace        RangeSpaceType;

    public:
      typedef MatrixObject< DomainSpaceType, RangeSpaceType >    MatrixObjectType;
      typedef typename MatrixObjectType :: MatrixType            MatrixType;
      typedef ISTLParallelMatrixAdapterInterface< MatrixType >   MatrixAdapterInterfaceType;

      //! return matrix adapter object that works with ISTL linear solvers
      static std::unique_ptr< MatrixAdapterInterfaceType >
      matrixAdapter( const MatrixObjectType& matrixObj,
                     const MatrixParameter& param = ISTLMatrixParameter() )
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
                           const MatrixParameter& param )
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        return new MatrixAdapterType( matrixObj.matrix(),
                                      matrixObj.domainSpace(), matrixObj.rangeSpace(), PreConType() );
      }

    };

#ifndef DISABLE_ISTL_PRECONDITIONING
    //! Specialization for domain space == range space
    template < class Space,
               template <class, class> class MatrixObject >
    class ISTLMatrixAdapterFactory< MatrixObject< Space, Space > >
    {
    public:
      enum ISTLPreConder_Id { none  = 0 ,      // no preconditioner
                              ssor  = 1 ,      // SSOR preconditioner
                              sor   = 2 ,      // SOR preconditioner
                              ilu_0 = 3 ,      // ILU-0 preconditioner
                              ilu_n = 4 ,      // ILU-n preconditioner
                              gauss_seidel= 5, // Gauss-Seidel preconditioner
                              jacobi = 6,      // Jacobi preconditioner
                              amg_ilu_0 = 7,   // AMG with ILU-0 smoother
                              amg_ilu_n = 8,   // AMG with ILU-n smoother
                              amg_jacobi = 9   // AMG with Jacobi smoother
      };

      typedef Space       DomainSpaceType ;
      typedef Space       RangeSpaceType;

      typedef MatrixObject< DomainSpaceType, RangeSpaceType >  MatrixObjectType;
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
                          std::size_t numIterations)
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix, numIterations, relaxFactor, preconditioning );
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
                             std::size_t numIterations)
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix, numIterations, relaxFactor, preconditioning, domainSpace.gridPart().comm() );
        return new MatrixAdapterType(matrix, domainSpace, rangeSpace, preconAdapter );
      }

      template <class MatrixAdapterType>
      static MatrixAdapterType*
      matrixAdapterObject( const MatrixObjectType& matrixObj,
                           const MatrixAdapterType*,
                           const MatrixParameter& param )
      {
        ISTLPreConder_Id preconditioning = (ISTLPreConder_Id)param.method() ;
        const double relaxFactor         = param.relaxation();
        const size_t numIterations       = param.numIterations();

        const DomainSpaceType& domainSpace = matrixObj.domainSpace();
        const RangeSpaceType&  rangeSpace  = matrixObj.rangeSpace();

        MatrixType& matrix = matrixObj.matrix();
        const auto procs = domainSpace.gridPart().comm().size();

        typedef typename MatrixType :: BaseType ISTLMatrixType;
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;

        typedef typename MatrixObjectType :: RowBlockVectorType      RowBlockVectorType;
        typedef typename MatrixObjectType :: ColumnBlockVectorType   ColumnBlockVectorType;

        // no preconditioner
        if( preconditioning == none )
        {
          return new MatrixAdapterType( matrix, domainSpace, rangeSpace, PreConType() );
        }
        // SSOR
        else if( preconditioning == ssor )
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqSSOR not working in parallel computations");

          typedef SeqSSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // SOR
        else if(preconditioning == sor )
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqSOR not working in parallel computations");

          typedef SeqSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // ILU-0
        else if(preconditioning == ilu_0)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILU0 not working in parallel computations");

          typedef FemSeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // ILU-n
        else if(preconditioning == ilu_n)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILUn not working in parallel computations");

          typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // Gauss-Seidel
        else if(preconditioning == gauss_seidel)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqGS not working in parallel computations");

          typedef SeqGS<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                      (PreconditionerType*)nullptr,
                                      matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // Jacobi
        else if(preconditioning == jacobi)
        {
          if( numIterations == 1 ) // diagonal preconditioning
          {
            typedef FemDiagonalPreconditioner< MatrixObjectType, RowBlockVectorType, ColumnBlockVectorType > PreconditionerType;
            typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
            PreConType preconAdapter( matrix, new PreconditionerType( matrixObj ) );
            return new MatrixAdapterType( matrix, domainSpace, rangeSpace, preconAdapter );
          }
          else if ( procs == 1 )
          {
            typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
            return createMatrixAdapter( (MatrixAdapterType *)nullptr,
                                        (PreconditionerType*)nullptr,
                                        matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
          }
          else
          {
            DUNE_THROW(InvalidStateException,"ISTL::SeqJac(Jacobi) only working with istl.preconditioning.iterations: 1 in parallel computations");
          }
        }
        // AMG ILU-0
        else if(preconditioning == amg_ilu_0)
        {
          // use original SeqILU0 because of some AMG traits classes.
          typedef SeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createAMGMatrixAdapter( (MatrixAdapterType *)nullptr,
                                         (PreconditionerType*)nullptr,
                                         matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // AMG ILU-n
        else if(preconditioning == amg_ilu_n)
        {
          typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createAMGMatrixAdapter( (MatrixAdapterType *)nullptr,
                                         (PreconditionerType*)nullptr,
                                         matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        // AMG Jacobi
        else if(preconditioning == amg_jacobi)
        {
          typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType,1> PreconditionerType;
          return createAMGMatrixAdapter( (MatrixAdapterType *)nullptr,
                                         (PreconditionerType*)nullptr,
                                         matrix, domainSpace, rangeSpace, relaxFactor, numIterations );
        }
        else
        {
          preConErrorMsg(preconditioning);
        }
        return new MatrixAdapterType(matrix, domainSpace, rangeSpace, PreConType() );
      }

      typedef ISTLParallelMatrixAdapterInterface< MatrixType >   MatrixAdapterInterfaceType;

      static void preConErrorMsg(int preCon)
      {
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
                     const MatrixParameter& param = ISTLMatrixParameter() )
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

    }; // end specialization of ISTLMatrixAdapterFactor for domainSpace == rangeSpace
#endif

#endif // end HAVE_DUNE_ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PRECONDITIONERWRAPPER_HH
