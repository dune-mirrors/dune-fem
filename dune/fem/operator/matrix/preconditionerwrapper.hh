#ifndef DUNE_FEM_PRECONDITIONERWRAPPER_HH
#define DUNE_FEM_PRECONDITIONERWRAPPER_HH

#include <memory>
#include <dune/common/shared_ptr.hh>

// standard diagonal preconditioner
#include <dune/fem/solver/diagonalpreconditioner.hh>


#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#endif

namespace Dune
{

  namespace Fem
  {

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
        typedef Dune::Amg::CoarsenCriterion<
          Dune::Amg::UnSymmetricCriterion<ISTLMatrixType,
                                    Dune::Amg::FirstDiagonal> > Criterion;

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
#endif // end HAVE_DUNE_ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PRECONDITIONERWRAPPER_HH
