#ifndef DUNE_PRECONDITIONERWRAPPER_HH
#define DUNE_PRECONDITIONERWRAPPER_HH

#include <memory>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#endif

namespace Dune {

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
            
    typedef Preconditioner<X,Y> PreconditionerInterfaceType;
    MatrixType& matrix_;
    mutable std::auto_ptr<PreconditionerInterfaceType> preconder_; 
    const bool preEx_;

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
      : matrix_(org.matrix_) 
      , preconder_(org.preconder_) 
      , preEx_(org.preEx_)
    {
    }
    
    //! constructor for creating none-preconditioning wrapper 
    PreconditionerWrapper (MatrixType& m) 
      : matrix_(m) 
      , preconder_()
      , preEx_(false)  
    {}
    
    //! create preconditioner of given type 
    template <class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, 
                          int iter,
                          field_type relax, 
                          const PreconditionerType* p )
      : matrix_(m)
      , preconder_( new PreconditionerType( matrix_, iter, relax ) )
      , preEx_(true) 
    {
    }
    
    //! create preconditioner of given type 
    template <class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, 
                          int iter,
                          field_type relax, 
                          const PreconditionerType* p ,
                          const bool )
      : matrix_(m)
      , preconder_( createAMGPreconditioner(iter, relax, p) )
      , preEx_(true) 
    {
    }
    
    //! \copydoc Preconditioner 
    virtual void pre (X& x, Y& b) 
    {
      /*
      // all the implemented Preconditioners do nothing in pre and post 
#ifndef NDEBUG 
      // apply preconditioner
      if( preEx_ ) 
      {
        X tmp (x);
        preconder_->pre(x,b);
        assert( std::abs( x.two_norm() - tmp.two_norm() ) < 1e-15);
      }
#endif
      */
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
#ifndef NDEBUG 
      // apply preconditioner
      if( preEx_ ) 
      {
        X tmp(x);
        preconder_->post(x);
        assert( std::abs( x.two_norm() - tmp.two_norm() ) < 1e-15);
      }
#endif
    }

  protected:  
    template <class S>
    PreconditionerInterfaceType* 
    createAMGPreconditioner( int iter, field_type relax, const S* ) 
    {
      // use BCRSMatrix type because of specializations in dune-istl
      typedef typename MatrixType :: BaseType ISTLMatrixType ;

      typedef Dune::Amg::CoarsenCriterion<
        Dune::Amg::UnSymmetricCriterion<ISTLMatrixType,
                                  Dune::Amg::FirstDiagonal> > Criterion;
      
      // preconditioner 
      typedef S Smoother ;
      //typedef SeqILUn< ISTLMatrixType, X, Y > Smoother;

      typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;

      SmootherArgs smootherArgs;

      smootherArgs.iterations = iter;
      smootherArgs.relaxationFactor = relax ;

      // matrix adapter 
      typedef MatrixAdapter< ISTLMatrixType, X, Y > OperatorType;
      OperatorType op( matrix_ );

      int coarsenTarget=1200;
      Criterion criterion(15,coarsenTarget);
      criterion.setDefaultValuesIsotropic(2);
      criterion.setAlpha(.67);
      criterion.setBeta(1.0e-8);
      criterion.setMaxLevel(10);

      // X or Y ???
      typedef Dune::Amg::AMG<OperatorType, X, Smoother> AMG;
      return new AMG(op, criterion, smootherArgs, 1, 1, 1, false);
    }
  };
#endif // end HAVE_DUNE_ISTL 

} // end namespace Dune 
#endif
