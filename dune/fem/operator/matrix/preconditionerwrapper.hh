#ifndef DUNE_PRECONDITIONERWRAPPER_HH
#define DUNE_PRECONDITIONERWRAPPER_HH

#include <memory>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#endif

namespace Dune {

#if HAVE_DUNE_ISTL
  //! wrapper class to store perconditioner 
  //! as the interface class does not have to category 
  //! enum 
  template<class MatrixImp>
  class PreconditionerWrapper 
    : public Preconditioner<typename MatrixImp :: RowBlockVectorType,
                            typename MatrixImp :: ColBlockVectorType>
  {
    typedef MatrixImp MatrixType;
    typedef typename MatrixImp :: RowBlockVectorType X;
    typedef typename MatrixImp :: ColBlockVectorType Y;
            
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

    //! set preconder to zero 
    PreconditionerWrapper (const PreconditionerWrapper& org) 
      : matrix_(org.matrix_) 
      , preconder_(org.preconder_) 
      , preEx_(org.preEx_)
    {
    }
    
    //! set preconder to zero 
    PreconditionerWrapper (MatrixType& m) 
      : matrix_(m) 
      , preconder_()
      , preEx_(false)  
    {}
    
    //! create preconditioner of given type 
    template <class PreconditionerType>
    PreconditionerWrapper(MatrixType & m,
                            int iter, field_type relax, const PreconditionerType*) 
      : matrix_(m)
      , preconder_(new PreconditionerType(m,iter,relax))
      , preEx_(true) 
    {
    }
    
    //! create preconditioner of given type 
    template <class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, 
                            field_type relax, const PreconditionerType*) 
      : matrix_(m)
      , preconder_(new PreconditionerType(m,relax))
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
  };
#endif // end HAVE_DUNE_ISTL 

} // end namespace Dune 
#endif
