#ifndef DUNE_FEM_ISTLMATRIXADAPTER_HH
#define DUNE_FEM_ISTLMATRIXADAPTER_HH

#if HAVE_DUNE_ISTL

#include <dune/common/timer.hh>
#include <dune/common/version.hh>

#include <dune/istl/operators.hh>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/padaptivespace/lagrange.hh>
#include <dune/fem/space/padaptivespace/discontinuousgalerkin.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{
  namespace Fem
  {

    template <class Matrix>
    class PreconditionerWrapper ;


    template <class MatrixImp>
    class ISTLParallelMatrixAdapterInterface
      : public AssembledLinearOperator< MatrixImp,
                 typename MatrixImp :: RowBlockVectorType,
                 typename MatrixImp :: ColBlockVectorType>
    {
      typedef ISTLParallelMatrixAdapterInterface< MatrixImp > ThisType;
    public:
      typedef MatrixImp MatrixType;
      typedef Fem::PreconditionerWrapper< MatrixType > PreconditionAdapterType;

      typedef typename MatrixType :: RowDiscreteFunctionType RowDiscreteFunctionType;
      typedef typename MatrixType :: ColDiscreteFunctionType ColumnDiscreteFunctionType;

      typedef typename RowDiscreteFunctionType :: DiscreteFunctionSpaceType RowSpaceType;

      typedef typename ColumnDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;
      typedef Fem::ParallelScalarProduct<ColumnDiscreteFunctionType> ParallelScalarProductType;

      typedef typename RowDiscreteFunctionType    :: DofStorageType   X;
      typedef typename ColumnDiscreteFunctionType :: DofStorageType   Y;

      //! export types
      typedef MatrixType  matrix_type;
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

#if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      //! define the category
      enum { category=SolverCategory::sequential };
#endif // #if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)

    public:
      //! copy constructor
      ISTLParallelMatrixAdapterInterface ( const ISTLParallelMatrixAdapterInterface &org )
      : matrix_( org.matrix_ ),
        rowSpace_( org.rowSpace_ ),
        colSpace_( org.colSpace_ ),
        scp_( colSpace_ ),
        preconditioner_( org.preconditioner_ ),
        averageCommTime_( org.averageCommTime_ ),
        threading_( org.threading_ )
      {}

      //! constructor: just store a reference to a matrix
      //! and the spaces as well as preconditioning method
      ISTLParallelMatrixAdapterInterface ( MatrixType &matrix,
                                           const RowSpaceType &rowSpace,
                                           const ColSpaceType &colSpace,
                                           const PreconditionAdapterType& precon,
                                           const bool threading = true )
      : matrix_( matrix ),
        rowSpace_( rowSpace ),
        colSpace_( colSpace ),
        scp_( colSpace_ ),
        preconditioner_( precon ),
        averageCommTime_( 0 ),
        threading_( threading )
      {}

      //! return communication time
      virtual double averageCommTime() const { return averageCommTime_ ; }

      //! return reference to preconditioner
      virtual PreconditionAdapterType &preconditionAdapter() { return preconditioner_; }

      //! return reference to preconditioner
      virtual const PreconditionAdapterType &preconditionAdapter() const { return preconditioner_; }

      //! return reference to preconditioner
      virtual ParallelScalarProductType &scp() { return scp_; }

      //! apply operator to x:  \f$ y = A(x) \f$
      void apply ( const X &x, Y &y ) const override
      {
        DUNE_THROW(NotImplemented,"interface method apply not overloaded!");
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      void applyscaleadd ( field_type alpha, const X &x, Y &y) const override
      {
        DUNE_THROW(NotImplemented,"interface method applyscaleadd not overloaded!");
      }

      //! apply operator and substract right-hand-side
      virtual double residuum(const Y& rhs, X& x) const
      {
        DUNE_THROW(NotImplemented,"interface method residuum not overloaded!");
        return 0.0;
      }

      //! get matrix
      const MatrixType& getmat () const override { return matrix_; }

      //! return true if matvec and preconditioning should use threading (depends on implementation of preconditioning)
      bool threading () const { return threading_; }

    protected:
      void communicate( Y &y ) const
      {
        if( rowSpace_.grid().comm().size() <= 1 )
          return;

        Dune::Timer commTime;
        ColumnDiscreteFunctionType tmp( "LagrangeParallelMatrixAdapter::communicate", colSpace_, y );
        colSpace_.communicate( tmp );
        averageCommTime_ += commTime.elapsed();
      }

      MatrixType &matrix_;
      const RowSpaceType &rowSpace_;
      const ColSpaceType &colSpace_;
      mutable ParallelScalarProductType scp_;
      PreconditionAdapterType preconditioner_;
      mutable double averageCommTime_;
      const bool threading_;
    };

    template <class MatrixImp>
    class LagrangeParallelMatrixAdapter
      : public ISTLParallelMatrixAdapterInterface< MatrixImp >
    {
      typedef LagrangeParallelMatrixAdapter< MatrixImp > ThisType;
      typedef ISTLParallelMatrixAdapterInterface< MatrixImp >     BaseType;
    public:
      typedef MatrixImp MatrixType;
      typedef Fem::PreconditionerWrapper<MatrixType> PreconditionAdapterType;

      typedef typename MatrixType :: RowDiscreteFunctionType RowDiscreteFunctionType;
      typedef typename MatrixType :: ColDiscreteFunctionType ColumnDiscreteFunctionType;

      typedef typename RowDiscreteFunctionType :: DiscreteFunctionSpaceType RowSpaceType;

      typedef typename ColumnDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;
      typedef Fem::ParallelScalarProduct<ColumnDiscreteFunctionType> ParallelScalarProductType;

      typedef typename RowDiscreteFunctionType :: DofStorageType     X;
      typedef typename ColumnDiscreteFunctionType :: DofStorageType  Y;

      //! export types
      typedef MatrixType  matrix_type;
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      using BaseType :: threading;

    protected:
      using BaseType :: matrix_;
      using BaseType :: scp_ ;
      using BaseType :: colSpace_ ;
      using BaseType :: rowSpace_ ;
      using BaseType :: averageCommTime_;

    public:
      //! copy constructor
      LagrangeParallelMatrixAdapter ( const LagrangeParallelMatrixAdapter &org )
      : BaseType( org )
      {}

      //! constructor: just store a reference to a matrix
      //! and the spaces as well as preconditioning method
      LagrangeParallelMatrixAdapter ( MatrixType &matrix,
                                      const RowSpaceType &rowSpace,
                                      const ColSpaceType &colSpace,
                                      const PreconditionAdapterType& precon,
                                      const bool threading = true )
      : BaseType( matrix, rowSpace, colSpace, precon, threading )
      {}

      //! apply operator to x:  \f$ y = A(x) \f$
      void apply ( const X &x, Y &y ) const override
      {
        if( threading() )
          matrix_.mvThreaded(x, y );
        else
          matrix_.mv( x, y );

        communicate( y );
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      void applyscaleadd ( field_type alpha, const X &x, Y &y) const override
      {
        if( rowSpace_.grid().comm().size() <= 1 )
        {
          if( threading() )
            matrix_.usmvThreaded(alpha, x, y );
          else
            matrix_.usmv(alpha,x,y);

          communicate( y );
        }
        else
        {
          Y tmp( y.size() );
          apply( x, tmp );
          y.axpy( alpha, tmp );
        }
      }

      double residuum(const Y& rhs, X& x) const override
      {
        Y tmp( rhs );

        this->apply(x,tmp);
        tmp -= rhs;

        // return global sum of residuum
        return scp_.norm(tmp);
      }

      SolverCategory::Category category () const override { return SolverCategory::sequential; }

    protected:
      void communicate( Y &y ) const
      {
        if( rowSpace_.grid().comm().size() <= 1 )
          return;

        Dune::Timer commTime;
        ColumnDiscreteFunctionType tmp( "LagrangeParallelMatrixAdapter::communicate", colSpace_, y );
        colSpace_.communicate( tmp );
        averageCommTime_ += commTime.elapsed();
      }
    };

    template <class MatrixImp>
    class DGParallelMatrixAdapter
      : public ISTLParallelMatrixAdapterInterface< MatrixImp >
    {
      typedef DGParallelMatrixAdapter< MatrixImp >   ThisType ;
      typedef ISTLParallelMatrixAdapterInterface< MatrixImp > BaseType;
    public:
      typedef MatrixImp MatrixType;
      typedef Fem::PreconditionerWrapper<MatrixType> PreconditionAdapterType;

      typedef typename MatrixType :: RowDiscreteFunctionType RowDiscreteFunctionType;
      typedef typename MatrixType :: ColDiscreteFunctionType ColumnDiscreteFunctionType;

      typedef typename RowDiscreteFunctionType :: DiscreteFunctionSpaceType RowSpaceType;

      typedef typename ColumnDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;
      typedef Fem::ParallelScalarProduct<ColumnDiscreteFunctionType> ParallelScalarProductType;

      typedef typename RowDiscreteFunctionType :: DofStorageType     X;
      typedef typename ColumnDiscreteFunctionType :: DofStorageType  Y;

      //! export types
      typedef MatrixType  matrix_type;
      typedef X domain_type;
      typedef Y range_type;
      typedef typename X::field_type field_type;

      using BaseType :: threading;

    protected:
      using BaseType :: matrix_;
      using BaseType :: scp_;
      using BaseType :: rowSpace_;
      using BaseType :: colSpace_;
      using BaseType :: averageCommTime_;

    public:
      //! constructor: just store a reference to a matrix
      DGParallelMatrixAdapter (const DGParallelMatrixAdapter& org)
        : BaseType( org )
      {}

      //! constructor: just store a reference to a matrix
      DGParallelMatrixAdapter (MatrixType& A,
                               const RowSpaceType& rowSpace,
                               const ColSpaceType& colSpace,
                               const PreconditionAdapterType& precon,
                               const bool threading = true)
        : BaseType( A, rowSpace, colSpace, precon, threading )
      {}

      //! apply operator to x:  \f$ y = A(x) \f$
      void apply (const X& x, Y& y) const override
      {
        // exchange data first
        communicate( x );

        // apply vector to matrix
        if( threading() )
          matrix_.mvThreaded( x, y );
        else
          matrix_.mv(x,y);

        // delete non-interior
        scp_.deleteNonInterior( y );
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      void applyscaleadd (field_type alpha, const X& x, Y& y) const override
      {
        // exchange data first
        communicate( x );

        // apply matrix
        if( threading() )
          matrix_.usmvThreaded(alpha, x, y );
        else
          matrix_.usmv(alpha,x,y);

        // delete non-interior
        scp_.deleteNonInterior( y );
      }

      double residuum(const Y& rhs, X& x) const override
      {
        // exchange data
        communicate( x );

        typedef typename ParallelScalarProductType :: AuxiliaryDofsType AuxiliaryDofsType;
        const AuxiliaryDofsType& auxiliaryDofs = scp_.auxiliaryDofs();

        typedef typename Y :: block_type LittleBlockVectorType;
        LittleBlockVectorType tmp;
        double res = 0.0;

        // counter for rows
        int i = 0;
        const int auxiliarySize = auxiliaryDofs.size();
        for(int auxiliary = 0; auxiliary<auxiliarySize; ++auxiliary)
        {
          const int nextAuxiliary = auxiliaryDofs[auxiliary];
          for(; i<nextAuxiliary; ++i)
          {
            tmp = 0;
            // get row
            typedef typename MatrixType :: row_type row_type;

            const row_type& row = matrix_[i];
            // multiply with row
            typedef typename MatrixType :: ConstColIterator ConstColIterator;
            ConstColIterator endj = row.end();
            for (ConstColIterator j = row.begin(); j!=endj; ++j)
            {
              (*j).umv(x[j.index()], tmp);
            }

            // substract right hand side
            tmp -= rhs[i];

            // add scalar product
            res += tmp.two_norm2();
          }
          ++i;
        }

        res = rowSpace_.grid().comm().sum( res );
        // return global sum of residuum
        return std::sqrt( res );
      }

      SolverCategory::Category category () const override { return SolverCategory::sequential; }

    protected:
      void communicate(const X& x) const
      {
        if( rowSpace_.grid().comm().size() <= 1 ) return ;

        Dune::Timer commTime;

        // create temporary discrete function object
        RowDiscreteFunctionType tmp ("DGParallelMatrixAdapter::communicate",
                                     rowSpace_, x );

        // exchange data by copying
        rowSpace_.communicate( tmp );

        // accumulate communication time
        averageCommTime_ += commTime.elapsed();
      }
    };

  } // end namespace Fem
} // end namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLMATRIXADAPTER_HH

