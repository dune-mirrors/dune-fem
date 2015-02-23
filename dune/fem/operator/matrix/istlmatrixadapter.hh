#ifndef DUNE_FEM_ISTLMATRIXADAPTER_HH
#define DUNE_FEM_ISTLMATRIXADAPTER_HH

#if HAVE_DUNE_ISTL

#include <dune/common/timer.hh>

#include <dune/istl/operators.hh>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/padaptivespace/lagrange.hh>
#include <dune/fem/space/padaptivespace/discontinuousgalerkin.hh>
#include <dune/fem/space/combinedspace/combineddiscretefunctionspace.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>

namespace Dune
{
  namespace Fem
  {
    template <class MatrixImp>
    class LagrangeParallelMatrixAdapter
      : public AssembledLinearOperator< MatrixImp,
                 typename MatrixImp :: RowBlockVectorType,
                 typename MatrixImp :: ColBlockVectorType>
    {
      typedef LagrangeParallelMatrixAdapter< MatrixImp > ThisType;
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

      //! define the category
      enum { category=SolverCategory::sequential };

    public:
      //! copy constructur
      LagrangeParallelMatrixAdapter ( const LagrangeParallelMatrixAdapter &org )
      : matrix_( org.matrix_ ),
        rowSpace_( org.rowSpace_ ),
        colSpace_( org.colSpace_ ),
        scp_( colSpace_ ),
        preconditioner_( org.preconditioner_ ),
        averageCommTime_( org.averageCommTime_ )
      {}

      //! constructor: just store a reference to a matrix
      //! and the spaces as well as preconditioning method
      LagrangeParallelMatrixAdapter ( MatrixType &matrix,
                                      const RowSpaceType &rowSpace,
                                      const ColSpaceType &colSpace,
                                      const PreconditionAdapterType& precon )
      : matrix_( matrix ),
        rowSpace_( rowSpace ),
        colSpace_( colSpace ),
        scp_( colSpace_ ),
        preconditioner_( precon ),
        averageCommTime_( 0 )
      {}

      //! return communication time
      double averageCommTime() const
      {
        return averageCommTime_ ;
      }

      //! return reference to preconditioner
      PreconditionAdapterType &preconditionAdapter()
      {
        return preconditioner_;
      }
      //! return reference to preconditioner
      const PreconditionAdapterType &preconditionAdapter() const
      {
        return preconditioner_;
      }

      //! return reference to preconditioner
      ParallelScalarProductType &scp()
      {
        return scp_;
      }

      //! apply operator to x:  \f$ y = A(x) \f$
      virtual void apply ( const X &x, Y &y ) const
      {
        matrix_.mv( x, y );
        communicate( y );
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      virtual void applyscaleadd ( field_type alpha, const X &x, Y &y) const
      {
        if( rowSpace_.grid().comm().size() <= 1 )
        {
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

      virtual double residuum(const Y& rhs, X& x) const
      {
        X tmp (x);

        this->apply(x,tmp);
        tmp -= rhs;

        // return global sum of residuum
        return scp_.norm(tmp);
      }

      //! get matrix via *
      virtual const MatrixType& getmat () const
      {
        return matrix_;
      }

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
    };

    template <class MatrixImp>
    class DGParallelMatrixAdapter
      : public AssembledLinearOperator< MatrixImp,
                 typename MatrixImp :: RowBlockVectorType,
                 typename MatrixImp :: ColBlockVectorType>
    {
      typedef DGParallelMatrixAdapter< MatrixImp > ThisType ;
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

      //! define the category
      enum { category=SolverCategory::sequential };

    protected:
      MatrixType& matrix_;
      const RowSpaceType& rowSpace_;
      const ColSpaceType& colSpace_;

      ParallelScalarProductType scp_;

      PreconditionAdapterType preconditioner_;
      mutable double averageCommTime_;

    public:
      //! constructor: just store a reference to a matrix
      DGParallelMatrixAdapter (const DGParallelMatrixAdapter& org)
        : matrix_(org.matrix_)
        , rowSpace_(org.rowSpace_)
        , colSpace_(org.colSpace_)
        , scp_(colSpace_)
        , preconditioner_(org.preconditioner_)
        , averageCommTime_( org.averageCommTime_ )
      {}

      //! constructor: just store a reference to a matrix
      DGParallelMatrixAdapter (MatrixType& A,
                               const RowSpaceType& rowSpace,
                               const ColSpaceType& colSpace,
                               const PreconditionAdapterType& precon )
        : matrix_(A)
        , rowSpace_(rowSpace)
        , colSpace_(colSpace)
        , scp_(colSpace_)
        , preconditioner_( precon )
        , averageCommTime_( 0.0 )
      {}

      //! return communication time
      double averageCommTime() const
      {
        return averageCommTime_ ;
      }

      //! return reference to preconditioner
      PreconditionAdapterType& preconditionAdapter() { return preconditioner_; }
      const PreconditionAdapterType& preconditionAdapter() const { return preconditioner_; }

      //! return reference to preconditioner
      ParallelScalarProductType& scp() { return scp_; }

      //! apply operator to x:  \f$ y = A(x) \f$
      virtual void apply (const X& x, Y& y) const
      {
        // exchange data first
        communicate( x );

        // apply vector to matrix
        matrix_.mv(x,y);

        // delete non-interior
        scp_.deleteNonInterior( y );
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        // exchange data first
        communicate( x );

        // apply matrix
        matrix_.usmv(alpha,x,y);

        // delete non-interior
        scp_.deleteNonInterior( y );
      }

      virtual double residuum(const Y& rhs, X& x) const
      {
        // exchange data
        communicate( x );

        typedef typename ParallelScalarProductType :: SlaveDofsType SlaveDofsType;
        const SlaveDofsType& slaveDofs = scp_.slaveDofs();

        typedef typename Y :: block_type LittleBlockVectorType;
        LittleBlockVectorType tmp;
        double res = 0.0;

        // counter for rows
        int i = 0;
        const int slaveSize = slaveDofs.size();
        for(int slave = 0; slave<slaveSize; ++slave)
        {
          const int nextSlave = slaveDofs[slave];
          for(; i<nextSlave; ++i)
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

      //! get matrix via *
      virtual const MatrixType& getmat () const
      {
        return matrix_;
      }
    protected:
      void communicate(const X& x) const
      {
        if( rowSpace_.grid().comm().size() <= 1 ) return ;

        Dune::Timer commTime;

        // create temporary discretet function object
        RowDiscreteFunctionType tmp ("DGParallelMatrixAdapter::communicate",
                                     rowSpace_, x );

        // exchange data by copying
        rowSpace_.communicate( tmp );

        // accumulate communication time
        averageCommTime_ += commTime.elapsed();
      }
    };

    template <class MatrixImp,class Space>
    struct ISTLParallelMatrixAdapter;

    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, LagrangeDiscreteFunctionSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef LagrangeParallelMatrixAdapter<MatrixImp> Type;
    };

    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, PAdaptiveLagrangeSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef LagrangeParallelMatrixAdapter<MatrixImp> Type;
    };
    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, DiscontinuousGalerkinSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DGParallelMatrixAdapter<MatrixImp> Type ;
    };
    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, LagrangeDiscontinuousGalerkinSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DGParallelMatrixAdapter<MatrixImp> Type ;
    };
    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, LegendreDiscontinuousGalerkinSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DGParallelMatrixAdapter<MatrixImp> Type ;
    };
    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DGParallelMatrixAdapter<MatrixImp> Type ;
    };
    template< class MatrixImp,
              class FunctionSpace, class GridPart, int polOrder, template< class > class Storage>
    struct ISTLParallelMatrixAdapter< MatrixImp, PAdaptiveDGSpace< FunctionSpace,GridPart,polOrder,Storage> >
    {
      typedef DGParallelMatrixAdapter<MatrixImp> Type ;
    };

    template<class MatrixImp, class Space1, class Space2>
    struct ISTLParallelMatrixAdapter<MatrixImp, CombinedDiscreteFunctionSpace<Space1, Space2> >
    {
      typedef LagrangeParallelMatrixAdapter<MatrixImp> Type;
    };

  } // end namespace Fem
} // end namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLMATRIXADAPTER_HH

