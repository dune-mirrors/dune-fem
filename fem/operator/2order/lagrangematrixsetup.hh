#ifndef DUNE_LAGRANGEMATRIXSETUP_HH
#define DUNE_LAGRANGEMATRIXSETUP_HH

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#endif

namespace Dune
{

  /** \class LagrangeMatrixSetup
   *  \brief Setup Matrix structure for Lagrange operators by including all
   *         Lagrange nodes of an element.
   */
  class LagrangeMatrixSetup
  {
  public:
    //! get number of entries per row for a block matrix, 
    //! i.e. here number of neighboring nodes + 1 
    template< class GridPart >
    static inline int stencilSizeEstimate ( const GridPart &gridPart )
    {
      return 15;//(GridPartImp :: GridType :: dimension * 2) + 1;
    }

    //! create entries for element and neighbors 
    template< class Space, class RowMapper, class ColMapper,
              class MatrixStructureMap, class DiscreteFunction >
    static inline void setup ( const Space &space, 
                               const RowMapper &rowMapper,
                               const ColMapper &colMapper,
                               MatrixStructureMap &indices,
                               const DiscreteFunction*)
    {
      typedef typename Space :: IteratorType IteratorType;

      indices.clear();

      const IteratorType end = space.end();
      for( IteratorType it = space.begin(); it != end; ++it )
        fill( space.gridPart(), *it, rowMapper, colMapper, indices );
    }

  protected:
    //! create entries for element and neighbors 
    template< class GridPart, class Entity,
              class RowMapper, class ColMapper >
    static inline void fill ( const GridPart &gridPart,
                              const Entity &entity,
                              const RowMapper &rowMapper,
                              const ColMapper &colMapper,
                              std :: map< int, std :: set< int > > &indices )
    {
      // type of local indices storage 
      typedef std :: set< int >  LocalIndicesType;

      typedef typename RowMapper :: DofMapIteratorType RowDofMapIteratorType;
      typedef typename ColMapper :: DofMapIteratorType ColDofMapIteratorType;

      const RowDofMapIteratorType rowEnd = rowMapper.end( entity );
      RowDofMapIteratorType rowIt = rowMapper.begin( entity );
      for( ; rowIt != rowEnd; ++rowIt )
      {
        LocalIndicesType &localIndices = indices[ rowIt.global() ];

        const ColDofMapIteratorType colEnd = colMapper.end( entity );
        ColDofMapIteratorType colIt = colMapper.begin( entity );
        for( ; colIt != colEnd; ++colIt )
          localIndices.insert( colIt.global() );
      }
    }
  };

  template <class TraitsImp>
  struct LagrangeMatrixTraits
  {
    typedef typename TraitsImp :: RowSpaceType RowSpaceType;
    typedef typename TraitsImp :: ColumnSpaceType ColumnSpaceType;

    typedef LagrangeMatrixSetup StencilType; 
    
    typedef ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;
  };
  
#if HAVE_DUNE_ISTL
  // forward 
  template <class MatrixImp>
  class LagrangeParallelMatrixAdapter;

  // specialization for ISTL matrices 
  template <class RowSpaceImp, class ColSpaceImp>
  struct LagrangeMatrixTraits<ISTLMatrixTraits<RowSpaceImp,ColSpaceImp> >
  {
    typedef RowSpaceImp RowSpaceType;
    typedef ColSpaceImp ColumnSpaceType;

    typedef LagrangeMatrixSetup StencilType; 
    
    typedef ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;

    template <class MatrixImp>
    struct Adapter
    {   
      // type of matrix adapter 
      typedef LagrangeParallelMatrixAdapter<MatrixImp> MatrixAdapterType;
    };
  };
  
  /*! 
    \brief Adapter to turn a matrix into a linear operator.
    Adapts a matrix to the assembled linear operator interface
  */
  template <class MatrixImp>
  class LagrangeParallelMatrixAdapter
    : public AssembledLinearOperator< MatrixImp,
               typename MatrixImp :: RowBlockVectorType,
               typename MatrixImp :: ColBlockVectorType>
  {
  public:
    typedef MatrixImp MatrixType;
    typedef PreconditionerWrapper<MatrixType> PreconditionAdapterType;
    
    typedef typename MatrixType :: RowDiscreteFunctionType RowDiscreteFunctionType;
    typedef typename MatrixType :: ColDiscreteFunctionType ColumnDiscreteFunctionType;

    typedef typename RowDiscreteFunctionType :: DiscreteFunctionSpaceType RowSpaceType;
    typedef CommunicationManager<RowSpaceType> CommunicationManagerType;

    typedef typename ColumnDiscreteFunctionType :: DiscreteFunctionSpaceType ColSpaceType;
    typedef ParallelScalarProduct<ColumnDiscreteFunctionType> ParallelScalarProductType;
    
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

    mutable CommunicationManagerType comm_;
    mutable ParallelScalarProductType scp_;

    PreconditionAdapterType preconditioner_;
    
  public:  
    //! constructor: just store a reference to a matrix
    LagrangeParallelMatrixAdapter (const LagrangeParallelMatrixAdapter& org)
      : matrix_(org.matrix_) 
      , rowSpace_(org.rowSpace_)
      , colSpace_(org.colSpace_)
      , comm_(rowSpace_)
      , scp_(colSpace_)
      , preconditioner_(org.preconditioner_)
    {}
    //! constructor: just store a reference to a matrix
    LagrangeParallelMatrixAdapter (MatrixType& A,
                             const RowSpaceType& rowSpace, 
                             const ColSpaceType& colSpace) 
      : matrix_(A) 
      , rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , comm_(rowSpace_)
      , scp_(colSpace)
      , preconditioner_(matrix_)
    {}

    //! constructor: just store a reference to a matrix
    template <class PreconditionerType>
    LagrangeParallelMatrixAdapter (MatrixType& A,
                             const RowSpaceType& rowSpace, 
                             const ColSpaceType& colSpace,
                             int iter, field_type relax, const PreconditionerType* dummy) 
      : matrix_(A) 
      , rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , comm_(rowSpace_)
      , scp_(colSpace_)
      , preconditioner_(matrix_,iter,relax,dummy)
    {}

    //! constructor: just store a reference to a matrix
    template <class PreconditionerType>
    LagrangeParallelMatrixAdapter (MatrixType& A,
                             const RowSpaceType& rowSpace, 
                             const ColSpaceType& colSpace, 
                             field_type relax, const PreconditionerType* dummy) 
      : matrix_(A) 
      , rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , comm_(rowSpace_)
      , scp_(colSpace_)
      , preconditioner_(matrix_,relax,dummy)
    {}

    //! return reference to preconditioner 
    PreconditionAdapterType& preconditionAdapter() { return preconditioner_; }

    //! return reference to preconditioner 
    ParallelScalarProductType& scp() { return scp_; }

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
      // apply matrix 
      matrix_.mv(x,y);

      // exchange data first 
      communicate( y );
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
      // apply matrix 
      matrix_.usmv(alpha,x,y);

      // exchange data first 
      communicate( y );
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
    void communicate(Y& y) const 
    {
      if( rowSpace_.grid().comm().size() <= 1 ) return ;
      
      // create temporary discretet function object 
      ColumnDiscreteFunctionType tmp ("LagrangeParallelMatrixAdapter::communicate",
                                   colSpace_, y );

      // exchange data 
      comm_.exchange( tmp , (DFCommunicationOperation :: Add *) 0 );
    }
  };
#endif
} // end namespace Dune 
#endif
