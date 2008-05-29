#ifndef DUNE_DGMATRIXSETUP_HH
#define DUNE_DGMATRIXSETUP_HH

#include <dune/fem/space/common/gridpartutility.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/space/common/commoperations.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/operators.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#endif

namespace Dune {

////////////////////////////////////////////////////////////
//
//  Setup of matrix structure 
//
////////////////////////////////////////////////////////////
/**  \brief Setup Matrix structure for DG operators by including
 * elements and it's neighbors. 
*/
class ElementAndNeighbors
{
public:
  //! get number of entries per row for a block matrix, 
  //! i.e. here number of neighbors + 1
  template <class GridPartImp> 
  static inline int stencilSizeEstimate(const GridPartImp& gridPart) 
  {
    return (GridPartImp :: GridType :: dimension * 2) + 1; 
  }

  //! create entries for element and neighbors 
  template <class SpaceImp,    
            class RowMapperType,
            class ColMapperType,
            class MatrixStructureMapImp,
            class DiscreteFunctionType>
  static inline void setup(const SpaceImp& space,    
                           const RowMapperType& rowMapper,
                           const ColMapperType& colMapper,
                           MatrixStructureMapImp& indices,
                           const DiscreteFunctionType* )
  {
    typedef typename SpaceImp :: GridPartType GridPartImp;
    GridPartImp& gridP = space.gridPart();

    typedef typename GridPartNewPartitionType<
      GridPartImp,All_Partition> :: NewGridPartType GridPartType;    

    const GridPartType gridPart ( gridP.grid() );
    
    //typedef typename SpaceImp :: GridPartType GridPartType;
    //const GridPartType& gridPart = space.gridPart();

    typedef ParallelScalarProduct<DiscreteFunctionType> ParallelScalarProductType;
    typedef typename ParallelScalarProductType :: BuildProxyType BuildProxyType;
    
    ParallelScalarProductType scp (space);

    std::auto_ptr<BuildProxyType> buildProxy = scp.buildProxy();

    // define used types 
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridPartType :: template Codim<0> :: IteratorType  IteratorType;

    // clear map 
    indices.clear();

    // we need All_Partition here to insert overlap entities 
    // only for diagonal 
    IteratorType endit = gridPart.template end<0>(); 
    for(IteratorType it = gridPart.template begin<0>(); it != endit; ++it)
    {
      const EntityType & en = *it;
      // add all column entities to row  
      fill(gridPart,en,rowMapper,colMapper,indices, *buildProxy);
    }

    // insert size as last ghost 
    buildProxy->insert( rowMapper.size() );
  }

protected:
  //! create entries for element and neighbors 
  template <class GridPartImp,
            class EntityImp,
            class RowMapperImp,
            class ColMapperImp,
            class ParallelScalarProductType>
  static inline void fill(const GridPartImp& gridPart,
                   const EntityImp& en,
                   const RowMapperImp& rowMapper,
                   const ColMapperImp& colMapper,
                   std::map< int , std::set<int> >& indices,
                   ParallelScalarProductType& slaveDofs)
  {
    assert( rowMapper.maxNumDofs () == 1 );
    // get index for entity 
    const int elRowIndex = rowMapper.mapToGlobal( en, 0 ); 

    // type of local indices storage 
    typedef std::set< int >  LocalIndicesType; 
    LocalIndicesType& localIndices = indices[elRowIndex];

    // insert diagonal for each element 
    localIndices.insert( elRowIndex );

    // if entity is not interior, insert into overlap entries 
    if(en.partitionType() != InteriorEntity)
    {
      slaveDofs.insert( elRowIndex );
    }

    // insert neighbors 
    typedef typename GridPartImp:: GridType :: template Codim<0>::EntityPointer EntityPointerType; 
    typedef typename GridPartImp:: IntersectionIteratorType IntersectionIteratorType;
    IntersectionIteratorType endnit = gridPart.iend(en);
    for(IntersectionIteratorType nit = gridPart.ibegin(en);
        nit != endnit; ++nit)
    {
      if(nit.neighbor())
      {
        // get neighbor 
        EntityPointerType ep = nit.outside();
        const EntityImp& nb = *ep;

        // get index of neighbor 
        const int nbColIndex = colMapper.mapToGlobal( nb , 0 );
        const int nbRowIndex = rowMapper.mapToGlobal( nb , 0 );

        // check whether to insert now 
        bool insertHere = (elRowIndex < nbRowIndex);
        bool nbInsert = true;
#if HAVE_MPI 
        // check partition type 
        if( nb.partitionType() != InteriorEntity )
        {
          insertHere = true;
          nbInsert = nb.partitionType() != GhostEntity;
          slaveDofs.insert( nbRowIndex );
        }
#endif
        // insert pair 
        if( insertHere )
        {
          // insert neighbor 
          localIndices.insert( nbColIndex );

          // insert symetric part with swaped row-col
          LocalIndicesType& nbIndices = indices[nbRowIndex];
          nbIndices.insert( nbColIndex );

          if( nbInsert )
          {
            const int elColIndex = colMapper.mapToGlobal( en , 0 );
            nbIndices.insert( elColIndex );  
          }
        }
      }
    }
  }
};

  template <class TraitsImp>
  struct DGMatrixTraits
  {
    typedef typename TraitsImp :: RowSpaceType RowSpaceType;
    typedef typename TraitsImp :: ColumnSpaceType ColumnSpaceType;

    typedef ElementAndNeighbors StencilType; 
    
    typedef ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;
  };
  
#if HAVE_DUNE_ISTL
  // forward 
  template <class MatrixImp>
  class DGParallelMatrixAdapter;

  // specialization for ISTL matrices 
  template <class RowSpaceImp, class ColSpaceImp>
  struct DGMatrixTraits<ISTLMatrixTraits<RowSpaceImp,ColSpaceImp> >
  {
    typedef RowSpaceImp RowSpaceType;
    typedef ColSpaceImp ColumnSpaceType;

    typedef ElementAndNeighbors StencilType; 
    
    typedef ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;

    template <class MatrixImp>
    struct Adapter
    {   
      // type of matrix adapter 
      typedef DGParallelMatrixAdapter<MatrixImp> MatrixAdapterType;
    };
  };
  
  /*! 
    \brief Adapter to turn a matrix into a linear operator.
    Adapts a matrix to the assembled linear operator interface
  */
  template <class MatrixImp>
  class DGParallelMatrixAdapter
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
    ParallelScalarProductType scp_;

    PreconditionAdapterType preconditioner_;
    
  public:  
    //! constructor: just store a reference to a matrix
    DGParallelMatrixAdapter (const DGParallelMatrixAdapter& org)
      : matrix_(org.matrix_) 
      , rowSpace_(org.rowSpace_)
      , colSpace_(org.colSpace_)
      , comm_(rowSpace_)
      , scp_(colSpace_)
      , preconditioner_(org.preconditioner_)
    {}
    //! constructor: just store a reference to a matrix
    DGParallelMatrixAdapter (MatrixType& A,
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
    DGParallelMatrixAdapter (MatrixType& A,
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
    DGParallelMatrixAdapter (MatrixType& A,
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
      
      // create temporary discretet function object 
      RowDiscreteFunctionType tmp ("DGParallelMatrixAdapter::communicate",
                                   rowSpace_, x );

      // exchange data by copying 
      comm_.exchange( tmp, (DFCommunicationOperation :: Copy*) 0);
    }
  };
#endif

} // end namespace Dune 
#endif
