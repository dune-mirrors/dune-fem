#ifndef DUNE_LAGRANGEMATRIXSETUP_HH
#define DUNE_LAGRANGEMATRIXSETUP_HH

#include <map>
#include <set>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/misc/functor.hh>

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
  template <bool addNonConformingNeighbors >
  class LagrangeMatrixSetup
  {
  public:
    //! get number of entries per row for a block matrix, 
    template< class Space >
    static inline int nonZerosEstimate( const Space &space )
    {
      enum { dimension = Space :: GridType :: dimension };
      return 8 * (1 << dimension) * space.order();
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
    template <class GK>
    struct FillFunctor
    {
      typedef std :: set< int >  LocalIndicesType;
      typedef GK GlobalKey;
      FillFunctor(std::map<int,LocalIndicesType> &indices) 
      : indices_(indices),
        localIndices_(0)
      {}
      void set(const std::size_t rowLocal, const GlobalKey &rowGlobal)
      {
        localIndices_ = &(indices_[ rowGlobal ]);
      }
      template <class ColGlobal>
      void operator() ( const std::size_t colLocal, const ColGlobal &colGlobal)
      {
        localIndices_->insert( colGlobal );
      }
      private:
      std::map<int,LocalIndicesType> &indices_;
      LocalIndicesType *localIndices_;
    };
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
      typedef FillFunctor<typename RowMapper::GlobalKeyType> Functor;
      typedef Fem::MatrixFunctor<ColMapper,Entity,Functor > MFunctor;
      rowMapper.mapEach(entity, MFunctor( colMapper, entity, Functor(indices) ) );
 
      // in case of non-conforming stabilization 
      if( addNonConformingNeighbors )
      {
        typedef typename GridPart :: IntersectionIteratorType IntersectionIteratorType;
        typedef typename IntersectionIteratorType :: Intersection IntersectionType;
        typedef typename GridPart :: template Codim<0> :: EntityPointerType EntityPointer;

        const IntersectionIteratorType endit = gridPart.iend( entity );
        for( IntersectionIteratorType it = gridPart.ibegin( entity );
             it != endit; ++it )
        {
          const IntersectionType& intersection = *it;
          if( intersection.neighbor() && ! intersection.conforming() ) 
          {
            EntityPointer ep = intersection.outside(); 
            const Entity& neighbor = *ep;
            rowMapper.mapEach(entity, MFunctor( colMapper, neighbor, Functor(indices) ) );
          }
        }
      }

#if 0
      // debug output 
      for(RowDofMapIteratorType rowIt = rowMapper.begin( entity );
          rowIt != rowEnd; ++rowIt )
      {
        LocalIndicesType &localIndices = indices[ rowIt.global() ];
        std::cout << "row[ " << rowIt.global() << "] contains: " <<std::endl;
        typedef typename LocalIndicesType :: iterator iterator;
        const iterator endlocal = localIndices.end();
        for(iterator it = localIndices.begin(); it != endlocal; ++it)
        {
          std::cout << *it << " ";
        }
        std::cout << std::endl;
      }
#endif
    }
  };

  /** \class LagrangeMatrixSetup
   *  \brief Setup Matrix structure for Lagrange operators by including all
   *         Lagrange nodes of an element for conforming Lagrange Elements.
   */
  struct LagrangeConformingMatrixSetup : public LagrangeMatrixSetup< false >
  {
  };


  /** \class LagrangeMatrixSetup
   *  \brief Setup Matrix structure for Lagrange operators by including all
   *         Lagrange nodes of an element for non-conforming Lagrange Elements.
   */
  struct LagrangeNonConformingMatrixSetup : public LagrangeMatrixSetup< true >
  {
  };

  /** \brief LagrangeMatrixTraits traits class for Matrix Objects 
      used to setup system matrices for FE operators 
      \note For standard Lagrange conforming Finite Elements the 
      addNonConformingNeighbors flag is false.
      For others is flag could be equal the inverse to GridPart :: conforming 
  */
  template <class TraitsImp, bool addNonConformingNeighbors = false>
  struct LagrangeMatrixTraits
  {
    typedef typename TraitsImp :: RowSpaceType RowSpaceType;
    typedef typename TraitsImp :: ColumnSpaceType ColumnSpaceType;

    typedef LagrangeMatrixSetup<addNonConformingNeighbors> StencilType; 
    
    typedef Fem::ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;
  };
  
#if HAVE_DUNE_ISTL
  // forward 
  template <class MatrixImp>
  class LagrangeParallelMatrixAdapter;

  // specialization for ISTL matrices 
  template <class RowSpaceImp, class ColSpaceImp, bool addNonConformingNeighbors>
  struct LagrangeMatrixTraits< Fem::ISTLMatrixTraits<RowSpaceImp,ColSpaceImp>, addNonConformingNeighbors >
  {
    typedef RowSpaceImp RowSpaceType;
    typedef ColSpaceImp ColumnSpaceType;

    typedef LagrangeMatrixSetup<addNonConformingNeighbors> StencilType; 
    
    typedef Fem::ParallelScalarProduct < ColumnSpaceType > ParallelScalarProductType;

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

      Timer commTime; 
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

#endif // #if HAVE_DUNE_ISTL

} // end namespace Dune 

#endif
