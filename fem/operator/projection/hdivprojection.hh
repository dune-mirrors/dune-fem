#ifndef DUNE_HDIV_PROJECTION_HH 
#define DUNE_HDIV_PROJECTION_HH 

//- Dune includes
#include <dune/grid/common/referenceelements.hh>

//- Dune-fem includes 
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh> 
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>

// make sure higher order Lagrange works (define USE_TWISTFREE_MAPPER)
#include <dune/fem/space/lagrangespace.hh>

#ifdef ENABLE_UG 
#include <dune/grid/uggrid.hh>
#endif

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class YaspGrid;



/** 
    \ingroup HdivProjectionOperator 
    \class HdivProjection 
    \brief H-div Projection for discontinuous discrete functions.
    The projection is described in detail in: 

    P. Bastian and B. Riviere. Superconvergence and H(div)-projection
    for discontinuous Galerkin methods. Int. J. Numer. Meth. Fluids.,
    42:1043-1057, 2003.

    (see homepage of Peter Bastian: 
      http://hal.iwr.uni-heidelberg.de/~peter/Papers/BDMpaper.pdf )

    Note:  

    This projection only works for polynomial order 1 and the following
    spaces: 
      - SimplexGrids + DiscontinuousGalerkinSpace
      - CubeGrids    + LegendreDiscontinuousGalerkinSpace 
    
*/
template <class DiscreteFunctionType>
class HdivProjection : public SpaceOperatorInterface<DiscreteFunctionType>
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType; 
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType; 
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType; 
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType; 
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType; 
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType; 
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType; 

  typedef typename GridPartType :: GridType GridType;
  
    enum { dimRange = 1 };
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain - 1 };
    enum { polOrdN = DiscreteFunctionSpaceType::polynomialOrder };

  typedef FunctionSpace<DomainFieldType,RangeFieldType,dimDomain,dimRange> FaceSpaceType;
  typedef FunctionSpace<DomainFieldType,RangeFieldType,DiscreteFunctionSpaceType::dimDomain,dimRange> ElSpaceType;

  enum { gradPolOrd = ((polOrdN - 1) < 0) ? 0 : (polOrdN - 1) };
  enum { bubblePolOrd = (polOrdN + 1) };
  
  template <class Space> struct Spaces; 
  
  template <class FunctionSpaceImp,
            class GridPartImp,
            int polOrd,
            template <class> class StorageImp,
            template <class,class,int,template <class> class> class DiscreteFunctionSpaceImp>
  struct Spaces<DiscreteFunctionSpaceImp<FunctionSpaceImp,GridPartImp,polOrd,StorageImp> >
  {
    //typedef LagrangeDiscreteFunctionSpace<ElSpaceType,GridPartImp,gradPolOrd,StorageImp> ElementGradientSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace<ElSpaceType,GridPartImp,gradPolOrd,StorageImp> ElementGradientSpaceType;
    typedef DiscreteFunctionSpaceImp<ElSpaceType,GridPartImp,gradPolOrd,StorageImp> ElementGradientSpaceType;
    //typedef DiscontinuousGalerkinSpace<ElSpaceType,GridPartImp,gradPolOrd,StorageImp> ElementGradientSpaceType;
    typedef DiscreteFunctionSpaceImp<FaceSpaceType,GridPartType,polOrdN,StorageImp> FaceDiscreteSpaceType; 
    //typedef LegendreDiscontinuousGalerkinSpace<FaceSpaceType,GridPartType,polOrdN,StorageImp> FaceDiscreteSpaceType; 
    //typedef DiscontinuousGalerkinSpace<FaceSpaceType,GridPartType,polOrdN,StorageImp> FaceDiscreteSpaceType; 
    typedef LagrangeDiscreteFunctionSpace< ElSpaceType, GridPartType,  bubblePolOrd,StorageImp> ElementDiscreteSpaceType; 
  };
  
  template <class DiscreteFunctionSpaceImp,
            int N, 
            Dune::DofStoragePolicy policy> 
  struct Spaces<CombinedSpace<DiscreteFunctionSpaceImp,N,policy> >
  {
    typedef typename DiscreteFunctionSpaceImp :: GridPartType GridPartImp;
    typedef Spaces<DiscreteFunctionSpaceImp> AllSpacesType; 
    typedef typename AllSpacesType :: ElementGradientSpaceType ElementGradientSpaceType;
    typedef typename AllSpacesType :: ElementDiscreteSpaceType ElementDiscreteSpaceType;
    typedef typename AllSpacesType :: FaceDiscreteSpaceType    FaceDiscreteSpaceType;
  };
  
  typedef Spaces<DiscreteFunctionSpaceType> AllSpacesType;
  typedef typename AllSpacesType :: ElementGradientSpaceType ElementGradientSpaceType;
  typedef typename AllSpacesType :: ElementDiscreteSpaceType ElementDiscreteSpaceType;
  typedef typename AllSpacesType :: FaceDiscreteSpaceType    FaceDiscreteSpaceType;

  const DiscreteFunctionSpaceType& space_;
  GridPartType & gridPart_;
  const FaceDiscreteSpaceType faceSpace_;
  const ElementDiscreteSpaceType elSpace_;
  const ElementGradientSpaceType gradSpace_;

 public:
  //! constructor taking space 
  HdivProjection(const DiscreteFunctionSpaceType& space) : 
    space_(space),
    gridPart_(space.gridPart()),
    faceSpace_( gridPart_ ),
    elSpace_( gridPart_ ),
    gradSpace_( gridPart_ )
  {}

  HdivProjection(const HdivProjection& org) : 
    space_(org.space_),
    gridPart_( org.gridPart_),
    faceSpace_( gridPart_ ),
    elSpace_( gridPart_ ),
    gradSpace_( gridPart_ )
  {}

  //! return reference to space 
  virtual const DiscreteFunctionSpaceType& space() const 
  {
    return space_; 
  }
  virtual void setTime(double) {
  }
  virtual double timeStepEstimate() const {
    return 0.;
  }

  //! return sum of jumps of discrete function normal to intersection 
  static double normalJump(const DiscreteFunctionType &discFunc, const int polyOrder = -1 ) 
  {
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType :: Intersection IntersectionType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 
   
    enum { dim = GridType::dimension };
    
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    const GridPartType & gridPart = space.gridPart();
    const int polOrd = (polyOrder <0) ? (2 * space.order() + 2) : polyOrder;
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    RangeType ret (0.0);
    RangeType neighRet (0.0);
    
    // define type of quadrature 
    typedef ElementQuadrature <GridPartType , 1> FaceQuadratureType; 

    double sum = 0.0;

    const LocalIdSetType & idSet = space.grid().localIdSet();
    
    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      EntityType & en = *it;
      const LocalFuncType lf = discFunc.localFunction(en);
     
      double localValue = 0.0;
      
      IntersectionIteratorType endnit = gridPart.iend(en);
      for(IntersectionIteratorType nit = gridPart.ibegin(en);
          nit != endnit; ++nit)
      {
        const IntersectionType& inter=*nit;
        // only interior faces are considered 
        if(inter.neighbor() )
        {
          EntityPointerType neighEp = inter.outside();
          EntityType&            nb = *neighEp;
          
          if(idSet.id( en ) < idSet.id( nb ))
          {
            FaceQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
            FaceQuadratureType faceQuadOuter(gridPart, inter, polOrd, FaceQuadratureType::OUTSIDE);

            const LocalFuncType neighLf = discFunc.localFunction(nb);
            
            const int quadNop = faceQuadInner.nop();
            for (int l = 0; l < quadNop ; ++l)
            {
              DomainType normal = 
                inter.unitOuterNormal(faceQuadInner.localPoint(l));

              lf.evaluate(faceQuadInner[l], ret);
              neighLf.evaluate(faceQuadOuter[l], neighRet);
              
              ret -= neighRet;
              
              double val = ret * normal; 
              val *= faceQuadInner.weight(l);

              localValue += std::abs(val);
            }
          }
        }
      } // end of intersection iterator 
      sum += std::abs(localValue);
    }

    return sum;
  }
  
private:  

  // only works for 2d right now 
  void curl(const DomainType & arg, DomainType & dest, const int d ) const 
  {
    if( DomainType :: dimension == 2 ) 
    {
      dest[0] =  arg[1];
      dest[1] = -arg[0];

      return ;
    }
    else if( DomainType :: dimension == 3 )
    {
      if( d == 0 ) 
      {
        dest[0] =  arg[1];
        dest[1] = -arg[0];
        dest[2] =  0; 
        return ;
      }
      else 
      {
        dest[0] =  arg[2];
        dest[1] =  0;
        dest[2] = -arg[0]; 
        return ;
      }
    }
    else 
    {
      assert( false );
      abort();
    }
  }

  template <class EntityType, 
            class LocalFunctionType, 
            class ArrayType, 
            class MatrixType, 
            class VectorType> 
  void bubblePart(const ElementDiscreteSpaceType& space, 
                  EntityType & en, 
                  const LocalFunctionType & uLF, const int startRow , 
                  ArrayType& uRets, 
                  MatrixType & matrix, VectorType& rhs ) const 
  {

    typedef typename ElementDiscreteSpaceType :: BaseFunctionSetType BaseFunctionSetType;
    typedef typename ElementDiscreteSpaceType :: LagrangePointSetType  LagrangePointSetType;
    typedef typename ElementDiscreteSpaceType :: Traits::GridPartType GridPartType;

    typedef typename ElementDiscreteSpaceType :: JacobianRangeType JacobianRangeType; 
    typedef typename ElementDiscreteSpaceType :: DomainType DomainType;
    
    enum { dim = GridType::dimension };
    
    const LagrangePointSetType& lagrangePointSet = space.lagrangePointSet( en ); 
    const BaseFunctionSetType bSet = space.baseFunctionSet( en ); 
    const int polOrd = 2 * space.order(); // + 2;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    const int cols = uLF.numDofs();
    assert( uRets.size() == cols );

    QuadratureType quad (en,polOrd);
    DomainType result;
    JacobianRangeType valTmp;
    JacobianRangeType val;
    DomainType bVal; 
    DomainType aVal; 

    // get geometry 
    typedef typename EntityType :: Geometry Geometry ;
    const Geometry& geo = en.geometry();
    // get geometry type 
    const GeometryType& type = geo.type(); 
    const int bubbleOffset = (type.isSimplex()) ? 0 : baseFunctionOffset( 0 );

    // get number of dofs for codim 0 (skip first for) 
    const int enDofs = numberOfBubbles( lagrangePointSet.numDofs( 0 ), type );

    // type of jacobian inverse 
    typedef typename Geometry :: ctype ctype; 
    enum { cdim  = Geometry :: coorddimension };
    enum { mydim = Geometry :: mydimension    };
    typedef FieldMatrix<ctype, cdim, mydim> JacobianInverseType;

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      // get jacobian inverse 
      const JacobianInverseType& inv = geo.jacobianInverseTransposed( quad.point(l) );
      
      // get integration element 
      const double intel = quad.weight(l) * 
        geo.integrationElement(quad.point(l));

      // evaluate u 
      uLF.evaluate(quad[l], result);

      // evaluate base functions of u 
      for(int j=0; j<cols; ++j)     
      {
        uLF.baseFunctionSet().evaluate(j, quad[l], uRets[j] );
      }

      // for all bubble functions 
      for( int i = 0 ; i<enDofs; ++i) 
      {
        // we might have other row 
        int row = startRow + i;

        // map to lagrange base function number 
        const int baseFct = 
          lagrangePointSet.entityDofNumber( 0, 0, i + bubbleOffset ); 
        
        // evaluate gradient 
        bSet.jacobian( baseFct, quad[l], valTmp );

        //apply inverse 
        inv.mv( valTmp[0], val[0] );

        for(int d = 0; d<mydim-1; ++d ) 
        {
          // apply curl 
          curl(val[0], aVal, d);

          double r = aVal * result; 
          r *= intel; 
          rhs[row] += r; 

          // for cols make matrix 
          for(int j=0; j<cols; ++j)     
          {
            double t = aVal * uRets[ j ];
            t *= intel; 
            matrix[ row ][ j ] += t;
          }

          // increase row 
          row += d;
        }
      }
    }
  }

  template <class EntityType, class LocalFunctionType, class ArrayType, 
            class MatrixType, class VectorType> 
  void gradientPart(const ElementGradientSpaceType & space, 
                    EntityType & en, 
                    const LocalFunctionType & uLF, const int startRow ,
                    ArrayType& uRets, MatrixType& matrix, VectorType& rhs ) const 
  {
    typedef typename ElementGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
    const BaseFunctionSetType bSet = space.baseFunctionSet( en ); 
    int polOrd = 2 * space.order() + 1;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    const int localRows = gradientBaseFct( bSet ); 
    const int cols = uLF.numDofs();

    QuadratureType quad (en,polOrd);

    RangeType result;
    RangeType uPhi;

    typedef typename ElementGradientSpaceType :: JacobianRangeType GradJacobianRangeType; 
    GradJacobianRangeType gradTmp;
    GradJacobianRangeType gradPhi;

    typedef typename EntityType :: Geometry Geometry ;
    const Geometry& geo = en.geometry();

    typedef typename Geometry :: ctype ctype; 
    enum { cdim  = Geometry :: coorddimension };
    enum { mydim = Geometry :: mydimension    };
    typedef FieldMatrix<ctype, cdim, mydim> JacobianInverseType;

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      // get jacobian inverse 
      const JacobianInverseType& inv = geo.jacobianInverseTransposed( quad.point( l ) );

      // get integration element 
      const double intel = quad.weight(l) * 
        geo.integrationElement( quad.point( l ) );
      
      // evaluate uLF 
      uLF.evaluate(quad[l], result);

      // evaluate base function on quadrature point 
      for(int j=0; j<cols; ++j)     
      {
        uLF.baseFunctionSet().evaluate(j, quad[l], uRets[ j ] );
      }

      for(int i=0; i<localRows; ++i)     
      {
        // we might have other row 
        const int row = startRow + i;

        // evaluate gradient (skip first function because this function
        // is constant and the gradient therefore 0 )
        bSet.jacobian( baseFunctionOffset( i ), quad[l], gradTmp);

        // apply jacobian Inverse 
        inv.mv( gradTmp[0], gradPhi[0] );

        const double uDGVal = result * gradPhi[0];
        rhs[row] += uDGVal * intel;  

        // for cols make matrix 
        for(int j=0; j<cols; ++j)     
        {
          //uLF.baseFunctionSet().evaluate(j, quad[l], uPhi);
          //const double val = uPhi * gradPhi[0]; 
          const double val = uRets[ j ] * gradPhi[0]; 
          matrix[row][j] += val * intel;
        }
      }
    }
  }

  template <class FaceBSetType, class GridType> 
  struct GetSubBaseFunctionSet
  {
    template <class EntityType, class SpaceType> 
    static inline FaceBSetType faceBaseSet(const EntityType& en, const SpaceType& space) 
    {
      return space.baseFunctionSet( *(en.template entity<1> (0) ));
    }
  };

  template <class FaceBSetType, int dim> 
  struct GetSubBaseFunctionSet<FaceBSetType, YaspGrid<dim,dim> >
  {
    template <class EntityType, class SpaceType> 
    static inline FaceBSetType faceBaseSet(const EntityType& en, const SpaceType& space) 
    {
      const GeometryType geoType (GeometryType::cube,dim-1);
      return space.baseFunctionSet( geoType ); 
    }
  };

#ifdef ENABLE_UG 
  template <class FaceBSetType, int dim> 
  struct GetSubBaseFunctionSet<FaceBSetType, UGGrid<dim> >
  {
    template <class EntityType, class SpaceType> 
    static inline FaceBSetType faceBaseSet(const EntityType& en, const SpaceType& space) 
    {
      const GeometryType geoType (en.geometry().type().basicType(),dim-1);
      return space.baseFunctionSet( geoType ); 
    }
  };
#endif

  enum { gradFuncOffset = 1 };
  template <class GradBaseFunctionSet>
  int gradientBaseFct(const GradBaseFunctionSet& gradSet) const 
  {
    return (gradPolOrd <= 0) ? 0 : gradSet.numBaseFunctions() - gradFuncOffset;
  }

  int baseFunctionOffset(const int i) const 
  {
    return i + gradFuncOffset; 
  }

  int numberOfBubbles( const int bubbles , const GeometryType& type) const 
  {
    return (type.isSimplex()) ? bubbles : bubbles - gradFuncOffset;
  }

  //! do projection of discrete functions  
  void project(const DiscreteFunctionType &uDG,
               DiscreteFunctionType & velo ) const 
  {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    enum { localBlockSize = FunctionSpaceType::localBlockSize };

    typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType; 
    typedef typename FunctionSpaceType :: DomainType DomainType;
    
    enum { dim = GridType::dimension };
    typedef typename GridType :: ctype coordType;
    
    const FunctionSpaceType& space = uDG.space();

    // for polOrd 0 this is not working 
    if(space.order() < 1 ) return ;

    const int polOrd = 2 * space.order() + 2;

    typedef typename FaceDiscreteSpaceType :: BaseFunctionSetType FaceBSetType  ; 
    typedef typename FaceDiscreteSpaceType :: RangeType FaceRangeType; 
    
    typedef typename ElementDiscreteSpaceType :: BaseFunctionSetType ElementBaseSetType  ; 
    typedef typename ElementGradientSpaceType :: BaseFunctionSetType GradientBaseSetType  ; 
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    Iterator start = space.begin();
    // if empty grid, do nothing 
    if( start == space.end() ) return; 

    // only working for spaces with one element type 
    if( space.multipleGeometryTypes() )
    {
      if( space.indexSet().geomTypes(0).size() > 1)
      {
        DUNE_THROW(NotImplemented,"H-div projection not implemented for hybrid grids"); 
      }
    }

    // colums are dofs of searched function 
    LocalFuncType lf = uDG.localFunction(*start); 
    const int numDofs = lf.numDofs();
    //std::cout << numDofs << " numDofs \n";

    const FaceBSetType faceSet = 
      GetSubBaseFunctionSet<FaceBSetType,GridType>::faceBaseSet( *start , faceSpace_ ); 
    // number of dofs on faces 
    const int numFaceDofs = faceSet.numBaseFunctions();
    
    const GeometryType startType = start->type();

    // get all element dofs from Lagrange space 
    const int numBubbleDofs = (space.order() <= 1) ? 0 : 
          numberOfBubbles( elSpace_.lagrangePointSet( *start ).numDofs ( 0 ) , startType );

    //std::cout << numBubbleDofs << " bubbleDofs \n";
  
    const GradientBaseSetType gradSet = gradSpace_.baseFunctionSet(*start);
    // in case of linear space the is zero 
    const int numGradDofs = gradientBaseFct( gradSet ); 
    //std::cout << numGradDofs << " numGradDofs \n";
  
    const ReferenceElement< coordType, dim > & refElem =
        ReferenceElements< coordType, dim >::general( startType );

    // get number of faces 
    const int overallFaceDofs = numFaceDofs * refElem.size(1);
    //std::cout << overallFaceDofs << " allFAceDofs \n";

    const int rows = (overallFaceDofs + numGradDofs + numBubbleDofs);
    
    // number of columns 
    const int cols = numDofs; 
    //std::cout << "cols " << cols << " | rows " << rows << "\n";

    // check rows == cols 
    if( space.order() == 1 )
    {
      if( (cols != rows) && dim > 2 ) 
        DUNE_THROW(InvalidStateException,"H-div for order 1 only works with symetric matrices in 3d"); 
    }

    MutableArray< RangeFieldType > rets(numDofs);
    MutableArray< RangeType > uRets(numDofs);

    typedef FieldMatrix<RangeFieldType,localBlockSize,localBlockSize> FieldMatrixType;

    // resulting system matrix 
    FieldMatrixType inv;
      
    typedef FieldVector<RangeFieldType,localBlockSize> VectorType; 
    VectorType fRhs(0.0);
    VectorType x(0.0);

    assert( numDofs == localBlockSize );
    if( numDofs != localBlockSize ) 
    {
      DUNE_THROW(InvalidStateException,"wrong sizes ");
    }

    // flag to say whether we have non-symetric or symetric situation 
    const bool nonSymetric = (cols != rows);

    // matrix type 
    typedef DenseMatrix<double> MatrixType; 
    MatrixType matrix(rows,cols);
    typedef typename MatrixType :: RowType RowType;
    // matrix type 
    MatrixType fakeMatrix(cols,cols);
    RowType rhs(rows,0.0);
    RowType fakeRhs(numDofs,0.0);

    // iterate over grid 
    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // get entity 
      const EntityType& en = *it;

      if( nonSymetric ) 
      {
        // reset values 
        matrix = 0.0;

        // reset rhs 
        for(int i=0; i<rows; ++i) 
        {
          rhs[i] = 0.0;
        }

        // fill non-symetric matrix 
        fillMatrix(gridPart_,en,uDG,
                   faceSpace_,
                   gradSpace_, overallFaceDofs, 
                   elSpace_, rows - numBubbleDofs, 
                   polOrd,numDofs,numFaceDofs,
                   rets,uRets, matrix,rhs);

        // apply least square 
        matrix.multTransposed(rhs, fakeRhs);
        fakeMatrix.multiply_AT_A(matrix);

        // copy values 
        for(int i=0; i<numDofs; ++i)
        {
          fRhs[i] = fakeRhs[i];
          for(int j=0; j<numDofs; ++j) 
          {
            inv[i][j] = fakeMatrix[i][j];
          }
        }
      }
      else 
      {
        // reset values 
        inv = 0.0;
        // reset rhs 
        fRhs = 0.0;

        assert( cols == rows );
        // fill inv and fRhs directly 
        fillMatrix(gridPart_,en,uDG,
                   faceSpace_,
                   gradSpace_, rows - numBubbleDofs - numGradDofs,
                   elSpace_, rows - numBubbleDofs, 
                   polOrd,numDofs,numFaceDofs,
                   rets, uRets, inv,fRhs);
      }

      // solve linear system 
      inv.solve(x,fRhs);
      
      // set new values to new velocity function 
      {
        LocalFuncType veloLF = velo.localFunction(en);
        for(int i=0; i<numDofs; ++i)
        {
          veloLF[i] = x[i];
        }
      }
    }
  }

  template <class GridPartType,
            class EntityType,
            class ArrayType, 
            class Array2Type,
            class MatrixType,
            class VectorType>
  void fillMatrix(const GridPartType& gridPart,
                  const EntityType& en,
                  const DiscreteFunctionType& uDG,
                  const FaceDiscreteSpaceType& faceSpace,
                  const ElementGradientSpaceType& gradSpace, const int startGradDofs, 
                  const ElementDiscreteSpaceType& elSpace, const int startBubbleDofs,
                  const int polOrd, const int numDofs, const int numFaceDofs,
                  ArrayType& rets, Array2Type& uRets, 
                  MatrixType& matrix, VectorType& rhs) const
  {
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;

    typedef typename FaceDiscreteSpaceType :: BaseFunctionSetType FaceBSetType  ; 
    typedef typename FaceDiscreteSpaceType :: RangeType FaceRangeType; 
    FaceRangeType faceVal;

    typedef typename DiscreteFunctionSpaceType::RangeType RangeType; 
    RangeType ret (0.0);
    RangeType neighRet (0.0);
    RangeType uPhi (0.0);

    // face quadrature type 
    //typedef CachingQuadrature<GridPartType, 1> FaceQuadratureType;
    typedef ElementQuadrature<GridPartType, 1> FaceQuadratureType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFuncType ;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      :: BaseFunctionSetType BaseFunctionSetType;

    // get uDg local on entity 
    const LocalFuncType uLF = uDG.localFunction(en);

    // get base functions set 
    const BaseFunctionSetType & bSet = uLF.baseFunctionSet(); 

    // iterate over intersections 
    IntersectionIteratorType endnit = gridPart.iend(en);
    for(IntersectionIteratorType nit = gridPart.ibegin(en);
        nit != endnit; ++nit)
    {
      const IntersectionType& inter=*nit;
      // get base function set of face 
      const FaceBSetType & faceSet =
        faceSpace.baseFunctionSet(inter.intersectionGlobal().type());
     
      const int firstRow = inter.numberInSelf() * numFaceDofs;
      
      // only interior faces are considered 
      if(inter.neighbor())
      {
        EntityPointerType neighEp = inter.outside();
        // get neighbor entity 
        const EntityType&   nb = *neighEp;
  
        // get local function of neighbor 
        const LocalFuncType uNeighLf = uDG.localFunction(nb);

        //typedef TwistUtility<GridType> TwistUtilityType;
        // for conforming situations apply Quadrature given
        //if( TwistUtilityType::conforming(gridPart.grid(),inter) )
        if( inter.conforming() )
        {
          // create quadratures 
          FaceQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
          FaceQuadratureType faceQuadOuter(gridPart, inter, polOrd, FaceQuadratureType::OUTSIDE);

          applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,
                             bSet,faceSet, uLF, uNeighLf,
                             firstRow, numFaceDofs,
                             rets,
                             ret,neighRet,faceVal,
                             matrix, rhs);
        }
        else 
        {
          // type of quadrature for non-conforming intersections 
          typedef typename FaceQuadratureType ::
            NonConformingQuadratureType NonConformingQuadratureType;
          // create quadratures 
          NonConformingQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
          NonConformingQuadratureType faceQuadOuter(gridPart, inter, polOrd, FaceQuadratureType::OUTSIDE);

          applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,
                             bSet,faceSet, uLF, uNeighLf, 
                             firstRow, numFaceDofs,
                             rets,
                             ret,neighRet,faceVal,
                             matrix, rhs);

        }
      }
     
      // only interior faces are considered 
      if(inter.boundary())
      {
        // create quadrature 
        FaceQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
        const int quadNop = faceQuadInner.nop();
        for (int l = 0; l < quadNop ; ++l)
        {
          DomainType unitNormal = 
            inter.integrationOuterNormal(faceQuadInner.localPoint(l));

          const double faceVol = unitNormal.two_norm();
          unitNormal *= 1.0/faceVol;

          // get integration weight 
          const double intel = faceVol * faceQuadInner.weight(l);

          // evaluate function 
          uLF.evaluate(faceQuadInner[l], ret);

          RangeFieldType val = ret * unitNormal; 
          val *= intel;
         
          // evaluate base functions 
          for(int i=0; i<numDofs; ++i) 
          {
            bSet.evaluate(i,faceQuadInner[l], uPhi); 
            rets[i]  = uPhi * unitNormal; 
            rets[i] *= intel;
          }

          int row = firstRow; 
          for(int j=0; j<numFaceDofs; ++j, ++row)
          {
            faceSet.evaluate(j,faceQuadInner.localPoint(l), faceVal);
            rhs[row] += val*faceVal[0];

            for(int i=0; i<numDofs; ++i) 
            {
              matrix[row][i] += (faceVal[0] * rets[i]);      
            }
          }
        }
      }
    }

    // add gradient part 
    if( gradPolOrd > 0 ) 
    {
      gradientPart(gradSpace, en, uLF, startGradDofs, uRets, matrix, rhs );
    }

    // add bubble part 
    if( bubblePolOrd > 2 ) 
    {
      bubblePart(elSpace, en, uLF, startBubbleDofs, uRets, matrix, rhs);
    }

    // printMatrix( matrix );
  }

  template <class MatrixType>
  void printMatrix(const MatrixType& matrix) const 
  {
    std::cout << "Print Matrix \n";
    for(size_t row = 0; row < matrix.N(); ++row) 
    {
      std::cout << row << ": ";
      for(size_t col = 0; col< matrix.M(); ++col) 
      {
        if( std::abs(  matrix[row][col] ) < 1e-12 ) 
          std::cout << "0 ";
        else 
          std::cout << matrix[row][col] << " ";
      }
      std::cout << std::endl; 
    }
    std::cout << "Finished print Matrix \n";
  }
      
  template <class IntersectionIteratorType, 
            class QuadratureType, 
            class BaseFunctionSetType,
            class FaceBaseFunctionSetType,
            class LocalFunctionType,
            class ArrayType,
            class RangeType,
            class FaceRangeType,
            class MatrixType,
            class RHSType>
  static void applyLocalNeighbor(const IntersectionIteratorType& nit,
                          const QuadratureType& faceQuadInner,
                          const QuadratureType& faceQuadOuter,
                          const BaseFunctionSetType& bSet,
                          const FaceBaseFunctionSetType& faceSet,
                          const LocalFunctionType& uLF,
                          const LocalFunctionType& uNeighLf,
                          const int firstRow,
                          const int numFaceDofs,
                          ArrayType& rets,
                          RangeType& ret, RangeType& neighRet,
                          FaceRangeType& faceVal,
                          MatrixType& matrix,
                          RHSType& rhs)
  {
    const typename IntersectionIteratorType::Intersection& inter = *nit;
    const int quadNop = faceQuadInner.nop();
    const int numDofs = uLF.numDofs();

    for (int l = 0; l < quadNop ; ++l)
    {
      DomainType unitNormal = 
        inter.integrationOuterNormal(faceQuadInner.localPoint(l));

      // get unit outer normal 
      const double faceVol = unitNormal.two_norm();
      unitNormal *= 1.0/faceVol;

      // integration weight 
      const double intel = faceVol * faceQuadInner.weight(l);

      // evaluate dg velocity 
      uLF.evaluate(faceQuadInner[l], ret);
      uNeighLf.evaluate(faceQuadOuter[l], neighRet);

      // take mean value 
      ret += neighRet;
      ret *= 0.5;

      RangeFieldType val = ret * unitNormal; 
      val *= intel;

      // evaluate base functions 
      for(int i=0; i<numDofs; ++i) 
      {
        bSet.evaluate(i, faceQuadInner[l], ret);
        rets[i]  = ret * unitNormal;
        rets[i] *= intel;
      }
     
      int row = firstRow; 
      for(int j=0; j<numFaceDofs; ++j, ++row )
      {
        faceSet.evaluate(j, faceQuadInner.localPoint(l), faceVal);
        rhs[row] += val * faceVal[0];

        for(int i=0; i<numDofs; ++i) 
        {
          matrix[row][i] += (faceVal[0] * rets[i]);      
        }
      }
    }
  }
 
public:
  //! application operator projection arg to H-div space 
  virtual void operator () (const DiscreteFunctionType &arg,
                            DiscreteFunctionType& dest) const 
  {
    // apply H-div projection 
    project(arg,dest);
  }

  template <class AdaptationType>
  void estimator(const DiscreteFunctionType &velo,
                 AdaptationType& adaptation) const 
  {
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;  
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType; 
    
    enum { dim = GridType :: dimension };
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 
    
    const DiscreteFunctionSpaceType& space = velo.space();
    GridPartType & gridPart = space.gridPart();
    int polOrd = space.order() + 1;

    // get local id set 
    const LocalIdSetType& localIdSet = gridPart.grid().localIdSet(); 
    
    // define type of face quadrature 
    typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;
    
    IteratorType endit = space.end();
    for(IteratorType it = space.begin(); it != endit ; ++it) 
    {
      const EntityType & en = *it;
      
      const LocalFuncType uLF = velo.localFunction(en);
      const double enVol = en.geometry().volume();
      
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      IntersectionIteratorType endnit = gridPart.iend(en);
      for(IntersectionIteratorType nit = gridPart.ibegin(en);
          nit != endnit; ++nit)
      {
        const IntersectionType& inter=*nit;
        double enError = 0.0;
        // only interior faces are considered 
        if(inter.neighbor())
        {
          EntityPointerType neighEp = inter.outside();
          EntityType&            nb = *neighEp;
          const double enVol_nbVol = 0.5 * (enVol + nb.geometry().volume());

#if HAVE_MPI 
          // get partition type 
          const bool interiorEntity = (nb.partitionType() == InteriorEntity);
#else 
          const bool interiorEntity = true;
#endif
          if( (localIdSet.id( en ) < localIdSet.id( nb )) 
#if HAVE_MPI 
             || ( ! interiorEntity ) 
#endif
            )
          {
            const LocalFuncType uNeighLf = velo.localFunction(nb);
            
            //typedef TwistUtility<GridType> TwistUtilityType;
            // for conforming situations apply Quadrature given
            //if( TwistUtilityType::conforming(gridPart.grid(),inter) )
            if( inter.conforming() )
            {
              // create quadratures 
              FaceQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
              FaceQuadratureType faceQuadOuter(gridPart, inter, polOrd, FaceQuadratureType::OUTSIDE);

              applyLocalNeighEstimator(nit,nb,faceQuadInner,faceQuadOuter,
                                 uLF, uNeighLf, enVol_nbVol, interiorEntity, 
                                 enError, adaptation);
            }
            else 
            {
              // type of quadrature for non-conforming intersections 
              typedef typename FaceQuadratureType ::
                NonConformingQuadratureType NonConformingQuadratureType;
              // create quadratures 
              NonConformingQuadratureType faceQuadInner(gridPart, inter, polOrd, FaceQuadratureType::INSIDE);
              NonConformingQuadratureType faceQuadOuter(gridPart, inter, polOrd, FaceQuadratureType::OUTSIDE);

              applyLocalNeighEstimator(nit,nb,faceQuadInner,faceQuadOuter,
                                 uLF, uNeighLf, enVol_nbVol, interiorEntity, 
                                 enError, adaptation);
            }
            
          } // end enId < nbId 
        } // end neighbor
        
        if(enError > 0.0)
        {
          adaptation.addToLocalIndicator( en , enError );
        }
      } // end intersection iterator
    }
  }

private:
  template <class IntersectionIteratorType, 
            class EntityType,
            class QuadratureType, 
            class LocalFunctionType,
            class AdaptationType>
  static inline void applyLocalNeighEstimator(const IntersectionIteratorType& nit,
                          const EntityType& nb,
                          const QuadratureType& faceQuadInner,
                          const QuadratureType& faceQuadOuter,
                          const LocalFunctionType& uLF,
                          const LocalFunctionType& uNeighLf,
                          const double enVol_nbVol,
                          const bool interiorEntity,
                          double& enError,
                          AdaptationType& adaptation)
  {
    const typename IntersectionIteratorType::Intersection& inter=*nit;
    enum { dim = GridType :: dimension };
    RangeType jump; 
    RangeType neighRet;
    double nbError = 0.0;
      
    const int quadNop = faceQuadInner.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      DomainType unitNormal = 
        inter.integrationOuterNormal(faceQuadInner.localPoint(l));
      
      double faceVol = unitNormal.two_norm();
      unitNormal *= 1.0/faceVol;

      // in case of power != 0 
      if(dim > 2) 
      {
        const double power = (1.0/(dim-1));
        faceVol = pow(faceVol, power);
      }            

      // evaluate | (u_l * n_l) + (u_r * n_r) | = | (u_l - u_r) * n_l |
      uLF.evaluate(faceQuadInner[l], jump);
      uNeighLf.evaluate(faceQuadOuter[l], neighRet);
      
      // get difference 
      jump -= neighRet; 

      const double weight = (enVol_nbVol) * faceQuadInner.weight(l);

      double error = weight * SQR(jump * unitNormal); 
      error = std::abs( error );
      
      enError += error;
      nbError += error;
    }

    if( (nbError > 0.0)
#if HAVE_MPI 
      // only add neihgbor on interior entities 
      && (interiorEntity)
#endif
      )
    {
      adaptation.addToLocalIndicator( nb , nbError );
    }
  }
};

} // end namespace 
#endif
