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

#ifdef ENABLE_UG 
#include <dune/grid/uggrid.hh>
#endif

namespace Dune {

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
  enum { bubblePolOrd = ((polOrdN - 2) < 0) ? 0 : (polOrdN - 2) };
  
  template <class Space> struct Spaces; 
  
  template <class FunctionSpaceImp,
            class GridPartImp,
            int polOrd,
            template <class> class StorageImp,
            template <class,class,int,template <class> class> class DiscreteFunctionSpaceImp>
  struct Spaces<DiscreteFunctionSpaceImp<FunctionSpaceImp,GridPartImp,polOrd,StorageImp> >
  {
    typedef DiscreteFunctionSpaceImp<ElSpaceType,GridPartImp,gradPolOrd,StorageImp> ElementGradientSpaceType;
    typedef DiscreteFunctionSpaceImp<FaceSpaceType,GridPartType,polOrdN,StorageImp> FaceDiscreteSpaceType; 
    typedef DiscreteFunctionSpaceImp<ElSpaceType,GridPartType, bubblePolOrd,StorageImp> ElementDiscreteSpaceType; 
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
  double bubbleQuad(const DomainType & point) const 
  {
    double ret = 1.0;
    for(int i=0; i<DomainType::dimension; ++i)
    {
      ret *= 4.0 * point[i] * (1.0-point[i]);
    }
    return ret;
  }
  
  // bubble function for simplex 
  double gradBubbleQuad(int deri, const DomainType & point) const 
  {
    assert( deri >= 0 );
    assert( deri < DomainType::dimension );
      
    double ret = 1.0;

    for(int i=0; i<DomainType::dimension; ++i) 
    {
      if( i != deri ) ret *= (4.0* point[i] *(1.0 - point[i]));
      else ret *= 4.0 * (1 - 2.0 *point[i]);
    }
    return ret; 
  }

  // bubble function for simplex 
  double bubbleSimplex(const DomainType & point) const 
  {
    double ret = 1.0;
    double lambdaNull = 1.0;
    for(int i=0; i<DomainType::dimension; ++i) 
    {
      ret *= point[i];
      lambdaNull -= point[i];
    }
    ret *= lambdaNull;
    return ret; 
  }
  
  // bubble function for simplex 
  double gradBubbleSimplex(int deri, const DomainType & point) const 
  {
    assert( deri >= 0 );
    assert( deri < DomainType::dimension );
      
    double ret = 1.0;
    double lambdaNull = 1.0;

    for(int i=0; i<DomainType::dimension; ++i) 
    {
      if( i != deri )  ret *= point[i];
      lambdaNull -= point[i];
    }
    ret *= lambdaNull;
    
    double first = -1.0;
    for(int i=0; i<DomainType::dimension; ++i)
    {
      first *= point[i]; 
    }
    return (ret+first); 
  }
  
  void gradientBubbleSimplex(const DomainType & point, DomainType & result) const 
  {
    for(int i=0; i<DomainType::dimension; ++i) 
    {
      result[i] = gradBubbleSimplex(i,point); 
    }
  }

  void gradientBubbleQuad(const DomainType & point, DomainType & result) const 
  {
    for(int i=0; i<DomainType::dimension; ++i) 
    {
      result[i] = gradBubbleQuad(i,point); 
    }
  }

  template <class BaseFunctionSetType, class QuadratureType, class JacoRangeType>
  void evalBubbleBase(const BaseFunctionSetType & bSet, int func, QuadratureType & quad,
      int quadPoint, JacoRangeType & result) const 
  {
    const DomainType & point = quad.point(quadPoint);
    
    bSet.jacobian(func,quad,quadPoint,result);
    //double bub = bubbleSimplex(point);
    double bub = bubbleQuad(point);
    result *= bub; 
    
    typedef typename ElementDiscreteSpaceType :: RangeType ElRangeType; 
    ElRangeType bTmp; 
    bSet.evaluate(func,quad[quadPoint], bTmp);
    RangeType tmp; 
    gradientBubbleQuad(point,tmp);
    //gradientBubbleSimplex(point,tmp);

    tmp *= bTmp[0];
    result[0]= tmp;
  }

  // only works for 2d right now 
  void curl(const DomainType & arg, DomainType & dest) const 
  {
    dest[0] =  arg[1];
    dest[1] = -arg[0];
  }

  void evalBubbleFunc(int comp, const DomainType & point, 
                      RangeType & result) const 
  {
    assert( comp == 0 );
    gradientBubbleSimplex(point,result);
  }

  template <class DiscreteFunctionSpaceImp, class EntityType, 
            class LocalFunctionType, class MatrixType, class VectorType> 
  void bubblePart(const DiscreteFunctionSpaceImp & space, 
                  EntityType & en, 
                  const LocalFunctionType & uLF, const int startRow , 
                  MatrixType & matrix, VectorType& rhs ) const 
  {
    typedef DiscreteFunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType; 

    typedef typename ElementDiscreteSpaceType::JacobianRangeType JacobianRangeType; 
    typedef typename FunctionSpaceType :: DomainType DomainType;
    
    enum { dim = GridType::dimension };
    
    const BaseFunctionSetType bSet = space.baseFunctionSet(en); 
    int polOrd = 2 * space.order() + 1;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    //const int localRows = 2 ;//bSet.numBaseFunctions(); 
    const int localRows = bSet.numBaseFunctions(); 
    //std::cout << localRows << " local bubbles \n";
    const int cols = uLF.numDofs();

    QuadratureType quad (en,polOrd);
    DomainType result;
    JacobianRangeType val;
    DomainType bVal; 
    DomainType aVal; 

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      const double intel = quad.weight(l) * 
        en.geometry().integrationElement(quad.point(l));
      
      uLF.evaluate(quad[l], result);
      for(int i=0; i<localRows; ++i)     
      {
        // we might have other row 
        int row = startRow + i;
        
        evalBubbleBase(bSet,i,quad,l,val);
        curl(val[0],aVal);

        double r = aVal * result; 
        r *= intel; 
        rhs[row] += r; 

        // for cols make matrix 
        for(int j=0; j<cols; ++j)     
        {
          uLF.baseFunctionSet().evaluate(j, quad[l], bVal);
          double t = aVal * bVal;
          t *= intel; 
          matrix[row][j] += t;
        }
      }
    }
  }

  template <class EntityType, class LocalFunctionType,
            class MatrixType, class VectorType> 
  void gradientPart(const ElementGradientSpaceType & space, 
                    EntityType & en, 
                    const LocalFunctionType & uLF, const int startRow , 
                    MatrixType& matrix, VectorType& rhs ) const 
  {
    if( space.order() <= 0 ) return ;

    typedef typename ElementGradientSpaceType::BaseFunctionSetType BaseFunctionSetType;
    const BaseFunctionSetType bSet = space.baseFunctionSet(en); 
    int polOrd = 2 * space.order() + 1;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    const int localRows = bSet.numBaseFunctions(); 
    const int cols = uLF.numDofs();

    QuadratureType quad (en,polOrd);

    RangeType result;
    RangeType uPhi;

    typedef typename ElementGradientSpaceType::JacobianRangeType GradJacobianRangeType; 
    GradJacobianRangeType gradPhi;

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      const double intel = quad.weight(l) * 
        en.geometry().integrationElement(quad.point(l));
      
      uLF.evaluate(quad[l], result);
      for(int i=0; i<localRows; ++i)     
      {
        // we might have other row 
        int row = startRow + i;
        
        bSet.jacobian(i,quad,l,gradPhi);
        double uDGVal = result * gradPhi[0];
        rhs[row] += uDGVal * intel;  

        // for cols make matrix 
        for(int j=0; j<cols; ++j)     
        {
          uLF.baseFunctionSet().evaluate(j, quad[l], uPhi);
          double val = uPhi * gradPhi[0]; 
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

    // only working for polOrd = 1 at the moment 
    //assert( space.order() == 1 );
    
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

    // only implemented for order 1 right now 
    if( space.order() > 1 )
    {
      std::cerr << std::endl;
      std::cerr << "WARNING: H-div projection not implemented for polOrd > 1 ! \n\n";
      // return doing nothing
      return ;
    }

    // colums are dofs of searched function 
    LocalFuncType lf = uDG.localFunction(*start); 
    const int numDofs = lf.numDofs();
    //std::cout << numDofs << " numDofs \n";

    const FaceBSetType faceSet = 
      GetSubBaseFunctionSet<FaceBSetType,GridType>::faceBaseSet( *start , faceSpace_ ); 
    // number of dofs on faces 
    const int numFaceDofs = faceSet.numBaseFunctions();
    
    const ElementBaseSetType elSet = elSpace_.baseFunctionSet(*start);
    const int numBubbleDofs = elSet.numBaseFunctions();
    //std::cout << numBubbleDofs << " bubbleDofs \n";
  
    const GradientBaseSetType gradSet = gradSpace_.baseFunctionSet(*start);
    // in case of linear space the is zero 
    const int numGradDofs = (space.order() <= 1) ? 0 : gradSet.numBaseFunctions();
    // std::cout << numGradDofs << " numGradDofs \n";
  
    const ReferenceElement< coordType, dim > & refElem =
        ReferenceElements< coordType, dim >::general(start->geometry().type());

    // get number of faces 
    const int overallFaceDofs = numFaceDofs * refElem.size(1);
    //std::cout << overallFaceDofs << " allFAceDofs \n";

    const int rows = (space.order() <= 1) ? (overallFaceDofs) : (overallFaceDofs + numGradDofs + numBubbleDofs);

    //const int rows = (overallFaceDofs + numGradDofs + numBubbleDofs);
    //const int rows = (space.order() <= 1) ? (overallFaceDofs) : (overallFaceDofs + numGradDofs);

    //std::cout << "faceDofs " << overallFaceDofs << " | rows " << rows << "\n";
    
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

    typedef FieldMatrix<RangeFieldType,localBlockSize,localBlockSize> FieldMatrixType;

    // resulting system matrix 
    FieldMatrixType inv;
      
    typedef FieldVector<RangeFieldType,localBlockSize> VectorType; 
    VectorType fRhs(0.0);
    VectorType x(0.0);

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
        for(int i=0; i<numDofs; ++i) 
        {
          rhs[i] = 0.0;
        }

        // fill non-symetric matrix 
        fillMatrix(gridPart_,en,uDG,faceSpace_,polOrd,numDofs,numFaceDofs,
                   rets,matrix,rhs);

        // apply least square 
        matrix.multTransposed(rhs,fakeRhs);
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
        fillMatrix(gridPart_,en,uDG,faceSpace_,polOrd,numDofs,numFaceDofs,
                   rets,inv,fRhs);
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
            class MatrixType,
            class VectorType>
  void fillMatrix(const GridPartType& gridPart,
                  const EntityType& en,
                  const DiscreteFunctionType& uDG,
                  const FaceDiscreteSpaceType& faceSpace,
                  const int polOrd, const int numDofs, const int numFaceDofs,
                  ArrayType& rets, MatrixType& matrix, VectorType& rhs) const
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

        typedef TwistUtility<GridType> TwistUtilityType;
        // for conforming situations apply Quadrature given
        if( TwistUtilityType::conforming(gridPart.grid(),inter) )
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
        bSet.evaluate(i,faceQuadInner[l], ret);
        rets[i]  = ret * unitNormal;
        rets[i] *= intel;
      }
     
      int row = firstRow; 
      for(int j=0; j<numFaceDofs; ++j, ++row )
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
            
            typedef TwistUtility<GridType> TwistUtilityType;
            // for conforming situations apply Quadrature given
            if( TwistUtilityType::conforming(gridPart.grid(),inter) )
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
