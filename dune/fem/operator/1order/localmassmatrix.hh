#ifndef DUNE_FEM_LOCALMASSMATRIX_HH
#define DUNE_FEM_LOCALMASSMATRIX_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/typeindex.hh>
#else
#include <dune/common/geometrytypeindex.hh>
#endif

#include <dune/fem/version.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/misc/checkgeomaffinity.hh>

namespace Dune {

namespace Fem 
{

/*! @addtogroup PassHyp
 *  ** @{
*/

/** \brief Local Mass Matrix inversion implementation, select the correct method in your
    implementation */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrixImplementation  
{
public:  
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType ctype;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  enum { dgNumDofs = DiscreteFunctionSpaceType :: localBlockSize };

  typedef Dune::FieldMatrix<ctype, dgNumDofs, dgNumDofs > DGMatrixType;
  typedef Dune::FieldVector<ctype, dgNumDofs >            DGVectorType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType; 
  typedef typename IndexSetType :: IndexType   IndexType;

  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType; 

  typedef typename GridPartType :: GridType GridType;
  typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;
  typedef typename EntityType :: Geometry  Geometry;

  typedef VolumeQuadratureImp VolumeQuadratureType;

  typedef Fem::GeometryAffinityCheck<VolumeQuadratureType>  GeometryAffinityCheckType;

  //! is true if grid is structured grid 
  enum { StructuredGrid = Capabilities::isCartesian< GridType >::v };

  typedef AllGeomTypes< typename GridPartType :: IndexSetType,GridType> GeometryInformationType;
  typedef typename GeometryInformationType :: DomainType DomainType;
  // use dynamic matrix from dune-common
  typedef Dune::DynamicMatrix< RangeFieldType > MatrixType; 

protected:  
  const DiscreteFunctionSpaceType& spc_;
  const IndexSetType& indexSet_;

  GeometryInformationType geoInfo_;
  const int volumeQuadOrd_;
  const bool affine_;

  mutable DGMatrixType dgMatrix_;
  mutable DGVectorType dgX_, dgRhs_;

  // use dynamic vector from dune-common
  mutable Dune::DynamicVector< RangeFieldType > rhs_;
  mutable MatrixType matrix_;

  mutable std::vector< RangeType > phi_;
  mutable std::vector< RangeType > phiMass_;

  typedef std::map< const int, MatrixType* > MassMatrixStorageType ;
  typedef std::vector< MassMatrixStorageType > LocalInverseMassMatrixStorageType;

  mutable LocalInverseMassMatrixStorageType localInverseMassMatrix_;

  // index of entity from index set, don't setup mass matrix for the same entity twice 
  mutable IndexType lastEntityIndex_; 
  mutable unsigned int lastTopologyId_ ;
  // sequence number (obtained from DofManager via the space)
  mutable int sequence_; 

  //! dummy caller 
  struct NoMassDummyCaller
  {
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    typedef Dune::FieldMatrix<ctype, dimRange, dimRange> MassFactorType;
    // return false since we don;t have a mass term
    bool hasMass() const { return false; }
    void mass(const EntityType&,
              const VolumeQuadratureType&,
              const int, 
              MassFactorType&) const
    {
    }
  };

  template <class BaseFunctionSetType> 
  MatrixType& getLocalInverseMassMatrix(const EntityType& en, 
                                        const Geometry& geo,
                                        const BaseFunctionSetType& baseSet,
                                        const int numBaseFct ) const 
  {
    const GeometryType geomType = geo.type();
    typedef typename MassMatrixStorageType :: iterator iterator ;
    MassMatrixStorageType& massMap = localInverseMassMatrix_[ GlobalGeometryTypeIndex :: index( geomType ) ];

    iterator it = massMap.find( numBaseFct ); 
    if( it == massMap.end() ) 
    {
      MatrixType* matrix = new MatrixType( numBaseFct, numBaseFct, 0.0 );
      massMap[ numBaseFct ] = matrix ;

      // create quadrature 
      VolumeQuadratureType volQuad(en, volumeQuadratureOrder( en ) );

      // setup matrix 
      buildMatrixNoMassFactor(en, geo, baseSet, volQuad, numBaseFct, *matrix, false );
      matrix->invert();

      return *matrix;
    }

    return *((*it).second );
  }

  //! return appropriate quadrature order, default is 2 * order(entity) 
  int volumeQuadratureOrder( const EntityType& entity ) const 
  {
    return ( volumeQuadOrd_ < 0 ) ? ( spc_.order( entity ) * 2 ) : volumeQuadOrd_ ;
  }

  //! return appropriate quadrature order, default is 2 * order() 
  int maxVolumeQuadratureOrder() const 
  {
    return ( volumeQuadOrd_ < 0 ) ? ( spc_.order() * 2 ) : volumeQuadOrd_ ;
  }

public:
  //! constructor taking space and volume quadrature order 
  LocalMassMatrixImplementation(const DiscreteFunctionSpaceType& spc, const int volQuadOrd = -1 ) 
    : spc_(spc) 
    , indexSet_( spc.indexSet() )  
    , geoInfo_( indexSet_ ) 
    , volumeQuadOrd_ ( volQuadOrd )
    , affine_ ( setup() )
    , rhs_(), matrix_() 
    , phi_( spc_.mapper().maxNumDofs() )
    , phiMass_( spc_.mapper().maxNumDofs() )
    , localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) )
    , lastEntityIndex_( -1 )
    , lastTopologyId_( ~0u )
    , sequence_( -1 )  
  {
  }

  //! copy constructor 
  LocalMassMatrixImplementation(const LocalMassMatrixImplementation& org) 
    : spc_(org.spc_),
      indexSet_( spc_.indexSet() ),  
      geoInfo_( indexSet_ ),
      volumeQuadOrd_( org.volumeQuadOrd_ ),
      affine_( org.affine_ ),
      rhs_( org.rhs_ ), matrix_( org.matrix_ ),
      phi_( org.phi_ ),
      phiMass_( org.phiMass_ ),
      localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) ),
      lastEntityIndex_( org.lastEntityIndex_ ),
      lastTopologyId_( org.lastTopologyId_ ),
      sequence_( org.sequence_ )  
  {
  }

  ~LocalMassMatrixImplementation()
  {
    typedef typename MassMatrixStorageType :: iterator iterator ;
    for (unsigned int i=0;i< localInverseMassMatrix_.size();++i)
    {
      const iterator end = localInverseMassMatrix_[i].end();
      for ( iterator it =  localInverseMassMatrix_[i].begin();it != end; ++it )
        delete (*it).second;
    }
  }

public:  
  //! returns true if geometry mapping is affine 
  bool affine () const { return affine_; }

  //! return mass factor for diagonal mass matrix 
  double getAffineMassFactor(const Geometry& geo) const 
  {
    return geoInfo_.referenceVolume( geo.type() ) / geo.volume(); 
  }

  //! apply local dg mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverse(const MassCallerType& caller, 
                    const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();

    // get appropriate affinty information 
    const bool affineGeometry = ( geo.affine() ) ? true : 
        GeometryAffinityCheckType :: checkGeometry( en, geo, volumeQuadratureOrder( en ) )  ;  

    if( affineGeometry && ! caller.hasMass() ) 
    {
      applyInverseLocally( caller, en, geo, lf );
    }
    else 
    {
      applyInverseDefault( caller, en, geo, lf );
    }
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    NoMassDummyCaller caller;
    applyInverse(caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(LocalFunctionType& lf) const 
  {
    applyInverse( lf.entity(), lf );
  }

  /////////////////////////////////////////////
  // end of public methods 
  /////////////////////////////////////////////
  
  
protected:  
  ///////////////////////////////////////////////////////////
  //  applyInverse for DG spaces 
  ///////////////////////////////////////////////////////////
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseDgOrthoNormalBasis(const MassCallerType& caller, 
                                      const EntityType& en, 
                                      LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();
    assert( dgNumDofs == lf.numDofs() );

    // affine_ can be a static information 
    const bool isAffine = ( affine() ) ? true : geo.affine() ;

    // in case of affine mappings we only have to multiply with a factor 
    if( isAffine && ! caller.hasMass() )
    {
      const double massVolInv = getAffineMassFactor( geo );

      // apply inverse mass matrix  
      for(int l=0; l<dgNumDofs; ++l) 
      {
        lf[l] *= massVolInv;
      }

      return; 
    }
    else 
    {
      // copy local function to right hand side 
      for(int l=0; l<dgNumDofs; ++l) 
      {
        // copy 
        dgRhs_[ l ] = lf[ l ];
      }

      // setup local mass matrix 
      buildMatrix( caller, en, geo, lf.baseFunctionSet(), dgNumDofs, dgMatrix_ );

      // solve linear system  
      dgMatrix_.solve( dgX_, dgRhs_ );

      // copy back to local function 
      for(int l=0; l<dgNumDofs; ++l)
      {
        lf[ l ] = dgX_[ l ];
      }

      return; 
    }
  }

  //! returns true if the entity has been changed 
  bool entityHasChanged( const EntityType& entity ) const 
  {
    // don't compute matrix new for the same entity 
    const int currentSequence   = spc_.sequence();
    const unsigned int topologyId = entity.type().id();
    const IndexType entityIndex = indexSet_.index( entity ) ;

    // check whether sequence has been updated 
    if( sequence_ != currentSequence ||
        lastEntityIndex_ != entityIndex || 
        lastTopologyId_  != topologyId ) 
    {
      // update identifiers 
      lastEntityIndex_ = entityIndex ;
      sequence_        = currentSequence;
      lastTopologyId_  = topologyId ;

      return true ;
    }
    else 
      // the entity did not change
      return false ;
  }

  ///////////////////////////////////////////////////////////
  //  standard applyInverse method 
  ///////////////////////////////////////////////////////////
  //! apply local mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseDefault(MassCallerType& caller, 
                           const EntityType& entity, 
                           const Geometry& geo,
                           LocalFunctionType& lf) const 
  {
    const int numDofs = lf.numDofs();

    // if sequence changed or entity index changed 
    // compute mass matrix new 
    if( entityHasChanged( entity ) || numDofs != int(matrix_.rows()) ) 
    {
      // resize temporary memory if necessary 
      if( numDofs != int(matrix_.rows()) ) 
      {
        // resize matrix 
        matrix_.resize( numDofs, numDofs );
      }

      // setup local mass matrix 
      buildMatrix( caller, entity, geo, lf.baseFunctionSet(), numDofs, matrix_ );

      // invert mass matrix 
      matrix_.invert();
    }

    // make sure that rhs_ has the correct size
    assert( numDofs == int(matrix_.rows()) );
    // resize temporary memory if necessary 
    if( numDofs != int(rhs_.size()) )
    {
      // resize vectors 
      rhs_.resize( numDofs );
    }

    assert( int(rhs_.size()) == numDofs );

    // copy local function to right hand side 
    for(int l=0; l<numDofs; ++l) 
    {
      rhs_[ l ] = lf[ l ];
    }

    // apply inverse to right hand side and store in lf 
    multiply( numDofs, matrix_, rhs_, lf );
  }

  //! implement matvec with matrix (mv of densematrix is too stupid)
  template <class Matrix, class Rhs, class X> 
  void multiply( const int size, 
                 const Matrix& matrix, 
                 const Rhs& rhs, 
                 X& x ) const 
  {
    assert( (int) matrix.rows() == size );
    assert( (int) matrix.cols() == size );
    assert( (int) rhs.size() == size );

    for( int row = 0; row < size; ++ row ) 
    {
      RangeFieldType sum = 0;
      // get matrix row 
      typedef typename Matrix :: const_row_reference  MatRow;
      MatRow& matRow = matrix[ row ];

      // multiply row with right hand side
      for( int col = 0; col < size; ++ col ) 
      {
        sum += matRow[ col ] * rhs[ col ];        
      }

      // set to result to result vector 
      x[ row ] = sum;
    }
  }

  ///////////////////////////////////////////////////////////
  //  local applyInverse method for affine geometries 
  ///////////////////////////////////////////////////////////
  //! apply local mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseLocally(MassCallerType& caller, 
                           const EntityType& en, 
                           const Geometry& geo,
                           LocalFunctionType& lf) const 
  {
    const int numDofs = lf.numDofs();

    // get local inverted mass matrix 
    MatrixType& invMassMatrix = 
      getLocalInverseMassMatrix( en, geo, lf.baseFunctionSet(), numDofs );

    // resize vectors 
    rhs_.resize( numDofs );

    const double massVolInv = getAffineMassFactor( geo );

    // copy local function to right hand side 
    // and apply inverse mass volume fraction 
    for(int l=0; l<numDofs; ++l) 
    {
      rhs_[ l ] = lf[ l ] * massVolInv ;
    }

    // apply inverse local mass matrix and store in lf 
    multiply( numDofs, invMassMatrix,  rhs_, lf );
  }

  //! setup and return affinity 
  bool setup() const 
  {
    // for structured grids this is always true 
    if( StructuredGrid ) return true;

    // get types for codim 0 
    const std::vector<Dune::GeometryType>& geomTypes = geoInfo_.geomTypes(0);

    // for simplices we also have affine mappings 
    if( (geomTypes.size() == 1) && geomTypes[0].isSimplex() ) 
    {
      return true;
    }

    bool affine = checkAffinity();
    //std::cout << "found affine grid? = " << affine << "\n";
    return affine;
  }

  //! check whether all geometry mappings are affine 
  bool checkAffinity() const 
  {
    return GeometryAffinityCheckType :: 
      checkAffinity( spc_.begin(), spc_.end(), maxVolumeQuadratureOrder() );
  }

  //! build local mass matrix 
  template <class MassCallerType, class Matrix> 
  void buildMatrix(const MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const std::size_t numDofs,
                   Matrix& matrix) const 
  {
    assert( numDofs == set.size() );

    // clear matrix 
    matrix = 0;

    // create quadrature 
    VolumeQuadratureType volQuad(en, volumeQuadratureOrder( en ) );

    if( caller.hasMass() )
    {
      // build matix with calling 
      buildMatrixWithMassFactor(caller, en, geo, set, volQuad, numDofs, matrix);
    }
    else 
    {
      buildMatrixNoMassFactor(en, geo, set, volQuad, numDofs, matrix);
    }
  }

  //! build local mass matrix with mass factor 
  template <class Matrix>
  void buildMatrixNoMassFactor(
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   const int numDofs,
                   Matrix& matrix,
                   const bool applyIntegrationElement = true ) const 
  {
    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = ( applyIntegrationElement ) ? 
          ( volQuad.weight(qp) * geo.integrationElement(volQuad.point(qp)) ) : volQuad.weight(qp) ;

      // eval base functions 
      set.evaluateAll(volQuad[qp], phi_);

      for(int m=0; m<numDofs; ++m)
      {
        const RangeType& phi_m = phi_[m];
        const ctype val = intel * (phi_m * phi_m);
        matrix[m][m] += val;
       
        for(int k=m+1; k<numDofs; ++k) 
        {
          const ctype val = intel * (phi_m * phi_[k]);
          matrix[m][k] += val;
          matrix[k][m] += val;
        }
      }
    }
  }

  //! build local mass matrix with mass factor 
  template <class MassCallerType, class Matrix> 
  void buildMatrixWithMassFactor(
                   const MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const VolumeQuadratureType& volQuad,
                   const int numDofs,
                   Matrix& matrix) const 
  {
    typedef typename MassCallerType :: MassFactorType MassFactorType;
    MassFactorType mass;

    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = volQuad.weight(qp)
         * geo.integrationElement(volQuad.point(qp));

      // eval base functions 
      set.evaluateAll( volQuad[qp], phi_);

      // call mass factor 
      caller.mass( en, volQuad, qp, mass);

      // apply mass matrix to all base functions 
      for(int m=0; m<numDofs; ++m)
      {
        mass.mv( phi_[m], phiMass_[m] );
      }

      // add values to matrix 
      for(int m=0; m<numDofs; ++m)
      {
        for(int k=0; k<numDofs; ++k) 
        {
          matrix[m][k] += intel * (phiMass_[m] * phi_[k]);
        }
      }
    }
  }
}; // end class LocalMassMatrixImplementation

///////////////////////////////////////////////////////////////////
//
//  default LocalMassMatrix Implementation 
//
///////////////////////////////////////////////////////////////////

/** \brief Local Mass Matrix for arbitrary spaces */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrix 
  : public LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp >
{
  typedef LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  LocalMassMatrix( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}
};


///////////////////////////////////////////////////////////////////
//
//  DG LocalMassMatrix Implementation 
//
///////////////////////////////////////////////////////////////////

/** \brief DG Local Mass Matrix for arbitrary spaces */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalMassMatrixImplementationDgOrthoNormal
  : public LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp >
{
  typedef LocalMassMatrixImplementation< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
public:
  typedef typename BaseType :: EntityType  EntityType;

  LocalMassMatrixImplementationDgOrthoNormal( const DiscreteFunctionSpaceImp& spc, const int volQuadOrd = -1 )
    : BaseType( spc, volQuadOrd )
  {}

  //! apply local dg mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverse(const MassCallerType& caller, 
                    const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    BaseType :: applyInverseDgOrthoNormalBasis( caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(const EntityType& en, 
                    LocalFunctionType& lf) const 
  {
    typename BaseType :: NoMassDummyCaller caller;
    applyInverse(caller, en, lf );
  }

  //! apply local dg mass matrix to local function lf without mass factor 
  template <class LocalFunctionType> 
  void applyInverse(LocalFunctionType& lf) const 
  {
    applyInverse( lf.entity(), lf );
  }

};

//! @}  

} // end namespace Fem  

// #if DUNE_FEM_COMPATIBILITY  
// put this in next version 1.4 

using Fem :: LocalMassMatrix ;

// #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune 
#endif
