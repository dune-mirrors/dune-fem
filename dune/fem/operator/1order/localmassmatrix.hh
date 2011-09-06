#ifndef DUNE_MASSMATRIXINVERTER_HH
#define DUNE_MASSMATRIXINVERTER_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/misc/checkgeomaffinity.hh>

namespace Dune {

/*! @addtogroup PassHyp
 *  ** @{
*/

/** \brief Local Mass Matrix for DG Operators */ 
template <class DiscreteFunctionSpaceImp, class VolumeQuadratureImp> 
class LocalDGMassMatrix 
{
public:  
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType ctype;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  enum { numDofs_ = DiscreteFunctionSpaceType :: localBlockSize };
  typedef FieldMatrix<ctype, numDofs_, numDofs_ > MatrixType;
  typedef FieldVector<ctype, numDofs_ > VectorType;

  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType BaseFunctionSetType; 
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity EntityType;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;

  typedef VolumeQuadratureImp VolumeQuadratureType;

  //! is true if grid is structured grid 
  enum { StructuredGrid = Capabilities::isCartesian< GridType >::v };

  typedef AllGeomTypes< typename GridPartType :: IndexSetType,GridType> GeometryInformationType;
  typedef typename GeometryInformationType :: DomainType DomainType;
protected:  
  const DiscreteFunctionSpaceType& spc_;

  GeometryInformationType geoInfo_;
  const int volumeQuadOrd_;
  const bool affine_;

  mutable MatrixType matrix_;
  mutable VectorType x_, rhs_;

  mutable std::vector< RangeType > phi_;
  mutable std::vector< RangeType > phiMass_;

  //! dummy caller 
  struct NoMassDummyCaller
  {
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    typedef FieldMatrix<ctype, dimRange, dimRange> MassFactorType;
    // return false since we don;t have a mass term
    bool hasMass() const { return false; }
    void mass(const EntityType&,
              const VolumeQuadratureType&,
              const int, 
              MassFactorType&) const
    {
    }
  };

public:
  //! constructor taking space and volume quadrature order 
  LocalDGMassMatrix(const DiscreteFunctionSpaceType& spc, const int volQuadOrd = -1 ) 
    : spc_(spc) 
    , geoInfo_( spc.indexSet() ) 
    , volumeQuadOrd_ ( ( volQuadOrd < 0 ) ? ( spc.order() * 2 ) : volQuadOrd )
    , affine_ ( setup() )
    , phi_( spc_.mapper().maxNumDofs() )
    , phiMass_( spc_.mapper().maxNumDofs() )
  {
  }

  //! copy constructor 
  LocalDGMassMatrix(const LocalDGMassMatrix& org) 
    : spc_(org.spc_),
      geoInfo_( org.geoInfo_),
      volumeQuadOrd_( org.volumeQuadOrd_ ),
      affine_( org.affine_ ),
      phi_( org.phi_ ),
      phiMass_( org.phiMass_ )
  {
  }

public:  
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
    if( spc_.continuous() ) 
      applyInverseContinuous( caller, en, lf );
    else 
      applyInverseDiscontinuous( caller, en, lf );
  }

  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseDiscontinuous(const MassCallerType& caller, 
                                 const EntityType& en, 
                                 LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();
    assert( numDofs_ == lf.numDofs() );

    // in case of affine mappings we only have to multiply with a factor 
    if( affine_ && ! caller.hasMass() )
    {
      const double massVolInv = getAffineMassFactor( geo );

      // apply inverse mass matrix  
      for(int l=0; l<numDofs_; ++l) 
      {
        lf[l] *= massVolInv;
      }

      return; 
    }
    else 
    {
      // copy local function to right hand side 
      for(int l=0; l<numDofs_; ++l) 
      {
        // copy 
        rhs_[l] = lf[l];
      }

      // setup local mass matrix 
      buildMatrix( caller, en, geo, lf.baseFunctionSet(), numDofs_, matrix_ );

      // solve linear system  
      matrix_.solve(  x_, rhs_ );

      // copy back to local function 
      for(int l=0; l<numDofs_; ++l)
      {
        lf[l] = x_[l];
      }

      return; 
    }
  }

  //! apply local mass matrix to local function lf
  //! using the massFactor method of the caller 
  template <class MassCallerType, class LocalFunctionType> 
  void applyInverseContinuous(const MassCallerType& caller, 
                              const EntityType& en, 
                              LocalFunctionType& lf) const 
  {
    // get geometry 
    const Geometry& geo = en.geometry();

    const int numDofs = lf.numDofs();
    DynamicVector< double > rhs( numDofs );
    DynamicVector< double > x( numDofs );

    // copy local function to right hand side 
    for(int l=0; l<numDofs; ++l) 
    {
      // copy 
      rhs[ l ] = lf[ l ];
    }

    // create matrix 
    DynamicMatrix< double > matrix( numDofs, numDofs, 0.0 );

    // setup local mass matrix 
    buildMatrix( caller, en, geo, lf.baseFunctionSet(), numDofs, matrix );

    // solve linear system  
    matrix.solve( x, rhs );

    for(int l=0; l<numDofs; ++l) 
    {
      // copy back
      lf[ l ] = x[ l ];
    }
    return; 
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

public:  
  //! returns true if geometry mapping is affine 
  bool affine () const { return affine_; }

protected:
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
    return GeometryAffinityCheck<VolumeQuadratureType> :: 
      checkAffinity( spc_.begin(), spc_.end(), volumeQuadOrd_);
  }

  //! build local mass matrix 
  template <class MassCallerType, class Matrix> 
  void buildMatrix(const MassCallerType& caller,
                   const EntityType& en,
                   const Geometry& geo, 
                   const BaseFunctionSetType& set,
                   const int numDofs,
                   Matrix& matrix) const 
  {
    assert( numDofs == set.numBaseFunctions() );

    // clear matrix 
    matrix = 0;

    // create quadrature 
    VolumeQuadratureType volQuad(en, volumeQuadOrd_ );

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
                   Matrix& matrix) const 
  {
    const int volNop = volQuad.nop();
    for(int qp=0; qp<volNop; ++qp) 
    {
      // calculate integration weight 
      const double intel = volQuad.weight(qp)
         * geo.integrationElement(volQuad.point(qp));

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
}; // end class LocalDGMassMatrix

//! @}  
} // end namespace Dune 
#endif
