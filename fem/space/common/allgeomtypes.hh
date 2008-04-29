#ifndef DUNE_ALLGEOMTYPES_HH
#define DUNE_ALLGEOMTYPES_HH

//- system includes 
#include <vector>
#include <map>

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/grid/common/referenceelements.hh>

#define USE_UG HAVE_UG 
#if USE_UG 
#include <dune/grid/uggrid.hh>
#endif

namespace Dune
{

  /**  \brief ReferenceVolume and local bary center keeper class. 
   */
  template< class GridImp , int codim>
  class GeometryInformation  
  {
  public:
    //! grid type 
    typedef GridImp GridType;

    //! dimension 
    enum { dim = GridType :: dimension - codim };

    //! coordinate type 
    typedef typename GridType :: ctype ctype;

    //! type of reference element 
    typedef ReferenceElement< ctype, dim > ReferenceElementType; 

    //! type of domain vector 
    typedef FieldVector<ctype, dim> DomainType;

    //! map that stores the volume of the reference element  
    typedef std::map<const Dune::GeometryType, double>  ReferenceVolumeMapType;

    //! map that stores the barycenter of the reference element  
    typedef std::map<const Dune::GeometryType, DomainType>  BaryCenterMapType;

    //! type of this class 
    typedef GeometryInformation < GridType , codim > ThisType;

  protected:
    mutable BaryCenterMapType localCenters_;
    mutable ReferenceVolumeMapType referenceVolumes_;
    
    //! constructor creating empty geometry information 
    inline explicit GeometryInformation()
      : localCenters_()
      , referenceVolumes_ ()
    {
    }

  public:
    //! creating geometry information due to given geometry types list 
    GeometryInformation( const std::vector<GeometryType> & geomTypes)
      : localCenters_() 
      , referenceVolumes_()
    {
      buildMaps( geomTypes );
    }

    //! copy constructor 
    GeometryInformation( const GeometryInformation& other ) 
      : localCenters_( other.localCenters_ ) 
      , referenceVolumes_( other.referenceVolumes_ )
    {
    }

  public:  
    //! return local bary center for geometry of type type 
    const DomainType& localCenter(const GeometryType& type) const 
    {
      assert( localCenters_.find( type ) != localCenters_.end() );
      return localCenters_[ type ];
    }

    //! return volume of reference element for geometry of type type 
    const double referenceVolume(const GeometryType& type) const 
    {
      assert( referenceVolumes_.find( type ) != referenceVolumes_.end() );
      return referenceVolumes_[ type ];
    }

    //! return reference element for type 
    static const ReferenceElementType& 
      referenceElement(const GeometryType& type ) 
    {
      return ReferenceElements<ctype , dim> ::general( type );
    }

  protected:  
    //! build maps 
    void buildMaps(const std::vector<GeometryType> & geomTypes)
    {
      for(size_t i=0; i<geomTypes.size(); ++i) 
      {
        // get local bary center 
        const ReferenceElementType& refElem = referenceElement( geomTypes[i] );
        localCenters_[ geomTypes[i] ] = refElem.position(0,0);
        referenceVolumes_[ geomTypes[i] ] = refElem.volume();
      }
    }

  };


  /**  \brief default implementation uses method geomTypes of given index
       set. Used in DiscreteFunctionSpaces.
   */
  template< class IndexSetImp, class GridImp >
  class AllGeomTypes : public GeometryInformation< GridImp , 0> 
  {
  public:
    typedef IndexSetImp IndexSetType;
    typedef GridImp GridType;

  private:
    typedef AllGeomTypes< IndexSetType, GridType > ThisType;

  protected:
    const IndexSetType &indexSet_;
    
  public:
    //! constructor storing index set reference 
    inline explicit AllGeomTypes( const IndexSetType &indexSet )
      : indexSet_( indexSet )
    {
      this->buildMaps( indexSet_.geomTypes(0) );
    }

    //! returns vector with geometry tpyes this index set has indices for
    const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      return indexSet_.geomTypes( codim );
    }

    //! UGGrid might have different geom types 
    static bool multipleGeomTypes ()
    {
      return false;
    }
  };



#if USE_UG 
  /** \brief specialisation fir UGGrid, because geomTypes method of index
      sets not usable in this case. 
  */
  template< class IndexSetImp, int dimworld >
  class AllGeomTypes< IndexSetImp, UGGrid< dimworld > > 
    : public GeometryInformation< UGGrid< dimworld > , 0 > 
  {
  public:
    typedef IndexSetImp IndexSetType;
    typedef UGGrid< dimworld > GridType;

  private:
    typedef AllGeomTypes< IndexSetType, GridType > ThisType;

    CompileTimeChecker< (dimworld == 2) || (dimworld == 3) >
      __UGGrid_Works_Only_for_Dimensions_2_and_3__;
   
  protected:
    enum { ncodim = dimworld + 1 };
    
    std::vector< GeometryType > geomTypes_[ ncodim ];

  public:
    //! constructor storing index set reference 
    inline explicit AllGeomTypes( const IndexSetType &indexSet )
    {
      // vertices 
      {
        geomTypes_[dimworld].push_back( GeometryType(GeometryType::simplex,0)); 
      }
      
      {
        if( dimworld == 2 )
        {
          // elements 
          geomTypes_[0].push_back( GeometryType(GeometryType::simplex,2)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::cube,2)); 

          // faces 
          geomTypes_[1].push_back( GeometryType(GeometryType::cube,1));
        }

        if( dimworld == 3 )
        {
          // elements 
          geomTypes_[0].push_back( GeometryType(GeometryType::simplex,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::cube,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::prism,3)); 
          geomTypes_[0].push_back( GeometryType(GeometryType::pyramid,3)); 

          // faces 
          geomTypes_[1].push_back( GeometryType(GeometryType::simplex,2)); 
          geomTypes_[1].push_back( GeometryType(GeometryType::cube,2)); 

          // edges 
          geomTypes_[2].push_back( GeometryType(GeometryType::cube,1));
        }
      }

      this->buildMaps( geomTypes_[0] );
    }                  
    
    //! returns vector with geometry types this index set has indices for
    const std :: vector< GeometryType > &geomTypes ( int codim ) const
    {
      assert( (codim >= 0) && (codim < ncodim) );
      return geomTypes_[ codim ];
    }

    //! UGGrid might have different geom types 
    static bool multipleGeomTypes ()
    {
      return true;
    }
  };
#endif

#undef USE_UG

} // end namespace Dune

#endif
