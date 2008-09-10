#ifndef DUNE_GEOMETRYCONVERSION_HH
#define DUNE_GEOMETRYCONVERSION_HH

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <dune/common/geometrytype.hh>

namespace Dune { 

  /** @ingroup HelperClasses 
     \brief Identifier for geometry types apart from dune::GeometryType */
  class GeometryIdentifier {
  public:

    // this enum specifies the highest number of available geometry types 
    enum { numTypes = 10 };

    // in GeometryType 
    // simplex = 0 
    // cube    = 1
    // pyramid = 2 
    // prism   = 3
    //
    // formula: id = (2 * type) + dimension 
    enum IdentifierType { 
      Vertex        = 0, 
      Line          = 1,
      Triangle      = 2 * (GeometryType::simplex) + 2, // = 2
      Tetrahedron   = 2 * (GeometryType::simplex) + 3, // = 3 
      Quadrilateral = 2 * (GeometryType::cube)    + 2, // = 4 
      Hexahedron    = 2 * (GeometryType::cube)    + 3, // = 5 
      Pyramid       = 2 * (GeometryType::pyramid) + 3, // = 7
      Prism         = 2 * (GeometryType::prism)   + 3, // = 9
      Unknown       = -1 };

    struct CheckNumbers
    {
      static void check() 
      {
        assert(Vertex        == 0);
        assert(Line          == 1);
        assert(Triangle      == 2);
        assert(Tetrahedron   == 3);
        assert(Quadrilateral == 4);
        assert(Hexahedron    == 5);
        assert(Pyramid       == 7);
        assert(Prism         == 9);
      }
    };

    // convert GeometryType::BasicType to IdentifierType by given formula above
    template <int cd,int dim> 
    struct FromGeometry
    {
      inline static IdentifierType convert (const GeometryType :: BasicType type)
      {
        return static_cast<IdentifierType> ((2 * type) + dim);       
      }
    };

    // dim = 0 is Vertex 
    template <int cd> 
    struct FromGeometry<cd,0> 
    { 
      inline static IdentifierType convert (const GeometryType :: BasicType type) 
      { 
        return Vertex; 
      }
    };
    
    // dim = 1 is Line 
    template <int cd> 
    struct FromGeometry<cd,1> 
    { 
      inline static IdentifierType convert (const GeometryType :: BasicType type) 
      { 
        return Line; 
      }
    };

  public:
    
    GeometryIdentifier(IdentifierType idType) :
      identifier_(idType) {}

    GeometryIdentifier(int dimension, const GeometryType & geoType) :
      identifier_(GeometryIdentifier::fromGeo(dimension, geoType)) {}

    inline operator GeometryType() const { 
      return GeometryIdentifier::toGeo(identifier_); 
    }
    inline operator IdentifierType() const { return identifier_; }

    static inline GeometryType toGeo(IdentifierType id) {
      switch(id) {
      case Vertex:
        return GeometryType(GeometryType::simplex,0);
      case Line:
      {
        //lines are simplices, more than cubes ;)
        return GeometryType(GeometryType::cube,1);
      }
      case Triangle:
        return GeometryType(GeometryType::simplex,2);
      case Tetrahedron:
        return GeometryType(GeometryType::simplex,3);
      case Quadrilateral:
        return GeometryType(GeometryType::cube,2);
      case Hexahedron:
        return GeometryType(GeometryType::cube,3);
      case Pyramid:
        return GeometryType(GeometryType::pyramid,3);
      case Prism:
        return GeometryType(GeometryType::prism,3);
      default:
        {
          assert(false);
          DUNE_THROW(NotImplemented,"GeometryType not available");
          // vertex is the new unknown ;) 
          return GeometryType(GeometryType::simplex,0);
        }
      }
    }

    //! return conversion from geometry to identifier 
    static inline IdentifierType fromGeo(const GeometryType & geo) {
      return fromGeo(geo.dim(),geo);
    }

    //! return conversion from geometry to identifier 
    static inline IdentifierType fromGeo(int dimension, const GeometryType & geo) 
    {
      
      switch(dimension) {
        case 0: return Vertex; 
        case 1: return Line; 
        case 2: 
                {
                  if( geo.isSimplex() ) return Triangle;
                  if( geo.isCube() )    return Quadrilateral;
                  
                  std::cerr<<"Wrong GeometryType in fromGeom in " << __FILE__ << " " << __LINE__ << "\n";
                  DUNE_THROW(NotImplemented,"GeometryType not implemented");
                  abort();
                  return Unknown;
                }
        case 3 : 
                {
                  if(geo.isSimplex()) return Tetrahedron;
                  if(geo.isCube()   ) return Hexahedron;
                  if(geo.isPyramid()) return Pyramid;
                  if(geo.isPrism()  ) return Prism;

                  std::cerr<<"Wrong GeometryType in fromGeom in " << __FILE__ << " " << __LINE__ << "\n";
                  DUNE_THROW(NotImplemented,"GeometryType not implemented");
                  abort();
                  return Unknown;
                }
      default:
        std::cerr<<"Wrong GeometryType in fromGeom in " << __FILE__ << " " << __LINE__ << "\n";
        DUNE_THROW(NotImplemented,"GeometryType not implemented");
        abort();
        return Unknown;
      }
    }

    template <class GeometryImp> 
    inline static IdentifierType fromGeometry (const GeometryImp & geo) 
    {
#ifndef NDEBUG
      CheckNumbers::check();
      IdentifierType id = FromGeometry<0,GeometryImp::mydimension>::convert(geo.type().basicType());
      assert(id >= 0);
      assert(id == fromGeo(geo.type())); 
#endif
      return FromGeometry<0,GeometryImp::mydimension>::convert(geo.type().basicType());
    }

  private:
    IdentifierType identifier_;
  };

} // end namespace Dune 

#endif
