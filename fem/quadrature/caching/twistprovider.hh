#ifndef DUNE_TWISTPROVIDER_HH
#define DUNE_TWISTPROVIDER_HH

//- System includes
#include <vector>
#include <map>
#include <memory>
#include <cassert>

//- Dune includes
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/alugrid/3d/topology.hh>
#include <dune/fem/quadrature/quadrature.hh>

//- Local includes
#include "pointmapper.hh"

namespace Dune {

  // Forward declaration
  template <class ct, int dim>
  class TwistProvider;
  template <class ct, int dim>
  class TwistMapperCreator;
  template <class ct, int dim>
  class TwistStorage;


  //! \brief Identifies quadrature points on faces with twists
  //! For a given quadrature type and a face with a given twist the TwistMapper
  //! provides a mapping from the quadrature point number on the twisted face
  //! to the quadrature point number on the (untwisted) reference face. (It
  //! removes the twist from the quadrature, so to say.) This is needed in
  //! unstructured grids when a quadrature point on a reference element's 
  //! face needs to be transformed to a quadrature point in the reference 
  //! element itself.
  //!
  //! The PointMapper objects are filled by the TwistMapperCreator and its
  //! subclasses.


  //! \brief Helper class which stores information about twist mapping for a
  //! given quadrature id.
  template <class ct, int dim>
  class TwistStorage 
  {
    typedef CachingTraits<ct, dim> Traits;
  public:
    typedef typename Traits::PointType PointType;
    typedef typename Traits::PointVectorType PointVectorType;
    typedef typename Traits::MapperType MapperType;

  public:
    //! Constructor
    //! \param minTwist Minimal possible twist for the face
    //! \param maxTwist Maximal possible twist for the face + 1
    explicit TwistStorage(int minTwist, int maxTwist);

    //! Add a new mapper for a given twist
    void addMapper(const MapperType& mapper, int twist);
    
    //! Add a point (in the case of asymmetric quadratures)
    size_t addPoint(const PointType& points);

    //! Access to a mapper
    const MapperType& getMapper(int twist) const;
    
    //! Access to the points for all twists (in the case of symmetric 
    //! quadratures, the points are identical with the quadrature points).
    const PointVectorType& getPoints() const;

    //! Minimal twist
    int minTwist() const;
    
    //! Maximal twist + 1
    int maxTwist() const;

  private:
    typedef typename Traits::MapperVectorType MapperVectorType;

  private:
    MapperVectorType mappers_;
    PointVectorType points_;

    int minTwist_;
    int maxTwist_;
  };

  //! \brief Access point for PointMapper objects with twist information
  //! PointMapper objects get created once and are reused as often as needed. 
  //! The TwistProvider serves in this context as the single point of access
  //! which is responsible for the creation and management of these objects.
  //! TwistProvider follows the monostate pattern.
  template <class ct, int dim>
  class TwistProvider 
  {
    typedef CachingTraits<ct, dim> Traits;
  public:
    //! Generic quadrature type
    typedef typename Traits::QuadratureType QuadratureType;
    //! Storage for the mappings of a specific quadrature id, alongside the
    //! resulting caching points (which may differ from the quadrature points
    //! in the case of asymmetric quadratures)
    typedef typename TwistMapperCreator<ct, dim>::TwistStorageType TwistStorageType;

  public:
    //! Delivers the PointMapper object for quadrature quad and twist twist.
    static const TwistStorageType& getTwistStorage(const QuadratureType& quad);

  private:
    typedef std::map<size_t, const TwistStorageType*> MapperContainerType;
    typedef typename MapperContainerType::iterator IteratorType;
    
  private:
    //! Gets called when a new mapper is created.
    static IteratorType createMapper(const QuadratureType& quad);

  private:
    // singleton class holding map with storages 
    class MapperContainer
    {
      // instance of map 
      MapperContainerType mappers_;
      
      //! cosntructor 
      MapperContainer() : mappers_() {}

      //! destructor  
      ~MapperContainer() 
      {
        IteratorType endit = mappers_.end();
        for(IteratorType it = mappers_.begin(); it != endit; ++it)
        {
          delete it->second;
          it->second = 0;
        }
        
      }
      
    public:
      //! return reference to mappers 
      static MapperContainerType& instance() 
      {
        // create singleton instance
        static MapperContainer mc;
        return mc.mappers_;
      }
    };
  };

  //! This class factors out all geometry dependent stuff in a strategy class
  template <class ct, int dim>
  class TwistMapperStrategy 
  {
  public:
    typedef FieldMatrix<ct, dim+1, dim> MatrixType;
    
  public:
    TwistMapperStrategy(int minTwist, int maxTwist) :
      minTwist_(minTwist),
      maxTwist_(maxTwist)
    {}

    //! virtual desctructor because of virtual functions 
    virtual ~TwistMapperStrategy() {}

    virtual const MatrixType& buildTransformationMatrix(int twist) const = 0;

    int minTwist() const { return minTwist_; }
    int maxTwist() const { return maxTwist_; }

  private:
    int minTwist_;
    int maxTwist_;
  };

  //! Helper class for TwistProvider which takes care of the creation process
  template <class ct, int dim>
  class TwistMapperCreator 
  {
    typedef CachingTraits<ct, dim> Traits;
  public:
    typedef typename Traits::QuadratureType QuadratureType;
    typedef typename Traits::PointType PointType;
    typedef typename Traits::MapperType MapperType;
    typedef FieldVector<ct, dim+1> CoordinateType;
    typedef TwistStorage<ct, dim> TwistStorageType;
    
  public:
    //! Constructor
    TwistMapperCreator(const QuadratureType& quad);

    //! Create the actual mapper for a given twist
    const TwistStorageType* createStorage() const;
    
    //! Lowest possible twist for the quadrature's geometry
    int minTwist() const {
      return helper_->minTwist();
    }

    //! Largest possible twist + 1 for the quadrature's geometry
    int maxTwist() const {
      return helper_->maxTwist();
    }
    
  private:
    typedef typename TwistMapperStrategy<ct, dim>::MatrixType MatrixType;

  private:
    TwistMapperCreator(const TwistMapperCreator&);
    TwistMapperCreator& operator=(const TwistMapperCreator&);
 
  private:
    const QuadratureType& quad_;    
    std::auto_ptr<TwistMapperStrategy<ct, dim> > helper_;

    static const ct eps_;
  };

  //! Implements the creator's functionality that depends on the underlying
  //! geometry. This is the special implementation for points.
  template <class ct, int dim>
  class PointTwistMapperStrategy : public TwistMapperStrategy<ct, dim> {
  public:
    typedef TwistMapperStrategy<ct, dim> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

  public:
    PointTwistMapperStrategy(GeometryType geo);
    
    //! virtual desctructor because of virtual functions 
    virtual ~PointTwistMapperStrategy() {}

    virtual const MatrixType& buildTransformationMatrix(int tiwst) const;
    
  private:
    const ReferenceCube<ct, dim>& refElem_;
    mutable MatrixType mat_;
  };
  //! Implements the creator's functionality that depends on the underlying
  //! geometry. This is the special implementation for line.
  template <class ct, int dim>
  class LineTwistMapperStrategy : public TwistMapperStrategy<ct, dim> {
  public:
    typedef TwistMapperStrategy<ct, dim> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

  public:
    LineTwistMapperStrategy(GeometryType geo);
    
    //! virtual desctructor because of virtual functions 
    virtual ~LineTwistMapperStrategy() {}

    virtual const MatrixType& buildTransformationMatrix(int tiwst) const;
    
  private:
    const ReferenceCube<ct, dim>& refElem_;
    mutable MatrixType mat_;
  };

  //! Implements the creator's functionality that depends on the underlying
  //! geometry. This is the special implementation for triangle.
  template <class ct, int dim>
  class TriangleTwistMapperStrategy : public TwistMapperStrategy<ct, dim> {
  public:
    typedef TwistMapperStrategy<ct, dim> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

  public:
    TriangleTwistMapperStrategy(GeometryType geo);
    //! virtual desctructor because of virtual functions 
    virtual ~TriangleTwistMapperStrategy() {}

    virtual const MatrixType& buildTransformationMatrix(int twist) const;

  private:
    typedef FaceTopologyMapping<tetra> FaceTopo;

  private:
    const ReferenceSimplex<ct, dim>& refElem_;
    mutable MatrixType mat_;
  };

  //! Implements the creator's functionality that depends on the underlying
  //! geometry. This is the special implementation for quadrilaterals
  template <class ct, int dim>
  class QuadrilateralTwistMapperStrategy : 
    public TwistMapperStrategy<ct, dim> {
  public:
    typedef TwistMapperStrategy<ct, dim> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

  public:   
    QuadrilateralTwistMapperStrategy(GeometryType geo);

    //! virtual desctructor because of virtual functions 
    virtual ~QuadrilateralTwistMapperStrategy() {}
    
    virtual const MatrixType& buildTransformationMatrix(int twist) const;
    
  private:
    typedef FaceTopologyMapping<hexa> FaceTopo;

  private:
    const ReferenceCube<ct, dim>& refElem_;
    mutable MatrixType mat_;
  };

} // end namespace Dune

#include "twistprovider.cc"
#endif
