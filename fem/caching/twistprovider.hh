#ifndef DUNE_TWISTPROVIDER_HH
#define DUNE_TWISTPROVIDER_HH

//- System includes
#include <vector>
#include <map>
#include <memory>
#include <cassert>

//- Dune includes
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/alu3dgrid/topology.hh>

//- Local includes
#include "../quadrature/quadrature.hh"

namespace Dune {

  // Forward declaration
  template <class ct, int dim>
  class TwistProvider;
  template <class ct, int dim>
  class TwistMapperCreator;

  //! \brief Identifies quadrature points on faces with twists
  //! For a given quadrature type and a face with a given twist the TwistMapper
  //! provides a mapping from the quadrature point number on the twisted face
  //! to the quadrature point number on the (untwisted) reference face. (It
  //! removes the twist from the quadrature, so to say.) This is needed in
  //! unstructured grids when a quadrature point on a reference element's 
  //! face needs to be transformed to a quadrature point in the reference 
  //! element itself.
  //!
  //! The TwistMapper objects are filled by the TwistMapperCreator and its
  //! subclasses.
  class TwistMapper {
  private:
    template <class ct, int dim>
    friend class TwistMapperCreator;

    //! Default constructor
    TwistMapper(size_t size) : indices_(size) {}

  public:
    //! The mapping from quadrature point indices on the twisted face to
    //! the quadrature point indices on the (untwisted) reference face
    size_t index(size_t quadPoint) const {
      assert(quadPoint >= 0 && quadPoint < indices_.size());
      return indices_[quadPoint];
    }

  private:
    std::vector<size_t> indices_;
  };

  //! \brief Access point for TwistMapper
  //! TwistMapper object get created once and are reused as often as needed. 
  //! The TwistProvider serves in this context as the single point of access
  //! which is responsible for the creation and management of these objects.
  //! TwistProvider follows the monostate pattern.
  template <class ct, int dim>
  class TwistProvider 
  {
  public:
    //! Generic quadrature type
    typedef Quadrature<ct, dim> QuadratureType;

  public:
    //! Delivers the TwistMapper object for quadrature quad and twist twist.
    static const TwistMapper& getTwistMapper(const QuadratureType& quad,
                                             int twist); 

  private:
    typedef std::map<size_t, std::vector<TwistMapper*> > MapperType;
    typedef typename MapperType::iterator IteratorType;
    
  private:
    //! Gets called when a new mapper is created.
    static IteratorType addMapper(const QuadratureType& quad);

  private:
    static MapperType mappers_;
    // Must be greater than the largest negative twist possible
    static const int offset_; 
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
  public:
    typedef Quadrature<ct, dim> QuadratureType;
    typedef typename QuadratureType::CoordinateType PointType;
    typedef FieldVector<ct, dim+1> CoordinateType;

  public:
    //! Constructor
    TwistMapperCreator(const QuadratureType& quad);

    //! Create the actual mapper for a given twist
    TwistMapper* createMapper(int twist) const;
    
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
 
    //! Are two points the same (floating point comparison with tolerance eps)
    bool samePoint(const PointType& first, const PointType& second) const;

  private:
    const QuadratureType& quad_;    
    std::auto_ptr<TwistMapperStrategy<ct, dim> > helper_;

    static const ct eps_;
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

    virtual const MatrixType& buildTransformationMatrix(int twist) const;

  private:
    typedef FaceTopologyMapping<hexa> FaceTopo;

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
