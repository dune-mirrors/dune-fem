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

  template <class ct, int dim>
  class TwistProvider;
  template <class ct, int dim>
  class TwistMapperCreator;

  class TwistMapper {
    template <class ct, int dim>
    friend class TwistMapperCreator;

  public:
    TwistMapper() : indices_() {}

    size_t index(size_t quadPoint) const {
      assert(quadPoint >= 0 && quadPoint < indices_.size());
      return indices_[quadPoint];
    }

  private:
    std::vector<size_t> indices_;
  };

  template <class ct, int dim>
  class TwistProvider 
  {
  public:
    typedef Quadrature<ct, dim> QuadratureType;

  public:
    static const TwistMapper& getTwistMapper(const QuadratureType& quad,
                                             int twist); 

  private:
    typedef std::map<size_t, std::vector<TwistMapper*> > MapperType;
    typedef typename MapperType::iterator IteratorType;
    
  private:
    static IteratorType addMapper(const QuadratureType& quad);

    static std::auto_ptr<TwistMapperCreator<ct, 1> > 
    newCreator(const Quadrature<ct, 1>& quad);
    
    static std::auto_ptr<TwistMapperCreator<ct, 2> > 
    newCreator(const Quadrature<ct, 2>& quad);

  private:
    static MapperType mappers_;
    // Must be greater than the largest negative twist possible
    static const int offset_; 
  };

  template <class ct, int dim>
  class TwistMapperCreator 
  {
  public:
    typedef Quadrature<ct, dim> QuadratureType;
    typedef FieldMatrix<ct, dim+1, dim> MatrixType;
    typedef typename QuadratureType::CoordinateType PointType;
    typedef FieldVector<ct, dim+1> CoordinateType;

  public:
    TwistMapperCreator(const QuadratureType& quad,
                       int minTwist,
                       int maxTwist);
         
    TwistMapper* createMapper(int twist) const;
    
    int minTwist() const {
      return minTwist_;
    }

    int maxTwist() const {
      return maxTwist_;
    }
    
  private:
    TwistMapperCreator(const TwistMapperCreator&);
    TwistMapperCreator& operator=(const TwistMapperCreator&);
 
    bool samePoint(const PointType& first, const PointType& second) const;

    virtual void buildTransformationMatrix(int twist, 
                                           MatrixType& mat) const = 0;

  private:
    const QuadratureType& quad_;
    mutable MatrixType mat_;
    
    int minTwist_;
    int maxTwist_;

    static const ct eps_;
  };

  template <class ct>
  class LineTwistMapperCreator : public TwistMapperCreator<ct, 1> {
  public:
    enum { dim = 1 };

    typedef TwistMapperCreator<ct, dim> BaseType;
    typedef typename BaseType::QuadratureType QuadratureType;
    typedef typename BaseType::MatrixType MatrixType;

  public:
    LineTwistMapperCreator(const QuadratureType& quad);
    
  private:
    virtual void buildTransformationMatrix(int twist, 
                                           MatrixType& mat) const;
    
  private:
    const ReferenceCube<ct, dim>& refElem_;
  };

  template <class ct>
  class TriangleTwistMapperCreator : public TwistMapperCreator<ct, 2> {
  public:
    enum { dim = 2 };

    typedef TwistMapperCreator<ct, dim> BaseType;
    typedef typename BaseType::QuadratureType QuadratureType;
    typedef typename BaseType::MatrixType MatrixType;

  public:
    TriangleTwistMapperCreator(const QuadratureType& quad);

  private:
    typedef FaceTopologyMapping<hexa> FaceTopo;

  private:
    virtual void buildTransformationMatrix(int twist, 
                                           MatrixType& mat) const;

  private:
    const ReferenceSimplex<ct, dim>& refElem_;
  };

  template <class ct>
  class QuadrilateralTwistMapperCreator : public TwistMapperCreator<ct, 2> {
  public:
    enum { dim = 2 };

    typedef TwistMapperCreator<ct, dim> BaseType;
    typedef typename BaseType::QuadratureType QuadratureType;
    typedef typename BaseType::MatrixType MatrixType;

  public:   
    QuadrilateralTwistMapperCreator(const QuadratureType& quad);
    
  private:
    typedef FaceTopologyMapping<hexa> FaceTopo;

  private:
    virtual void buildTransformationMatrix(int twist, 
                                           MatrixType& mat) const;
    
  private:
    const ReferenceCube<ct, dim>& refElem_;
  };

} // end namespace Dune

#include "twistprovider.cc"

#endif
