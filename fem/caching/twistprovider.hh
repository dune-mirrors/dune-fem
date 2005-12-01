#ifndef DUNE_TWISTPROVIDER_HH
#define DUNE_TWISTPROVIDER_HH

//- System includes
#include <vector>
#include <map>

//- Dune includes
#include <dune/grid/common/referenceelements.hh>

//- Local includes
#include "../quadrature/quadrature.hh"

namespace Dune {

  template <class ct, int dim>
  class TwistProvider;

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
    static MapperType mappers_;
    // Must be greater than the largest negative twist possible
    static const int offset_; 
  };

  template <class ct, int dim>
  class TwistMapperCreator 
  {
  public:
    typedef Quadarture<ct, dim> QuadratureType;
    typedef FieldMatrix<ct, dim+1, dim> MatrixType;
    typedef typename Quadrature::CoordinateType PointType;
    typedef typename FieldVector<ct, dim+1> CoordinateType;

  public:
    TwistMapperCreator(const QuadratureType& quad);
         
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
 
    

  private:
    const QuadratureType& quad_;
    mutable MatrixType mat_;
    
    int minTwist_;
    int maxTwist_;

    std::vector<PointType> corners_;

    static const double eps_;
  };

} // end namespace Dune

#include "twistprovider.cc"

#endif
