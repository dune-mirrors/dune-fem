#ifndef DUNE_SIMPLEXPOINTS_HH
#define DUNE_SIMPLEXPOINTS_HH

#include <cassert>

//- Dune includes
#include <dune/common/fvector.hh>

//- local includes 
#include "ugquadratures.hh"

// include pardg quadratures 
#ifdef ENABLE_PARDG 
#include <quadrature.hpp>
#endif

namespace Dune {

  //! Adapter to the quadratures defined by UG.
  //! There is a 1-1 relationship between UGSimplexPointsAdapter object and
  //! UG quadrature implementation. 
  template <int dim>
  class UGSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    //! Constructor.
    UGSimplexPointsAdapter(int order) :
      quad_(0),
      refVol_(0.0)
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);

      switch (dim) {
      case 1:
        refVol_ = 1.0;
        break;
      case 2:
        refVol_ = 0.5;
        break;
      case 3:
        refVol_ = 1.0/6.0;
        break;
      default:
        DUNE_THROW(NotImplemented,
                   "Quadrature points for simplices only up to dim == 3");
      }
    }

    //! Number of quadrature points.
    int numPoints() const 
    {
      return quad_->nip;
    }

    //! The actual order of the quadrature.
    int order() const 
    {
      return quad_->order;
    }

    //! Access to the ith quadrature point.
    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) {
        result[j] = quad_->local[i][j];
      }
      return result;
    }

    //! Access to the ith quadrature weight.
    double weight(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      // scale with volume of reference element!
      return refVol_*quad_->weight[i];
    }
    
  private:
    const UGQuadratureType* quad_;
    double refVol_;
  };


  //! check this implementation, does not work with 
  //! cached quadrature 
  template <int dim>
  class AlbertaSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
#ifdef HAVE_ALBERTA_FOUND
    typedef ALBERTA QUAD AlbertaQuadratureType;
#else 
    struct QUAD {
      int degree;
      int n_points; 
      const double ** lambda; 
      const double * w; 
    };
    typedef QUAD AlbertaQuadratureType;
#endif
    
    typedef FieldVector<double, dim> CoordinateType;

  public:
    AlbertaSimplexPointsAdapter(int order) :
      quad_(0)
    {
#ifdef HAVE_ALBERTA_FOUND
      quad_ = get_quadrature(dim,order);  
#endif
      assert( quad_ );
    }

    int numPoints() const 
    {
      assert( quad_ );
      return quad_->n_points;
    }

    int order() const 
    {
      assert( quad_ );
      return quad_->degree;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      assert( quad_ );
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) 
      {
        // lambda is of dim+1 length, 
        // but we jsut drop the last coordinate
        result[j] = quad_->lambda[i][j];
      }
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      // no scaling, because already scaled 
      return quad_->w[i];
    }
    
  private:
    const AlbertaQuadratureType* quad_;
  };

#ifdef ENABLE_PARDG 
  //! Adapter to the quadratures defined by UG.
  //! There is a 1-1 relationship between UGSimplexPointsAdapter object and
  //! UG quadrature implementation. 
  template <int dim>
  class ParDGSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
    typedef typename pardg::Quadrature<dim> ParDGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    //! Constructor.
    ParDGSimplexPointsAdapter(int order) :
      quad_(ParDGQuadratureType::quadrature(order)),
      order_(order)
    {
      /*
      double sum = 0.0;
      for(int i=0; i<quad_.number_of_points(); ++i)
      {
        sum += quad_.w(i);
      }
      std::cout << "sum of weights " << sum << "\n"; 
      */
    }

    //! Number of quadrature points.
    int numPoints() const 
    {
      return quad_.number_of_points();
    }

    //! The actual order of the quadrature.
    int order() const 
    {
      return order_;
    }

    //! Access to the ith quadrature point.
    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) 
      {
        result[j] = quad_.x(i)[j];
      }
      return result;
    }

    //! Access to the ith quadrature weight.
    double weight(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      // scale with volume of reference element!
      return quad_.w(i);
    }
    
  private:
    const ParDGQuadratureType& quad_;
    const int order_;
  };
#endif

} // end namespace Dune

#endif
