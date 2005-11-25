#ifndef DUNE_SIMPLEXPOINTS_HH
#define DUNE_SIMPLEXPOINTS_HH

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/quadrature/fixedorder/ugquadratures.hh>

namespace Dune {

  template <int dim>
  class UGSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
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

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) {
        result[j] = quad_->local[i][j];
      }
      return result;
    }

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

  /* Not used: UGSimplexPointsAdapter does the same
  class UGTrianglePointsAdapter {
  public:
    enum { dim = 2 };
    enum { numCorners = 3 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    UGTrianglePointsAdapter(int order) :
      quad_(0) 
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);
    }

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0; i < numPoints());
      CoordinateType result;
      result[0] = quad_->local[i][0];
      result[1] = quad_->local[i][1];
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0; i < numPoints());
      // scale with volume of reference element!
      return 0.5*quad_->weight[i];
    }
    
  private:
    const UGQuadratureType* ugQuad_;
  };

  class UGTetrahedronPointsAdapter {
  public:
    enum { dim = 2 };
    enum { numCorners = 3 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    UGTetrahedronPointsAdapter(int order) :
      quad_(0) 
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);
    }

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0; i < numPoints());
      CoordinateType result;
      result[0] = quad_->local[i][0];
      result[1] = quad_->local[i][1];
      result[2] = quad_->local[i][2];
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0; i < numPoints());
      // scale with volume of reference element!
      return 0.166666666666666666666667*quad_->weight[i];
    }
  private:    
    const UGQuadratureType* ugQuad_;
  };
  */

  template <int dim>
  class AlbertaSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
#if HAVE_ALBERTA 
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
#if HAVE_ALBERTA 
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
      //std::cout << "Order = " << quad_->degree << "\n";
      return quad_->degree;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      assert( quad_ );
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) {
        // lambda is of dim+1 length, 
        // but we jsut drop the last coordinate
        result[j] = quad_->lambda[i][j];
      }
      //std::cout << result << " result \n";
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

} // end namespace Dune

#endif
