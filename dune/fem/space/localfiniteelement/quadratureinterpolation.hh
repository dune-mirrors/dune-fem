#ifndef DUNE_LOBATTOBASIS_HH
#define DUNE_LOBATTOBASIS_HH

#include <fstream>

#include <dune/geometry/type.hh>
#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/fem/common/forloop.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>
#include <dune/localfunctions/lagrange/emptypoints.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>


namespace Dune
{
  namespace Impl
  {
    template <class Field>
    struct Builder
    {
      template <class Points1DType>
      static int size(const GeometryType gt, const Points1DType &points1D)
      {
        if (gt.dim()==0) return 1;
        else if (gt.isPrismatic())
          return Builder<Field>::size(Impl::getBase(gt),points1D)*points1D.size(); // (order-1);
        else
        {
          std::cout << "Not implemented for pyramid geometries still missing!\n";
          std::abort();
        }
      }
      template <unsigned int dim, class Points1DType>
      static void setup(const GeometryType gt, const Points1DType &points1D,
                        LagrangePoint< Field, dim > *points )
      {
        if (dim==0)
        {
          points->weight_ = 1.;
          return;
        }
        if (gt.dim()==0)
        {
          points->point_[0] = Zero<Field>();
          points->weight_ = 1.;
        }
        else if (gt.isPrismatic())
        {
          GeometryType baseGt = Impl::getBase(gt);
          assert(dim>=gt.dim());
          Builder<Field>::template setup<dim>(baseGt,points1D,points);
          const unsigned int baseSize = Builder::size(baseGt,points1D);
          for (unsigned int i=0;i<baseSize;++i)
          {
            Field weight = points[i].weight_;
            for (unsigned int q=0;q<points1D.size();q++)
            {
              const unsigned int pos = q*baseSize+i;
              for (unsigned int d=0;d<gt.dim()-1;++d)
                points[pos].point_[d] = points[i].point_[d];
              points[pos].point_[gt.dim()-1]=points1D[q].first;
              points[pos].weight_ = weight*points1D[q].second;
            }
          }
        }
        else
        {
          std::cout << "Not implemented for pyramid geometries still missing!\n";
          std::abort();
        }
      }
    };
  } // namespace Impl

  template< class F, unsigned int dim >
  struct PointSetFromQuadrature : public EmptyPointSet<F,dim>
  {
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef EmptyPointSet<F,dim> Base;
    typedef typename Base::LagrangePoint Point;
    typedef std::vector<std::pair<Field,Field>> Points1DType;
    PointSetFromQuadrature(unsigned int order)
      : Base(order), quadOrder_(-1)
    {}

    template <GeometryType::Id geometryId, class Quad>
    bool build (const Quad& quadFactory)
    {
      constexpr GeometryType gt(geometryId);
      unsigned int order = Base::order();
      const auto &quad = quadFactory(order);
      quadOrder_ = quad.order();
      assert( quad.size() == order+1 );
      bool withEndPoints = false;
      Points1DType points1D;
      Field vertexWeight = 1;
      for (unsigned int i=0;i<=order;++i)
      {
        // remove corner points if part of the quadrature - are added in by
        // topology construction of points
        Field p = field_cast<Field>(quad[i].position());
        Field q = p-1.;
        if (std::abs(p)<1e-12 || std::abs(q)<1e-12)
        {
          withEndPoints = true;
          vertexWeight = quad[i].weight(); // assuming weight is identical for both end points
        }
        else
          points1D.push_back(std::make_pair(p,quad[i].weight()));
      }
      if (withEndPoints)
        Dune::Fem::ForLoop<Setup::template InitCodim,0,dimension>::
          apply(gt,order,points1D,vertexWeight,points_);
      else
        Setup::template InitCodim<dimension>::
          apply(gt,order,points1D,vertexWeight,points_);
      return true;
    }
    static bool supports ( GeometryType gt, int order )
    {
      return gt.isCube();
    }
    template< GeometryType::Id geometryId>
    static bool supports ( std::size_t order ) {
      return supports( GeometryType( geometryId ), order );
    }
    unsigned int quadOrder() const
    {
      return quadOrder_;
    }
    protected:
    using Base::points_;
    unsigned int quadOrder_;
    private:
    struct Setup
    {
      template <int pdim>
      struct InitCodim
      {
        static const unsigned int codim = dimension-pdim;
        static void apply(GeometryType gt, const unsigned int order,
                          const Points1DType &points1D,
                          const Field &vertexWeight,
                          std::vector<Point> &points)
        {
          const unsigned int subEntities = Dune::Geo::Impl::size(gt.id(),gt.dim(),codim);
          for (unsigned int subEntity=0;subEntity<subEntities;++subEntity)
          {
            GeometryType subGt(Impl::baseTopologyId(gt.id(),gt.dim(),codim),gt.dim()-codim);
            unsigned int oldSize = points.size();
            unsigned int size = Impl::Builder<Field>::size(subGt,points1D);
            if (size==0) continue;
            points.resize(oldSize+size);
            std::vector< LagrangePoint<Field,dimension-codim> > subPoints(size);
            Impl::Builder<Field>::template setup<dimension-codim>( subGt, points1D, &(subPoints[0]) );

            const auto &refElement = referenceElement<Field,dimension>(gt);
            const auto &mapping = refElement.template geometry< codim >( subEntity );

            LagrangePoint<Field,dimension> *p = &(points[oldSize]);
            for ( unsigned int nr = 0; nr<size; ++nr, ++p)
            {
              p->point_    = mapping.global( subPoints[nr].point_ );
              p->weight_   = subPoints[nr].weight_ * std::pow(vertexWeight,codim)*refElement.volume();
              p->localKey_ = LocalKey( subEntity, codim, nr );
            }
          }
        }
      };
    };
  };

  ///////////////////////////////////////////////////////
  //  a point set from dune-localfunctions
  ///////////////////////////////////////////////////////
  template< class F, unsigned int dim >
  struct EquidistantPointSetDerived : public Dune::EquidistantPointSet< F, dim >
  {
    typedef Dune::EquidistantPointSet< F, dim > Base;

    EquidistantPointSetDerived(unsigned int order)
      : Base(order)
    {}

    // needed for InterpolationQuadrature
    static const int pointSetId = Dune::QuadratureType::size;

    static unsigned int quad2PolOrder(int order) { return order; }
    unsigned int quadOrder() const { return Base::order(); }

    static auto buildCubeQuadrature(unsigned int quadOrder)
    {
      using namespace Impl;
      EquidistantPointSetDerived ps(quad2PolOrder(quadOrder));
      ps.template build<GeometryTypes::cube(dim)>();
      return ps;
    }
  };


  ///////////////////////////////////////////////////////
  // Some point sets from dune-geometry
  ///////////////////////////////////////////////////////
  template< class F, unsigned int dim >
  struct GaussLobattoPointSet : public PointSetFromQuadrature<F,dim>
  {
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef PointSetFromQuadrature<F,dim> Base;
    typedef typename Base::LagrangePoint Point;

    // enum identifier from dune-geometry QuadratureRules
    static const int pointSetId = Dune::QuadratureType::GaussLobatto;

    GaussLobattoPointSet(unsigned int order)
      : Base(order)
    {}
    template< GeometryType::Id geometryId >
    bool build ()
    {
      // get LobattoQuad with order+1 points
      auto quadFactory = [](int order)
      { return Dune::QuadratureRules<Field,1>::rule(
          Dune::GeometryTypes::line, pol2QuadOrder(order),
                    Dune::QuadratureType::GaussLobatto);
      };
      return Base::template build<geometryId>(quadFactory);
    }
    static unsigned int pol2QuadOrder(int order)
    {
      return (order>0)? 2*order-1 : 0;
    }
    static unsigned int quad2PolOrder(int order)
    {
      return order/2 + 1;
    }

    static auto buildCubeQuadrature(unsigned int quadOrder)
    {
      using namespace Impl;
      GaussLobattoPointSet ps(quad2PolOrder(quadOrder));
      ps.template build<GeometryTypes::cube(dim)>();
      return ps;
    }
  };


  template< class F, unsigned int dim >
  struct GaussLegendrePointSet : public PointSetFromQuadrature<F,dim>
  {
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef PointSetFromQuadrature<F,dim> Base;
    typedef typename Base::LagrangePoint Point;

    // enum identifier from dune-geometry QuadratureRules
    static const int pointSetId = Dune::QuadratureType::GaussLegendre;

    GaussLegendrePointSet(unsigned int order)
      : Base(order)
    {}
    template< GeometryType::Id geometryId >
    bool build ()
    {
      // get LobattoQuad with order+1 points
      auto quadFactory = [](int order)
      { return Dune::QuadratureRules<Field,1>::rule(
          Dune::GeometryTypes::line, pol2QuadOrder(order), Dune::QuadratureType::GaussLegendre);
      };
      return Base::template build<GeometryType(geometryId)>(quadFactory);
    }

    static unsigned int pol2QuadOrder(int order)
    {
      return 2*order+1;
    }
    static unsigned int quad2PolOrder(int order)
    {
      return order/2;
    }

    static auto buildCubeQuadrature(unsigned int quadOrder)
    {
      using namespace Impl;
      GaussLegendrePointSet ps(quad2PolOrder(quadOrder));
      ps.template build<GeometryTypes::cube(dim)>();
      return ps;
    }
  };

  ///////////////////////////////////////////////////////
  // Some point sets for preconditioning providing
  // interpolation based on cell centers of refined cells
  ///////////////////////////////////////////////////////
  template< class F, unsigned int dim >
  struct CellCentersPointSet : public PointSetFromQuadrature<F,dim>
  {
    typedef PointSetFromQuadrature<F,dim> Base;
    static const unsigned int dimension = dim;
    typedef F Field;

    // needed for InterpolationQuadrature
    static const int pointSetId = Dune::QuadratureType::size+1;

    CellCentersPointSet(unsigned int order)
      : Base(order)
    {}

    /** A list of points formed by the mid points of
     *  p+1 equispaced intervals between 0 and 1.
     */
    template<typename ct>
    class CellCenters : public Dune::QuadratureRule<ct,1>
    {
    public:
      /** brief The highest quadrature order available */
      enum { highest_order=7 };

      friend class QuadratureRuleFactory<ct,1>;
      CellCenters(int p)
      {
        this->delivered_order = p;

        const int n = p+1; // we have p+1 intervals
        const Field h = 1./Field(n);
        // compute cell centers of the interval
        for( int i=0; i<n; ++i )
        {
          const Field x = 0.5*h + h*i;
          assert( x > 0.0 && x < 1.0 );
          this->push_back(QuadraturePoint<ct,1>(x, h));
        }
      }
    };

    static unsigned int quad2PolOrder(int order) { return order; }

    unsigned int quadOrder() const { return Base::order(); }

    template< GeometryType::Id geometryId >
    bool build ()
    {
      // get FV centers as points
      auto quadFactory = [](int order)
      {
        return CellCenters<Field>(order);
      };
      return Base::template build<geometryId>(quadFactory);
    }
    static auto buildCubeQuadrature(unsigned int quadOrder)
    {
      using namespace Impl;
      CellCentersPointSet ps(quad2PolOrder(quadOrder));
      ps.template build<GeometryTypes::cube(dim)>();
      return ps;
    }
  };

}  // namespace DUNE

#endif // HAVE_DUNE_LOCALFUNCTIONS

#endif // DUNE_LOBATTOBASIS_HH
