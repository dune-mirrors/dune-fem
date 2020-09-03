#ifndef DUNE_LOBATTOBASIS_HH
#define DUNE_LOBATTOBASIS_HH

#include <fstream>

#include <dune/geometry/type.hh>
#include <dune/geometry/topologyfactory.hh>
#include <dune/geometry/quadraturerules.hh>
// #include <dune/fem/space/localfiniteelement/lobattoquadrature.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>
#include <dune/localfunctions/lagrange/emptypoints.hh>

#include <dune/fem/common/forloop.hh>

namespace Dune
{
  namespace Impl
  {
    template< class Field >
    struct InterpolationPoints
    {
      template <class Quad>
      InterpolationPoints ( unsigned int order, const Quad &quadFactory )
      : points_( 0 ), withEndPoints_(false)
      {
        const auto &quad = quadFactory(order);
        assert( quad.size() == order+1 );
        // check if corner points are part of the quadrature - expexted to be
        // first and last point
        for (unsigned int i=0;i<=order;++i)
        {
          Field p = field_cast<Field>(quad[i].position());
          Field q = p-1.;
          if (std::abs(p)<1e-12 || std::abs(q)<1e-12)
            withEndPoints_ = true;
          else
            points_.push_back(p);
        }
      }
      std::vector<Field> points()
      {
        return points_;
      }
      bool withEndPoints() const
      {
        return withEndPoints_;
      }
      const unsigned int size() const
      {
        return points_.size();
      }
      const unsigned int order() const
      {
        return points_.size()+1;
      }
      const Field &point(int i) const
      {
        return points_[i];
      }
      private:
      std::vector<Field> points_;
      bool withEndPoints_;
    };

    template <class Field>
    struct Builder
    {
      static int size(GeometryType gt, const std::vector<Field> &points1D)
      {
        if (gt.dim()==0) return 1;
        else if (Impl::isTopology(Impl::prismConstruction,gt.id(),gt.dim()))
        {
          GeometryType baseGt(Impl::baseTopologyId(gt.id(),gt.dim()),gt.dim()-1);
          return Builder<Field>::size(baseGt,points1D)*points1D.size(); // (order-1);
        }
        else
        {
          std::cout << "Not implemented for pyramid geometries still missing!\n";
          abort();
        }
      }
      template <unsigned int dim>
      static void setup(GeometryType gt, const std::vector<Field> &points1D,
                        LagrangePoint< Field, dim > *points )
      {
        if (dim==0) return;
        if (gt.dim()==0) points->point_[0] = Zero<Field>();
        else if (Impl::isTopology(Impl::prismConstruction,gt.id(),gt.dim()))
        {
          GeometryType baseGt(Impl::baseTopologyId(gt.id(),gt.dim()),gt.dim()-1);
          const unsigned int order = points1D.size()+1;
          assert(dim>=gt.dim());
          assert(points1D.size()==order-1);
          Builder<Field>::template setup<dim>(baseGt,points1D,points);
          const unsigned int baseSize = Builder::size(baseGt,points1D);
          for (unsigned int q=0;q<points1D.size();q++)
          {
            for (unsigned int i=0;i<baseSize;++i)
            {
              const unsigned int pos = q*baseSize+i;
              for (unsigned int d=0;d<gt.dim()-1;++d)
                points[pos].point_[d] = points[i].point_[d];
              points[pos].point_[gt.dim()-1]=points1D[q];
            }
          }
        }
        else
        {
          std::cout << "Not implemented for pyramid geometries still missing!\n";
          abort();
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
    PointSetFromQuadrature(unsigned int order)
      : Base(order)
    {}

    template <class Topology, class Quad>
    bool build (const Quad& quadFactory)
    {
      unsigned int order = Base::order();
      Impl::InterpolationPoints<Field> points1D(order,quadFactory);
      if (points1D.withEndPoints())
        Dune::Fem::ForLoop<Setup<Topology>::template InitCodim,0,dimension>::
          apply(GeometryType(Topology()),order,points1D.points(),points_);
      else
        Setup<Topology>::template InitCodim<dimension>::
          apply(GeometryType(Topology()),order,points1D.points(),points_);
      return true;
    }
    static bool supports ( GeometryType gt, int order )
    {
      return gt.isCube();
    }
    template <class Topology>
    static bool supports (int order)
    {
      return supports( GeometryType(Topology()), order );
    }
    protected:
    using Base::points_;
    private:
    template <class Topology>
    struct Setup
    {
      template <int pdim>
      struct InitCodim
      {
        static const unsigned int codim = dimension-pdim;
        static void apply(GeometryType gt, const unsigned int order,
                          const std::vector<Field> &points1D,
                          std::vector<Point> &points)
        {
          const unsigned int size = Dune::Geo::Impl::size(gt.id(),gt.dim(),codim);
          for (unsigned int subEntity=0;subEntity<size;++subEntity)
          {
            GeometryType subGt(Impl::baseTopologyId(gt.id(),gt.dim(),codim),gt.dim()-codim);
            unsigned int oldSize = points.size();
            unsigned int size = Impl::Builder<Field>::size(subGt,points1D);
            if (size==0) continue;
            points.resize(oldSize+size);
            std::vector< LagrangePoint<Field,dimension-codim> > subPoints(size);
            Impl::Builder<Field>::template setup<dimension-codim>( subGt, points1D,&(subPoints[0]) );

            const GeometryType geoType( Topology::id, dimension );
            const auto &refElement = referenceElement<Field,dimension>(gt);
            const auto &mapping = refElement.template geometry< codim >( subEntity );

            LagrangePoint<Field,dimension> *p = &(points[oldSize]);
            for ( unsigned int nr = 0; nr<size; ++nr, ++p)
            {
              p->point_ = mapping.global( subPoints[nr].point_ );
              p->localKey_ = LocalKey( subEntity, codim, nr );
              #if 0 // DEBUG
              bool test = Impl::ReferenceElement<Topology,Field>::checkInside(p->point_);
              if (!test)
                std::cout << "not inside" << std::endl;
              #endif
            }
          }
        }
      };
    };
  };
  template< class F, unsigned int dim >
  struct LobattoPointSet : public PointSetFromQuadrature<F,dim>
  {
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef PointSetFromQuadrature<F,dim> Base;
    typedef typename Base::LagrangePoint Point;
    LobattoPointSet(unsigned int order)
      : Base(order)
    {}
    template <class Topology>
    bool build ()
    {
      // get LobattoQuad with order+1 points
      auto quadFactory = [](int order)
      { return Dune::QuadratureRules<Field,1>::rule(
          Dune::GeometryTypes::line, 2*order-1, Dune::QuadratureType::GaussLobatto);
      };
      return Base::template build<Topology>(quadFactory);
    }
  };
  template< class F, unsigned int dim >
  struct GaussPointSet : public PointSetFromQuadrature<F,dim>
  {
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef PointSetFromQuadrature<F,dim> Base;
    typedef typename Base::LagrangePoint Point;
    GaussPointSet(unsigned int order)
      : Base(order)
    {}
    template <class Topology>
    bool build ()
    {
      // get LobattoQuad with order+1 points
      auto quadFactory = [](int order)
      { return Dune::QuadratureRules<Field,1>::rule(
          Dune::GeometryTypes::line, 2*order+1, Dune::QuadratureType::GaussLegendre);
      };
      return Base::template build<Topology>(quadFactory);
    }
  };
}  // namespace DUNE
#endif // DUNE_LOBATTOBASIS_HH
