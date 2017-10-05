#include <config.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include <dune/fem/test/testgrid.hh>

#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/bindable.hh>

typedef Dune:: GridSelector::GridType HGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;

template <class GridPart>
struct A : public Dune::Fem::BindableGridFunction< GridPart, Dune::Dim<1> >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, Dune::Dim<1> > Base;
  using Base::Base;

  template <class Point>
  void evaluate(const Point &x, typename Base::RangeType &ret) const
  {
    auto xGl = Base::global(x);
    ret[0] = xGl[0]*xGl[0]*(1.-xGl[1])*(1.-xGl[1]);
  }
  unsigned int order() const { return 2; }
  std::string name() const { return "A"; } // needed for output
};
template <class GridPart>
struct B : public Dune::Fem::BindableGridFunction< GridPart, Dune::Dim<1> >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, Dune::Dim<1> > Base;
  using Base::Base;

  template <class Point>
  void evaluate(const Point &x, typename Base::RangeType &ret) const
  {
    auto xGl = Base::global(x);
    ret[0] = xGl[0]+xGl[1];
  }
  template <class Point>
  void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
  {
    ret[0][0] = 1.;
    ret[0][1] = 1.;
  }
  unsigned int order() const { return 2; }
};
template <class GFA, class GFB>
struct Difference : public Dune::Fem::BindableGridFunction< typename GFA::GridPartType, typename GFA::RangeType >
{
  typedef Dune::Fem::BindableGridFunction<GridPartType, typename GFA::RangeType> Base;
  using EntityType = typename Base::EntityType;

  Difference(const GFA &gfa, const GFB &gfb)
    : Base(gfa.gridPart()), lgfa_(gfa), lgfb_(gfb) {}
  void bind(const EntityType &entity)
  {
    Base::bind(entity);
    lgfa_.bind(entity);
    lgfb_.bind(entity);
  }
  void unbind()
  {
    Base::unbind();
    lgfa_.unbind();
    lgfb_.unbind();
  }
  template <class Point>
  void evaluate(const Point &x, typename Base::RangeType &ret) const
  {
    ret = lgfa_.evaluate(x) - lgfb_.evaluate(x);
  }
  unsigned int order() const { return std::max(lgfa_.order(),lgfb_.order()); }
  private:
  Dune::Fem::ConstLocalFunction<GFA> lgfa_;
  Dune::Fem::ConstLocalFunction<GFB> lgfb_;
};
template <class GF, int component>
struct Jacobian : public Dune::Fem::BindableGridFunction< typename GF::GridPartType, typename GF::DomainType >
{
  typedef Dune::Fem::BindableGridFunction<typename GF::GridPartType, typename GF::DomainType> Base;
  using EntityType = typename Base::EntityType;

  Jacobian(const GF &gf)
    : Base(gf.gridPart()), lgf_(gf) {}
  void bind(const EntityType &entity)
  {
    Base::bind(entity);
    lgf_.bind(entity);
  }
  void unbind()
  {
    Base::unbind();
    lgf_.unbind();
  }
  template <class Point>
  void evaluate(const Point &x, typename Base::RangeType &ret) const
  {
    ret = lgf_.jacobian(x)[component];
  }
  unsigned int order() const { return std::max(lgf_.order()-1,1u); }
  private:
  Dune::Fem::ConstLocalFunction<GF> lgf_;
};

// main program
int main(int argc, char ** argv)
{
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try
  {
    HGridType &grid = Dune::Fem::TestGrid :: grid();

    GridPartType gridPart( grid );
    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart, 2 );
    Dune::Fem::H1Norm< GridPartType > h1norm( gridPart, 2 );

    double error1 = l2norm.distance( A<GridPartType>(gridPart), B<GridPartType>(gridPart) );
    double error2 = l2norm.norm( Difference< A<GridPartType>,B<GridPartType> >
         ( A<GridPartType>(gridPart), B<GridPartType>(gridPart) ) );

    if( std::abs( error1 - error2 )  > 1e-8 )
      DUNE_THROW(Dune::InvalidStateException,"functon adapter check failed!");

    double other = l2norm.norm( Jacobian<B<GridPartType>,0>(B<GridPartType>(gridPart)) );
    double jac   = h1norm.norm( B<GridPartType>(gridPart) );
    std::cout << error1 << " " << error2 << " " << other << " " << jac << std::endl;

    A<GridPartType> a(gridPart);
    Dune::Fem::VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addVertexData(a);
    vtkWriter.pwrite("test_lfa", Dune::Fem::Parameter::commonOutputPath().c_str(),".");

    return 0;
  }
  catch( const Dune::Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
