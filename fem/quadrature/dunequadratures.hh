#ifndef DUNE_DUNEQUADRATURES_HH
#define DUNE_DUNEQUADRATURES_HH

//- Dune includes 
#include <dune/grid/common/quadraturerules.hh>

namespace Dune {

  template <class ct, int dim>
  class QuadratureImp;


  template< typename ct, int dim >
  class QuadratureRulesFactory :
     public QuadratureImp<ct,dim>
  {     
    typedef QuadratureImp<ct,dim>         BaseType;
    typedef QuadratureRulesFactory<ct,dim> ThisType;

    typedef QuadratureRule<ct,dim> DuneQuadratureRuleType;

    enum { highest_order_cube    = CubeQuadratureRule<ct,dim>::highest_order };
    enum { highest_order_simplex = SimplexQuadratureRule<ct,dim>::highest_order };

    enum { highest_order = 44 };//(highest_order_cube < highest_order_simplex) ? 
//                highest_order_cube : highest_order_simplex };

    const GeometryType elementGeometry_;
    int order_;
  public:
   QuadratureRulesFactory(
              const GeometryType& geo,
              const int order,
              const size_t id)
      : BaseType(id)
      , elementGeometry_(geo)
    {
      // get gauss quadrature 
      const DuneQuadratureRuleType& rule = 
        QuadratureRules<ct,dim> :: rule(geo,order, QuadratureType::Gauss); 
      assert( order <= rule.order());

      order_ = rule.order();

      typedef typename DuneQuadratureRuleType :: iterator iterator;
      iterator endit = rule.end();
      for(iterator it = rule.begin(); it != endit; ++it)
      {
        this->addQuadraturePoint((*it).position(),(*it).weight());
      }
    }

    //! return order of points set 
    int order () const { return order_; }
    //! return geometry type set was created for 
    GeometryType geometry() const { return elementGeometry_; }

    // return max order 
    static int maxOrder () { return highest_order; }

  };

} // end namespace Dune 
#endif
