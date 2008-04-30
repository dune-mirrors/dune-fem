#ifndef DUNE_FEM_BOUNDARYMODEL_HH
#define DUNE_FEM_BOUNDARYMODEL_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  class DefaultBoundaryModelProperties
  {
  public:
    enum { hasDirichletValues = true };
    enum { hasNeumannValues = true };
    enum { hasRobinValues = true };
  };



  template< class FunctionSpaceImp,
            class BoundaryModelImp,
            class PropertiesImp >
  class BoundaryModelInterface
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef PropertiesImp Properties;

  private:
    typedef BoundaryModelInterface< FunctionSpaceType, BoundaryModelImp, Properties >
      ThisType;

  public:
    typedef ThisType BoundaryModelInterfaceType;

    enum BoundaryType { Dirichlet, Neumann, Robin };

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    /** \brief obtain type of boundary intersection
     *
     *  \param[in] intersection intersection whose boundary type shall be returned
     *
     *  \returns boundary type of the intersection
     */
    template< class IntersectionIteratorType >
    inline BoundaryType boundaryType ( const IntersectionIteratorType &intersection ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().boundaryType( intersection ) );
      return asImp().boundaryType( intersection );
    }

    /** \brief evaluate Dirichlet boundary value
     *
     *  Given an intersection and a quadrature-point pair, evaluate the
     *  Dirichlet boundary data in the specified point.
     *
     *  \note The entity the boundary belongs to can be obtained from the
     *        intersection by the intersecion.inside method.
     *
     *  \note The quadrature is not guaranteed to be a codim-1-quadrature
     *
     *  \param[in]  intersection intersection on which the quadrature point lies
     *
     *  \param[in]  quadrature   quadrature to use
     *
     *  \param[in]  point        index of the evaluation point within the
     *                           quadrature
     *
     *  \param[out] phi          Dirichlet value in the quadrature point
     */
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues ( const IntersectionIteratorType &intersection,
                                  const QuadratureType &quadrature,
                                  int point,
                                  RangeType &phi ) const
    {
      assert( Properties :: hasDirichletValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().dirichletValues( intersection, quadrature, point, phi ) );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues ( const IntersectionIteratorType &intersection,
                                const QuadratureType &quadrature,
                                int point,
                                RangeType &phi ) const
    {
      assert( Properties :: hasNeumannValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().neumannValues( intersection, quadrature, point, phi ) );
    }
    
    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues ( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int point,
                              RangeType &phi ) const
    {
      assert( Properties :: hasRobinValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().robinValues( intersection, quadrature, point, phi ) );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType robinAlpha ( const IntersectionIteratorType &intersection,
                                       const QuadratureType &quadrature,
                                       int point ) const
    {
      assert( Properties :: hasRobinValues );
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().robinAlpha( intersection, quadrature, point ) );
      return asImp().robinAlpha( intersection, quadrature, point );
    }

  protected:
    inline const BoundaryModelImp &asImp () const
    {
      return static_cast< const BoundaryModelImp& >( *this );
    }

    inline BoundaryModelImp &asImp ()
    {
      return static_cast< BoundaryModelImp& >( *this );
    }
  };



  template< class FunctionSpaceImp,
            class BoundaryModelImp,
            class PropertiesImp = DefaultBoundaryModelProperties >
  class BoundaryModelDefault
  : public BoundaryModelInterface< FunctionSpaceImp, BoundaryModelImp, PropertiesImp >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef PropertiesImp Properties;

  private:
    typedef BoundaryModelDefault< FunctionSpaceType, BoundaryModelImp, Properties >
      ThisType;
    typedef BoundaryModelInterface< FunctionSpaceType, BoundaryModelImp, Properties >
      BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    
  public:
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues ( const IntersectionIteratorType &intersection,
                                  const QuadratureType &quadrature,
                                  int point,
                                  RangeType &phi ) const
    {
      assert( Properties :: hasDirichletValues );
      phi = 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues ( const IntersectionIteratorType &intersection,
                                const QuadratureType &quadrature,
                                int point,
                                RangeType &phi ) const
    {
      assert( Properties :: hasNeumannValues );
      phi = 0;
    }
    
    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues ( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int point,
                              RangeType &phi ) const
    {
      assert( Properties :: hasRobinValues );
      phi = 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType robinAlpha ( IntersectionIteratorType &intersection,
                                       QuadratureType &quadrature,
                                       int point ) const
    {
      assert( Properties :: hasRobinValues );
      return 0;
    }
  };

}

#endif
