#ifndef DUNE_BASEFUNCTIONPROXY_HH
#define DUNE_BASEFUNCTIONPROXY_HH

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

namespace Dune
{
  namespace Fem
  {

  /** \addtogroup BaseFunction
   *  \{
   */
    
  template < class BaseFunctionSetImp > 
  class SimpleBaseFunctionProxy; 

  template <class BaseFunctionSetImp>
  struct SimpleBaseFunctionProxyTraits  
  {
    //! Export function space type
    typedef typename BaseFunctionSetImp :: Traits :: FunctionSpaceType FunctionSpaceType;

    //! Exact type of the base function
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
  };

  /** \class SimpleBaseFunctionProxy
   *  \brief proxy object for base function sets
   */
  template< class BaseFunctionSetImp > 
  class SimpleBaseFunctionProxy
  : public BaseFunctionSetDefault< SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > >
  {
  public:
    typedef SimpleBaseFunctionProxyTraits< BaseFunctionSetImp > Traits;
    
  private:
    typedef SimpleBaseFunctionProxy< BaseFunctionSetImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;

  public:
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

    enum { dimRange = FunctionSpaceType :: dimRange };
    enum { dimDomain = FunctionSpaceType :: dimDomain };

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  public:
    SimpleBaseFunctionProxy ()
    : baseSet_( 0 )
    {}
      
    SimpleBaseFunctionProxy ( const BaseFunctionSetImp* baseSet )
    : baseSet_( baseSet ) 
    {}

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::size */
    size_t size () const
    {
      return baseFunctionSet().size();
    }
    
    /** \copydoc Dune::Fem::BaseFunctionSetInterface::geometryType */
    GeometryType geometryType () const
    {
      return baseFunctionSet().geometryType();
    }
   
    /** \copydoc Dune::Fem::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrder> &diffVariable,const Point &x,RangeType &value) const */
    template< int diffOrder, class Point >
    DUNE_VERSION_DEPRECATED(1,4,remove)
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrder > &diffVariable,
                    const Point &x,
                    RangeType &value ) const
    {
      baseFunctionSet().evaluate( baseFunction, diffVariable, x, value );
    }

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::evaluate(const int baseFunction,const Point &x,RangeType &value) const */
    template< class Point >
    DUNE_VERSION_DEPRECATED(1,4,remove)
    void evaluate ( const int baseFunction,
                    const Point &x,
                    RangeType &value ) const
    {
      baseFunctionSet().evaluate( baseFunction, x, value );
    }

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::jacobian(const int baseFunction,const Point &x,JacobianRangeType &jacobian) const */
    template< class Point >
    DUNE_VERSION_DEPRECATED(1,4,remove)
    void jacobian( const int baseFunction,
                   const Point &x,
                   JacobianRangeType &jacobian ) const
    {
      baseFunctionSet().jacobian( baseFunction, x, jacobian );
    }

    template< class PointType, 
              class LocalDofVectorType >
    void evaluateAll ( const PointType &x,
                       const LocalDofVectorType& dofs,
                       RangeType& ret) const
    {
      baseFunctionSet().evaluateAll( x, dofs, ret );
    }

    template< int diffOrder, 
              class PointType, 
              class LocalDofVectorType >
    void evaluateAll ( const FieldVector<int, diffOrder>& diffVariable,
                       const PointType &x,
                       const LocalDofVectorType& dofs,
                       RangeType& ret) const
    {
      baseFunctionSet().evaluateAll( diffVariable, x, dofs, ret );
    }

    template< class PointType, class RangeVectorType >
    void evaluateAll ( const PointType &x, RangeVectorType &ret ) const
    {
      baseFunctionSet().evaluateAll( x, ret );
    }


    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit,   
                       const LocalDofVectorType& dofs,
                       GlobalJacobianRangeType& ret) const
    {
      baseFunctionSet().jacobianAll( x, gjit, dofs, ret );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeVectorType >
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       GlobalJacobianRangeVectorType &ret ) const
    {
      baseFunctionSet().jacobianAll( x, gjit, ret );
    }

    template< class PointType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const RangeType& factor,
                LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, factor, dofs );
    }

    template< class QuadratureType, 
              class LocalDofVectorType,
              class RangeVectorType>
    void evaluateRanges ( const QuadratureType& quad,
                          const LocalDofVectorType& dofs,
                          RangeVectorType& rangeVector) const
    {
      baseFunctionSet().evaluateRanges( quad, dofs, rangeVector );
    }

    template< class QuadratureType, 
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    void evaluateJacobians( const QuadratureType& quad,
                            const Geometry& geometry, 
                            const LocalDofVectorType& dofs,
                            JacobianRangeVectorType& jacVector) const
    {
      baseFunctionSet().evaluateJacobians( quad, geometry, dofs, jacVector );
    }

    template< class QuadratureType, 
              class RangeVectorType,
              class LocalDofVectorType >
    void axpyRanges ( const QuadratureType& quad,
                      const RangeVectorType& factorVec,
                      LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpyRanges( quad, factorVec, dofs );
    }

    template< class QuadratureType, 
              class Geometry,
              class JacobianRangeVectorType,
              class LocalDofVectorType >
    void axpyJacobians( const QuadratureType& quad,
                        const Geometry& geometry,
                        const JacobianRangeVectorType& factorVec,
                        LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpyJacobians( quad, geometry, factorVec, dofs );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit,   
                const GlobalJacobianRangeType& factor,
                LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, gjit, factor, dofs );
    }

    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit,
                const RangeType& factor1,
                const GlobalJacobianRangeType& factor2,
                LocalDofVectorType& dofs) const
    {
      baseFunctionSet().axpy( x, gjit, factor1, factor2, dofs );
    }

  protected:
    const BaseFunctionSetImp &baseFunctionSet () const
    {
      assert( baseSet_ );
      return *baseSet_;
    }

  protected:
    // base function set 
    const BaseFunctionSetImp *baseSet_; 
  }; // end SimpleBaseFunctionProxy 

  /** \} */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_BASEFUNCTIONPROXY_HH
