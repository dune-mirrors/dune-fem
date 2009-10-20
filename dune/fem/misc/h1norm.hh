#ifndef DUNE_FEM_H1NORM_HH
#define DUNE_FEM_H1NORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

  template< class GridPart >
  class H1Norm
  : public L2Norm< GridPart >
  {
  public:
    typedef GridPart GridPartType;

  private:
    typedef H1Norm< GridPartType > ThisType;
    typedef L2Norm< GridPartType > BaseType;

  protected:
    template< class Function >
    class FunctionJacobianSquare;

  protected:
    typedef typename BaseType :: GridIteratorType GridIteratorType;
    typedef typename BaseType :: IntegratorType IntegratorType;

    typedef typename GridIteratorType :: Entity EntityType;

  protected:
    using BaseType :: gridPart;
    using BaseType :: comm;

  public:
    inline explicit H1Norm ( const GridPartType &gridPart );
    inline H1Norm ( const ThisType &other );

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );

  public:
    template< class DiscreteFunctionType >
    inline typename DiscreteFunctionType :: RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;
    
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename UDiscreteFunctionType :: RangeFieldType
    distance ( const UDiscreteFunctionType &u,
               const VDiscreteFunctionType &v ) const;
  };

}

#include "h1norm_inline.hh"

#endif
