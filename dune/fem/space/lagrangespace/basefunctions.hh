#ifndef DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH
#define DUNE_LAGRANGESPACE_BASEFUNCTIONS_HH

#include "genericbasefunctions.hh"

namespace Dune
{

  // LagrangeBaseFunction
  // --------------------
   
  template< class FunctionSpace, GeometryType :: BasicType type,
            unsigned int dim, unsigned int pOrder >
  class LagrangeBaseFunction
  : public BaseFunctionInterface< FunctionSpace >
  {
    typedef LagrangeBaseFunction< FunctionSpace, type, dim, pOrder > ThisType;
    typedef BaseFunctionInterface< FunctionSpace > BaseType;

  public:
    typedef typename GeometryWrapper< type, dim > :: GenericGeometryType
      GenericGeometryType;
    typedef GenericLagrangeBaseFunction
      < FunctionSpace, GenericGeometryType, pOrder >
      GenericBaseFunctionType;
      
    typedef typename GenericBaseFunctionType :: DomainType DomainType;
    typedef typename GenericBaseFunctionType :: RangeType RangeType;

  protected:
    GenericBaseFunctionType baseFunction_;

  public:
    LagrangeBaseFunction( unsigned int baseNum );

    virtual void evaluate ( const FieldVector< deriType, 0 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;
    
    virtual void evaluate ( const FieldVector< deriType, 1 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;

    virtual void evaluate ( const FieldVector< deriType, 2 > &diffVariable,
                            const DomainType &x,
                            RangeType &phi ) const;

    virtual int order () const;
  };



#ifndef COMPILE_FEMLIB
  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 1, 1 >, type, 1, 1 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 2, 1 >, type, 2, 1 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 3, 1 >, type, 3, 1 >;


  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 1, 1 >, type, 1, 1 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 2, 1 >, type, 2, 1 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 3, 1 >, type, 3, 1 >;


  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 1, 1 >, type, 1, 2 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 2, 1 >, type, 2, 2 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 3, 1 >, type, 3, 2 >;


  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 1, 1 >, type, 1, 2 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 2, 1 >, type, 2, 2 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 3, 1 >, type, 3, 2 >;


  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 1, 1 >, type, 1, 3 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 2, 1 >, type, 2, 3 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< double, double, 3, 1 >, type, 3, 3 >;


  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 1, 1 >, type, 1, 3 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 2, 1 >, type, 2, 3 >;

  template< GeometryType :: BasicType type >
  class LagrangeBaseFunction< FunctionSpace< float, float, 3, 1 >, type, 3, 3 >;
#endif



  // LagrangeBaseFunctionFactory
  // ---------------------------

  template< class ScalarFunctionSpace, unsigned int dim, unsigned int pOrder >
  class LagrangeBaseFunctionFactory
  : public BaseFunctionFactory< ScalarFunctionSpace >
  {
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpace, dim, pOrder >
      ThisType;
    typedef BaseFunctionFactory< ScalarFunctionSpace > BaseType;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template< class ScalarFunctionSpace, unsigned int pOrder >
  class LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
  : public BaseFunctionFactory< ScalarFunctionSpace >
  {
    typedef LagrangeBaseFunctionFactory< ScalarFunctionSpace, 3, pOrder >
      ThisType;
    typedef BaseFunctionFactory< ScalarFunctionSpace > BaseType;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry );

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



#ifndef COMPILE_FEMLIB
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > >
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > >
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > >
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 1 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 1, 1 > >
  {
    typedef FunctionSpace< float, float, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 1 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 2, 1 > >
  {
    typedef FunctionSpace< float, float, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 1 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 3, 1 > >
  {
    typedef FunctionSpace< float, float, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > >
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > >
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > >
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 2 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 1, 1 > >
  {
    typedef FunctionSpace< float, float, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 2 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 2, 1 > >
  {
    typedef FunctionSpace< float, float, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 2 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 3, 1 > >
  {
    typedef FunctionSpace< float, float, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > >
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > >
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > >
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };



  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1, 3 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 1, 1 > >
  {
    typedef FunctionSpace< float, float, 1, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2, 3 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 2, 1 > >
  {
    typedef FunctionSpace< float, float, 2, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class LagrangeBaseFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3, 3 >
  : public BaseFunctionFactory< FunctionSpace< float, float, 3, 1 > >
  {
    typedef FunctionSpace< float, float, 3, 1 > ScalarFunctionSpace;

  public:
    explicit LagrangeBaseFunctionFactory ( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual ~LagrangeBaseFunctionFactory ();

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };
#endif
  
}

#include "basefunctions_inline.hh"

#endif
