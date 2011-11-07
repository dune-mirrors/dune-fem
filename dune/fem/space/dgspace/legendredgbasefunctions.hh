#ifndef DUNE_LEGENDREDGBASEFUNCTIONS_HH
#define DUNE_LEGENDREDGBASEFUNCTIONS_HH

#include <dune/common/misc.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/array.hh>

#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>

#include <dune/fem/space/dgspace/legendrepoly.hh>

namespace Dune
{
  
  namespace Fem {  


  //! number of legendre base functions for given polord and dim 
  template <int p, int dim>
  struct NumLegendreBaseFunctions
  {
    enum { numBaseFct = Power_m_p<p+1,dim>::power };
  };

  template < class FunctionSpaceType, int polOrd>
  class LegendreDGBaseFunction :
    public BaseFunctionInterface<FunctionSpaceType>
  {
    //Template Meta Programm for evaluating tensorproduct polynomial in arbitrary dimensions
    template<int dim,int i,int PolOrd>
    struct Eval
    {
      template <class FieldVectorType>
      static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int idx)
      {
        assert( int(x.dimension) == int(dim) );
        const int num = idx % (PolOrd+1);
        // LegendrePoly lp=LegendrePoly(num);
        return lp.evaluate(num,x[i-1]) * Eval<dim, i-1, PolOrd>::apply(lp,x,(idx-num)/(PolOrd+1));
      }
    };

    template<int dim,int PolOrd>
    struct Eval<dim,0,PolOrd>
    {
      template <class FieldVectorType>
      static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int idx)
      {
        return 1.0;
      }
    };

    //Template Meta Programm for evaluating the partial derivative of tensorproduct polynomials in arbitrary dimensions
    template<int dim,int i,int PolOrd>
    struct EvalD
    {
      template <class FieldVectorType>
      static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int j, const int idx)
      {
        assert( int(x.dimension) == int(dim) );
        const int num=idx%(PolOrd+1);
        if( (i-1) != j )
          return lp.evaluate(num,x[i-1]) * EvalD<dim,i-1,PolOrd>::apply(lp,x, j, (idx-num)/(PolOrd+1));
        else
          return lp.jacobian(num,x[i-1]) * EvalD<dim,i-1,PolOrd>::apply(lp,x, j, (idx-num)/(PolOrd+1));
      }
    };

    template<int dim,int PolOrd>
    struct EvalD<dim,0,PolOrd>
    {
      template <class FieldVectorType>
      static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int j, const int idx)
      {
        return 1.0;
      }
    };

  protected:
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    enum{ dim=DomainType::dimension };

  public:
    LegendreDGBaseFunction( const int baseNum ) 
      : lp_( LegendrePoly :: instance() ),
        baseNum_( baseNum )
    {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }
    
    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = Eval<dim,dim,polOrd>::apply( lp_, x, baseNum_ );
    }

    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = EvalD<dim,dim,polOrd>::apply( lp_, x, diffVariable[0], baseNum_);
    
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      assert(false); // Not implemented
      abort();
    }

    static int numBaseFunctions() 
    {
      return NumLegendreBaseFunctions<polOrd,dim>::numBaseFct;
    }

  protected:
    // reference to legendre polynomials 
    const LegendrePoly& lp_;

    // my number in the set of basis functions
    const int baseNum_;
  };

  template <class ScalarFunctionSpaceImp, int polOrd>
  class LegendreDGBaseFunctionFactory : 
    public BaseFunctionFactory<ScalarFunctionSpaceImp> 
  {
  public:
    typedef ScalarFunctionSpaceImp FunctionSpaceType;
    dune_static_assert(FunctionSpaceType::dimRange == 1,"only scalar functions spaces allowed" ) ;
    enum { dim = FunctionSpaceType::dimDomain };

    // true if hierarchical ordering of basis functions 
    // should be used 
    static const bool hierarchical = true ; 

    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;

    // number of basis functions 
    enum { numBaseFct = NumLegendreBaseFunctions<polOrd,dim>::numBaseFct };

  protected:  
    typedef array<int, 2*numBaseFct> BaseNumberMapType;

    template<int pOrd, int d>
    struct DimLoop 
    {
      template < class A >
      static void loop( const int term, 
                        A& polOrders,
                        BaseNumberMapType& baseMap, 
                        int& oldBaseNum,
                        int& newBaseNum )
      {
        // go through all p because of old base counter
        for( int p=0; p<=pOrd; ++p) 
        {
          polOrders[ d ] = p ;
          DimLoop< pOrd, d-1> :: loop( term, polOrders, baseMap, oldBaseNum, newBaseNum );
        }
      }
    };

    template<int pOrd>
    struct DimLoop<pOrd, 0>
    {
      template< class A >
      static bool checkEntries( const A& array, const int p ) 
      {
        bool found = false ;
        for( int i=0; i<array.size(); ++i) 
        {
          if( array[ i ] == p ) found = true ;
          if( array[ i ] > p ) return false ;
        }
        return found ;
      }

      template< class A > 
      static void loop( const int term, 
                        A& polOrders,
                        BaseNumberMapType& baseMap, 
                        int& oldBaseNum,
                        int& newBaseNum )
      {
        // go through all p because of old base counter
        for( int p=0; p<=pOrd; ++p, ++oldBaseNum) 
        {
          polOrders[ 0 ] = p ;

          // at least one of the entries has to be of order p
          // non of the entries must be larger than p 
          if( checkEntries( polOrders, term ) )
          {
            assert( newBaseNum < int(numBaseFct) );
            baseMap[ newBaseNum ] = oldBaseNum;
            ++newBaseNum;
          }
        }
      }
    };

    // map for renumbering of basis functions to have a hierarchical basis 
    BaseNumberMapType baseFunctionMap_;
  public:
    //! constructor creating basis function factory 
    LegendreDGBaseFunctionFactory(const GeometryType& geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>( geo ),
      baseFunctionMap_()
    {
      baseFunctionMap_.fill( -1 );

      // switch numbering 
      if( hierarchical ) 
      {

        // small array to store the curently used pol ords 
        array<int, dim> polStorage ;
        polStorage.fill( -1 );

        int newBaseNum = 0;
        // check for all terms the number of base functions 
        for(int term=0; term <= polOrd; ++term )
        {
          // initialize oldBaseNum counter new for each term 
          int oldBaseNum = 0;
          // construct mapping for this term 
          DimLoop<polOrd, dim-1> :: 
            loop( term, polStorage, baseFunctionMap_, oldBaseNum, newBaseNum );
        }
      }

#ifndef NDEBUG 
      //std::cout << numBaseFct << std::endl;
      for( int i=0; i<numBaseFct; ++i ) 
      {
        //std::cout << "bsae " << i << " is " << baseFunctionMap_[ i ] << std::endl;
        assert( baseFunctionMap_[ i ] >= 0 );
      }
#endif
    }

    virtual BaseFunctionType* baseFunction(int i) const 
    {
      // only for cubes we have LegendreBaseFunctions 
      if( ! this->geometry().isCube() )
      {
        DUNE_THROW(NotImplemented,"LegendreBaseFunctions only implemented for cubes!");
      }

      // if hierarchical is enabled then switch numbering 
      const int baseNum = ( hierarchical ) ? baseFunctionMap_[ i ] : i;
      return new LegendreDGBaseFunction<FunctionSpaceType, polOrd> ( baseNum );
    }
    
    virtual int numBaseFunctions() const 
    {
      return numBaseFct;
    }
  };

  } // end namespace Fem 

} // end namespace Dune
#endif
