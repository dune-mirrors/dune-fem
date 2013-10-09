#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_DEFAULT_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_DEFAULT_HH

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/misc/engineconcept.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    /** \addtogroup LocalFunction
     * \{
     */

    // LocalFunctionDefault
    // --------------------

    /** \class LocalFunctionDefault
     *  \brief default implementation of a LocalFunction
     */
    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    class LocalFunctionDefault
    : public EngineDefault< LocalFunctionImp >
    {
      typedef EngineDefault< LocalFunctionImp > BaseType;

    public:
      //! type of  discrete function space the local function belongs to
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      //! type of the entity, the local function lives on is given by the space 
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

      //! field type of the domain
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      //! field type of the range
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
      //! type of domain vectors, i.e., type of coordinates
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      //! type of range vectors, i.e., type of function values
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      //! type of Jacobian, i.e., type of evaluated Jacobian matrix
      typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! type of the Hessian
      typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

      //! type of base function set  
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType
        BasisFunctionSetType;

      //! dimension of the domain
      enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
      //! dimension of the range
      enum { dimRange = DiscreteFunctionSpaceType::dimRange };

    protected:
      //! type of entity's geometry
      typedef typename EntityType::Geometry GeometryType;

    public:
      /** \copydoc Dune::Fem::LocalFunction::operator+=(const LocalFunction<T> &lf) */
      template< class T >
      void operator+= ( const LocalFunction< T > &lf );

      /** \copydoc Dune::Fem::LocalFunction::operator-=(const LocalFunction<T> &lf) */
      template< class T >
      void operator-= ( const LocalFunction< T > &lf );

      template< class T >
      void assign ( const LocalFunction< T > &lf );

      void clear ( );
      
      /** \copydoc Dune::Fem::LocalFunction::axpy(const RangeFieldType s,const LocalFunction<T> &lf) */
      template< class T >
      void axpy ( const RangeFieldType s, const LocalFunction< T > &lf );

      /** \copydoc Dune::Fem::LocalFunction::evaluate(const PointType &x,RangeType &ret) const */
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const;

      /** \copydoc Dune::Fem::LocalFunction::jacobian(const PointType &x,JacobianRangeType &ret) const */
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const;

      template< class PointType >
      void hessian ( const PointType &x, HessianRangeType &hessian ) const;

      /** \copydoc Dune::Fem::LocalFunction::axpy(const PointType &x,const RangeType &factor) */
      template< class PointType >
      void axpy( const PointType &x, const RangeType &factor );

      /** \copydoc Dune::Fem::LocalFunction::axpy(const PointType &x,const JacobianRangeType &factor) */
      template< class PointType >
      void axpy( const PointType &x, const JacobianRangeType &factor );

      /** \copydoc Dune::Fem::LocalFunction::axpy(const PointType &x,const RangeType &factor1,const JacobianRangeType &factor2) */
      template< class PointType >
      void axpy( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 );

      /** \copydoc Dune::Fem::LocalFunction::order() const */
      int order () const { return asImp().basisFunctionSet().order(); }

      int numScalarDofs () const
      {
        const int numDofs = asImp().numDofs();
        assert( numDofs % dimRange == 0 );
        return numDofs / dimRange;
      }

    private:
      template< class PointType >
      void doEvaluate ( const FieldVector< int, 0 > &,
                        const PointType &x, RangeType &ret ) const
      {
        return evaluate( x, ret );
      }

      template< class PointType >
      void doEvaluate ( const FieldVector< int, 1 > &diffVariable,
                        const PointType &x, RangeType &ret ) const
      {
        JacobianRangeType tmp;
        jacobian( x, tmp );
        const int j = diffVariable[ 0 ];
        for( int i = 0; i < dimRange; ++i )
          ret[ i ] = tmp[ i ][ j ];
      }

      template< int diffOrder, class PointType >
      void doEvaluate ( const FieldVector< int, diffOrder > &diffVariable,
                        const PointType &x, RangeType &ret ) const
      {
        DUNE_THROW( NotImplemented, "Method evaluate() not implemented yet." );
      }

    protected:
      using BaseType::asImp;
    };

    /** \} */



    // Implementation of LocalFunctionDefault
    // --------------------------------------

    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class T >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::operator+= ( const LocalFunction< T > &lf )
    {
      asImp().axpy( RangeFieldType(1), lf );
    }



    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class T >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::operator-= ( const LocalFunction< T > &lf )
    {
      asImp().axpy( RangeFieldType(-1), lf );
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class T >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::assign ( const LocalFunction< T > &lf )
    {
      const int numDofs = asImp().numDofs();
      assert( numDofs == lf.numDofs() );

      for( int i = 0; i < numDofs; ++i )
        asImp()[ i ] = lf[ i ];
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::clear ()
    {
      const int numDofs = asImp().numDofs();
      for( int i = 0; i < numDofs; ++i )
        asImp()[ i ] = 0.0;
    }
    
    
    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class T >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::axpy ( const RangeFieldType s, const LocalFunction< T > &lf )
    {
      const int numDofs = asImp().numDofs();
      assert( numDofs == lf.numDofs() );

      for( int i = 0; i < numDofs; ++i )
        asImp()[ i ] += s * lf[ i ];
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< int diffOrder, class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                   const PointType &x, RangeType &ret ) const
    {
      // asImp().baseFunctionSet().evaluateAll( diffVariable, x, asImp(), ret );
      doEvaluate( diffVariable, x, ret );
    }

    
    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      :: evaluate ( const PointType &x,
                    RangeType &ret ) const
    {
      asImp().basisFunctionSet().evaluateAll( x, asImp(), ret );
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::jacobian ( const PointType &x, JacobianRangeType &ret ) const
    {
      asImp().basisFunctionSet().jacobianAll( x, asImp(), ret );
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::hessian ( const PointType &x, HessianRangeType &hessian ) const
    {
      asImp().basisFunctionSet().hessianAll( x, asImp(), hessian );
    }
   
    
    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::axpy ( const PointType &x, const RangeType &factor )
    {
      asImp().basisFunctionSet().axpy( x, factor, asImp() );
    }


    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::axpy ( const PointType &x, const JacobianRangeType &factor )
    {
      asImp().basisFunctionSet().axpy( x, factor, asImp() );
    }

    
    template< class DiscreteFunctionSpace, class LocalFunctionImp >
    template< class PointType >
    inline void LocalFunctionDefault< DiscreteFunctionSpace, LocalFunctionImp >
      ::axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
    {
      asImp().basisFunctionSet().axpy( x, factor1, factor2, asImp() );
    }

  } // end namespace Fem

} // end namespace Dune 

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_DEFAULT_HH
