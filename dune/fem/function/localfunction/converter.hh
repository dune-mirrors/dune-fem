#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONVERTER_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CONVERTER_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/fem/function/common/function.hh>
#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/space/common/functionspace.hh>


namespace Dune
{

  namespace Fem
  {

    // LocalFunctionConverter
    // --------------------

    /** \brief implementation of a Dune::Fem::LocalFunction on a FunctionSpace V restircted/prolongated from an other
     *         function space W.
     *
     *  The HostLocalFunction is assumed to fulfill the LocalFunctioninterface.
     *  Basically the following functions are implemented on the HostLocalFunction:
     *  \code
     *    template< class Point >
     *    void evaluate ( const Point &x, RangeType &ret ) const;
     *
     *    template< class Point >
     *    void jacobian ( const Point &x, JacobianRangeType &jac ) const;
     *
     *    tempalte< class Point >
     *    void hessian ( const Point &x, HessianRangeType &hess ) const;
     *
     *    const EntityType &entity () const;
     *
     *    int size() const;
     *
     *    void init ( const EntityType & entity );
     *  \endcode
     *
     *  The template paramter Converter, is used to get the restriction/prolongation onto the space V.
     *  Converter is expected to provide the method
     *  {Hessian,Jacobian,.}RangeType converter( Host{Hessian,Jacobian,. }RangeType );
     *  which does the acctual mapping onto V.
     *  The dimension of the new Range is obtained from the method Convertor::operator( HostRangeType )::dimension.
     *
     *  Users may prescribe how the parameter localFunction is stored by providing a
     *  fourth template parameter, the storage policy. Further informations on the storage policy can be found in
     *  the file dune/fem/function/common/instationary.hh.
     *
     *  The free-standing function
     *  \code
     *    Dune::Fem::localFunctionConverter
     *  \endcode
     *  may be used to conveniently create a new instance of a LocalFunctionConverter. Use
     *  \code
     *    auto g = localFunctionConverter( localFunction, Converter )
     *  \endcode
     *  to create an converted local function.
     *
     *  \tparam  HostLocalFunction original local function
     *  \tparam  Converter structure which provides the restirction/prolongation from W to V
     *  \tparam  StoragePolicy  storage policy
     *
     *  \ingroup LocalFunction
     */

    template< class HostLocalFunction, class Converter, template< class > class Storage = __InstationaryFunction::HoldCopy >
    class LocalFunctionConverter
      : private Storage< HostLocalFunction >
    {
      typedef LocalFunctionConverter< HostLocalFunction, Converter, Storage > ThisType;
      typedef Storage< HostLocalFunction > BaseType;

      typedef typename HostLocalFunction::RangeType HostRangeType;
      typedef typename HostLocalFunction::JacobianRangeType HostJacobianRangeType;
      typedef typename HostLocalFunction::HessianRangeType HostHessianRangeType;

    public:
      // obtain new dimRange from Converter
      static const int dimRange = decltype( std::declval< Converter >() ( std::declval< HostRangeType >() ) ) ::dimension;

      // define new FunctionSpace
      typedef typename ToNewDimRangeFunctionSpace< typename HostLocalFunction::FunctionSpaceType, dimRange >::Type FunctionSpaceType;

      // types from HostLocalFunction
      typedef typename HostLocalFunction::EntityType EntityType;

      // types from FunctionSpace
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      static const int dimDomain = FunctionSpaceType::dimDomain;

      struct Traits
      {
        typedef typename FunctionSpaceType::DomainType DomainType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
        typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;
        typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      };

      LocalFunctionConverter ( const HostLocalFunction &hostLocalFunction, const Converter &converter = Converter() )
        : BaseType( hostLocalFunction ), converter_( converter )
      {}

      LocalFunctionConverter ( HostLocalFunction &&hostLocalFunction, const Converter &converter = Converter() )
        : BaseType( std::move( hostLocalFunction ) ), converter_( converter )
      {}

      template< class Point >
      void evaluate ( const Point &p, RangeType &ret ) const
      {
        HostRangeType hRet;
        this->get().evaluate( p, hRet );
        ret = converter_( hRet );
      }
      template< class Point >
      RangeType operator()( const Point &p ) const
      {
        RangeType ret;
        evaluate(p,ret);
        return ret;
      }

      template< class Point >
      void jacobian ( const Point &p, JacobianRangeType &jac ) const
      {
        HostJacobianRangeType hJac;
        this->get().jacobian( p, hJac );
        jac = converter_( hJac );
      }

      template< class Point >
      void hessian ( const Point &p, HessianRangeType &hes ) const
      {
        HostHessianRangeType hHes;
        this->get().hessian( p, hHes );
        hes = converter_( hHes );
      }

      template< class Quadrature, class ... Vectors >
      void evaluateQuadrature ( const Quadrature &quad, Vectors& ... vector ) const
      {
        std::ignore = std::make_tuple(
            ( evaluateQuadratureImp( quad, vector, vector[ 0 ] ), 1 )... );
      }

      int order () const { return this->get().order(); }

      EntityType entity () const { return this->get().entity(); }

      void init ( const EntityType &entity ) { this->get().init( entity ); }

    protected:
      template< class QuadratureType, class VectorType >
      void evaluateQuadratureImp ( const QuadratureType &quadrature, VectorType &values, const RangeType & ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluate( quadrature[ qp ], values[ qp ] );
      }

      template< class QuadratureType, class VectorType >
      void evaluateQuadratureImp ( const QuadratureType &quadrature, VectorType &values, const JacobianRangeType & ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobian( quadrature[ qp ], values[ qp ] );
      }

      template< class QuadratureType, class VectorType >
      void evaluateQuadratureImp ( const QuadratureType &quadrature, VectorType &values, const HessianRangeType & ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          hessian( quadrature[ qp ], values[ qp ] );
      }

      Converter converter_;
    };



    // localFunctionConverter
    // ----------------------

    template< class HostLocalFunction, class Converter >
    LocalFunctionConverter< HostLocalFunction, Converter, __InstationaryFunction::HoldCopy >
    localFunctionConverter ( HostLocalFunction hostLocalFunction, const Converter &converter = Converter() )
    {
      typedef LocalFunctionConverter< HostLocalFunction, Converter, __InstationaryFunction::HoldCopy > LocalFunctionConverterType;
      return LocalFunctionConverterType( std::move( hostLocalFunction ), converter );
    }

    template< class HostLocalFunction, class Converter >
    LocalFunctionConverter< typename std::remove_const< HostLocalFunction >::type, Converter, __InstationaryFunction::HoldReference >
    localFunctionConverter ( std::reference_wrapper< HostLocalFunction > hostLocalFunction, const Converter &converter = Converter() )
    {
      typedef LocalFunctionConverter< typename std::remove_const< HostLocalFunction >::type, Converter, __InstationaryFunction::HoldReference > LocalFunctionConverterType;
      return LocalFunctionConverterType( hostLocalFunction.get(), converter );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONVERTER_HH
