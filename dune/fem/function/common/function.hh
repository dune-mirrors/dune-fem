#ifndef DUNE_FEM_FUNCTION_HH
#define DUNE_FEM_FUNCTION_HH

// dune-common includes
#include <dune/common/fvector.hh>

// dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/operator/common/mapping.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    /** @addtogroup Functions
        Functions are mapping from one finite dimensional
        vector space into another, e.g.,
        \f$K^n\f$ into \f$L^m\f$.
        They are element of a
        \ref FunctionSpaceInterface "function space".

        \remark The interface is given by Function.
        @{
    **/


    /*! \brief
        Abstract class representing a function


        Template parameters are:
        -  FunctionSpaceImp      type of the function space where the function
                                 belongs to.
        -  FunctionImp           type of the implemented function (Barton-Nackman)

        @interfaceclass
    **/
    template< class FunctionSpaceImp, class FunctionImp >
    class Function
    : public BartonNackmanInterface< Function< FunctionSpaceImp, FunctionImp >,
                                     FunctionImp >,
      public Mapping < typename FunctionSpaceImp :: DomainFieldType,
                       typename FunctionSpaceImp :: RangeFieldType,
                       typename FunctionSpaceImp :: DomainType,
                       typename FunctionSpaceImp :: RangeType >
    {
      typedef Function< FunctionSpaceImp, FunctionImp > ThisType;
      typedef BartonNackmanInterface< ThisType, FunctionImp > BaseType;

    public:
      //! type of function space this function belongs to
      typedef FunctionSpaceImp FunctionSpaceType;

      //! type of the implementation (Barton-Nackman)
      typedef FunctionImp FunctionType;

      //! field type of domain
      typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
      //! field type of range
      typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
      //! domain type
      typedef typename FunctionSpaceType :: DomainType DomainType;
      //! range type
      typedef typename FunctionSpaceType :: RangeType RangeType;
      //! jacobian type
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
      //! hessian type
      typedef typename FunctionSpaceType :: HessianRangeType HessianRangeType;

      //! type of mapping base class
      typedef Mapping< DomainFieldType, RangeFieldType, DomainType, RangeType >
        MappingType;

    protected:
      using BaseType::asImp;

      /** \brief default constructor */
      Function () = default;

      Function ( const ThisType& ) = default;

    public:
      ThisType& operator= ( const ThisType& ) = delete;

      //! destructor
      virtual ~Function ()
      {}

      /** \brief application operator call evaluate
          \param[in] arg argument
          \param[out] dest destination, i.e. f(arg)
      */
      virtual void operator()(const DomainType & arg, RangeType & dest) const
      {
        evaluate(arg,dest);
      }

      /** \brief evaluate the function
       *
       *  \param[in]  x      evaluation point
       *  \param[out] value  value of the function in x
       */
      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().evaluate( x, value ) );
      }

      /** \brief evaluate the Jacobian of the function
       *
       *  \param[in]  x         evaluation point
       *  \param[out] jacobian  value of the Jacobian in x
       */
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().jacobian( x, jacobian ) );
      }

      /** \brief evaluate the hessian of the function
       *
       *  \param[in]  x        evaluation point
       *  \param[out] hessian  value of the hessian in x
       */
      void hessian ( const DomainType &x, HessianRangeType &hessian ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().hessian( x, hessian ) );
      }

    private:
      //! Helper function for Mapping
      //! With this function, a combined mapping can choose the right application
      //! operator (i.e. the one from Mapping itself, or from Function/Operator)
      //! \note: Do not override this definition
      virtual void apply (const DomainType& arg, RangeType& dest) const
      {
        operator()(arg, dest);
      }
    };

    ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HH
