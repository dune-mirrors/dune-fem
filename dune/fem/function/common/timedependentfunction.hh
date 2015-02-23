#ifndef DUNE_FEM_TIMEDEPENDENTFUNCTION_HH
#define DUNE_FEM_TIMEDEPENDENTFUNCTION_HH

// dune-fem includes
#include <dune/fem/function/common/function.hh>

namespace Dune
{

  namespace Fem
  {
    /** @addtogroup Functions
        \remark The interface is given by Function.
        @{
    **/

    /*! \brief
        Abstract class representing a function f(t,.)

        Template parameters are:
        -  F   type of the time dependent function that we want turn into a function by fixing time

        @interfaceclass
    **/
    template <class F>
    class TimeDependentFunction
      : public Function< typename F :: FunctionSpaceType,
                         TimeDependentFunction< F > >
    {
      typedef F FunctionType ;
    public:
      typedef typename FunctionType :: FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType :: DomainType                   DomainType;
      typedef typename FunctionSpaceType :: RangeType                    RangeType;
      typedef typename FunctionSpaceType :: DomainFieldType              DomainFieldType;
      typedef typename FunctionSpaceType :: RangeFieldType               RangeFieldType;
      typedef typename FunctionSpaceType :: JacobianRangeType            JacobianRangeType;

    public:
      TimeDependentFunction( const FunctionType& f, const double time )
        : function_( f ), time_( time )
      {}

      //! forward call to internal time dependent function
      void evaluate( const DomainType& x, RangeType& result ) const
      {
        function_.evaluate(x, time_, result );
      }

    protected:
      const FunctionType& function_;
      const double time_;
    };

    ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_HH
