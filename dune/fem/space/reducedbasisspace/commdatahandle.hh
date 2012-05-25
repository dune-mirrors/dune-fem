#ifndef DUNE_FEM_REDUCEDBASISSPACE_COMMHANDLE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_COMMHANDLE_HH

#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  template< class DiscreteFunction, class Operation >
  class ReducedBasisCommDataHandle
  : public CommDataHandleIF
    < ReducedBasisCommDataHandle< DiscreteFunction, Operation >,
      typename DiscreteFunction :: RangeFieldType >
  {
  public:
    typedef DiscreteFunction DiscreteFunctionType;

    typedef Operation OperationType;

    typedef typename DiscreteFunctionType :: RangeFieldType DataType;

  private:
    typedef ReducedBasisCommDataHandle< DiscreteFunctionType, OperationType >
      ThisType;
    typedef CommDataHandleIF< ThisType, DataType > BaseType;

  public:
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

  protected:
    DiscreteFunctionType *const function_;

  public:
    inline explicit
    ReducedBasisCommDataHandle ( DiscreteFunctionType &function )
    : function_( &function )
    {}

    inline
    ReducedBasisCommDataHandle ( const ThisType &other )
    : function_( other.function_ )
    {}

  private:
    ThisType &operator= ( const ThisType & );

  public:
    inline bool contains ( int dim, int codim ) const
    {
      return false;
    }

    inline bool fixedsize ( int dim, int codim ) const
    {
      return true;
    }

    template< class MessageBuffer, class Entity >
    void gather ( MessageBuffer &buffer, const Entity &entity ) const
    {
    }

    template< class MessageBuffer, class Entity >
    void scatter ( MessageBuffer &buffer, const Entity &entity, size_t n )
    {
    }
    
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      return 0;
    }
  };
  
}

#endif
