#ifndef DUNE_FEM_SHAPEFUNCTIONSET_MAPPED_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_MAPPED_HH

//- C++ includes
#include <cstddef>
#include <iostream>

/**
  @file
  @brief Interface for shape function sets
*/
namespace Dune
{
  namespace Fem
  {
      
    template<class OriginalShapeFunctionSet, class Mapping>
    class MappedShapeFunctionSet
    {
      typedef MappedShapeFunctionSet<OriginalShapeFunctionSet,Mapping> ThisType;
   
    public:
    
      typedef OriginalShapeFunctionSet OriginalShapeFunctionSetType;
      typedef Mapping MappingType;

    protected:
      template<class Functor>
      struct MappedFunctor;
     
    public:
      typedef typename OriginalShapeFunctionSetType::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType; 
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;  
    public:
      explicit MappedShapeFunctionSet(const OriginalShapeFunctionSetType & originalShapeFunctionSet, 
                                      const MappingType& mapping = MappingType() )
        : originalShapeFunctionSet_(originalShapeFunctionSet),
          mapping_(mapping)
      { 
      }

      int order() const {return originalShapeFunctionSet_.order();}

      std::size_t size() const {return originalShapeFunctionSet_.size();}
     
      const MappingType& map() const { return mapping_; }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
       
      const OriginalShapeFunctionSetType& originalShapeFunctionSet() const { return originalShapeFunctionSet_; }

    protected:
      OriginalShapeFunctionSetType originalShapeFunctionSet_;
      MappingType mapping_;

    };
    
    
    
    //MappedShapeFunctionSet::MappedFunctor
    //------------------------------------
    template<class OriginalShapeFunctionSet, class Mapping>
    template<class Functor>
    struct MappedShapeFunctionSet<OriginalShapeFunctionSet, Mapping>::MappedFunctor
    {
      
      explicit MappedFunctor(Functor functor, const Mapping& map)
      :functor_(functor),
       map_(map)
      {}
      
      template< class Value >
      void operator()(const std::size_t i, const Value &v)
      {
        size_t reali=map_[i];
        functor_(reali,v);
      }  
    private:
      Functor functor_;
     const MappingType& map_;
    };


    //Implementation of MappedShapeFunctionSet
    //---------------------------------------
   
    template<class OriginalShapeFunctionSet, class Mapping>
    template< class Point, class Functor >
    inline void MappedShapeFunctionSet< OriginalShapeFunctionSet,Mapping>
      ::evaluateEach(const Point &x, Functor functor) const
      {
        originalShapeFunctionSet().evaluateEach(x,MappedFunctor<Functor>(functor,map()));
      }
    
    template<class OriginalShapeFunctionSet, class Mapping>
    template< class Point, class Functor >
    inline void MappedShapeFunctionSet< OriginalShapeFunctionSet,Mapping>
      ::jacobianEach(const Point &x, Functor functor) const
      {
        originalShapeFunctionSet().jacobianEach(x,MappedFunctor<Functor>(functor,map()));
      }

    template<class OriginalShapeFunctionSet, class Mapping>
    template< class Point, class Functor >
    inline void MappedShapeFunctionSet< OriginalShapeFunctionSet,Mapping>
      ::hessianEach(const Point &x, Functor functor) const
      {
        originalShapeFunctionSet().hessianEach(x,MappedFunctor<Functor>(functor,map()));
      }



  }// namespace FEM

}// namespace Dune



#endif //DUNE_FEM_SHAPEFUNCTIONSET_MAPPED.HH

