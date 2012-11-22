#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALRESTRICTPROLONG_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_LOCALRESTRICTPROLONG_HH

// C++ includes
#include <cassert>

// dune-fem includes
#include <dune/fem/space/common/localrestrictprolong.hh>

//local includes
#include "declaration.hh"


namespace Dune
{

  namespace Fem
  {

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    class DefaultLocalRestrictProlong< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
    {
      typedef DefaultLocalRestrictProlong< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > > ThisType;

    public:
      typedef RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > DiscreteFunctionSpaceType; 

      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType::RangeType RangeType;

      static const int dimension = GridPart::dimension;

      DefaultLocalRestrictProlong ( const DiscreteFunctionSpaceType &space )
      : space_( space )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight ) {}

      template< class FT, class ST, class LocalGeometry >
      void restrictLocal ( LocalFunction< FT > &lfFather, const LocalFunction< ST > &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const;

      template< class FT, class ST, class LocalGeometry >
      void prolongLocal ( const LocalFunction< FT > &lfFather, LocalFunction< ST > &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const;

      bool needCommunication () const { return false; }

    private:
      const DiscreteFunctionSpaceType &space_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALRESTRICTPROLONG_HH
