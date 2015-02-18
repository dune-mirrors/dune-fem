#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_POWERLOCALRESTRICPROLONG_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_POWERLOCALRESTRICPROLONG_HH

#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dune/fem/common/subvector.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

#include <dune/fem/space/common/localrestrictprolong.hh>


namespace Dune
{

  namespace Fem
  {

    // PowerLocalRestricProlong
    // ------------------------

    template< class DiscreteFunctionSpace, int N >
    class PowerLocalRestrictProlong
    {
      typedef PowerLocalRestrictProlong< DiscreteFunctionSpace, N > ThisType;

      // type of contained DefaultLocalRestrictProlong
      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpace > LocalRestrictProlongType;

    public:
      // type of DomainField
      typedef typename LocalRestrictProlongType::DomainFieldType DomainFieldType;

      PowerLocalRestrictProlong ( const DiscreteFunctionSpace &space )
        : localRestrictProlong_( space )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        localRestrictProlong_.setFatherChildWeight( weight );
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        typedef DenseSubVector< const typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;
        typedef DenseSubVector< typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;

        typedef typename LFSon::BasisFunctionSetType::ScalarBasisFunctionSetType SubSonBasisFunctionSetType;
        typedef typename LFFather::BasisFunctionSetType::ScalarBasisFunctionSetType SubFatherBasisFunctionSetType;

        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().scalarBasisFunctionSet();
        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().scalarBasisFunctionSet();

        for( std::size_t i = 0; i < N; ++i )
        {
          std::size_t sonBasisSetSize = subSonBasisFunctionSet.size();
          std::size_t fatherBasisSetsize = subFatherBasisFunctionSet.size();

          SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), sonBasisSetSize, sonBasisSetSize * i );
          SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), fatherBasisSetsize, fatherBasisSetsize * i );

          BasicConstLocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon >
          subLFSon( subSonBasisFunctionSet, sonSubDofVector );
          LocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather >
          subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );

          localRestrictProlong_.restrictLocal( subLFFather, subLFSon, geometryInFather, initialize );
        }
      }


      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        typedef DenseSubVector< typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;
        typedef DenseSubVector< const typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;

        typedef typename LFSon::BasisFunctionSetType::ScalarBasisFunctionSetType SubSonBasisFunctionSetType;
        typedef typename LFFather::BasisFunctionSetType::ScalarBasisFunctionSetType SubFatherBasisFunctionSetType;

        SubSonBasisFunctionSetType subSonBasisFunctionSet = lfSon.basisFunctionSet().scalarBasisFunctionSet();
        SubFatherBasisFunctionSetType subFatherBasisFunctionSet = lfFather.basisFunctionSet().scalarBasisFunctionSet();

        for( std::size_t i = 0; i < N; ++i )
        {
          std::size_t sonBasisSetSize = subSonBasisFunctionSet.size();
          std::size_t fatherBasisSetsize = subFatherBasisFunctionSet.size();

          SubDofVectorTypeSon sonSubDofVector( lfSon.localDofVector(), sonBasisSetSize, sonBasisSetSize * i );
          SubDofVectorTypeFather fatherSubDofVector( lfFather.localDofVector(), fatherBasisSetsize, fatherBasisSetsize * i );

          LocalFunction< SubSonBasisFunctionSetType, SubDofVectorTypeSon >
          subLFSon( subSonBasisFunctionSet, sonSubDofVector );
          BasicConstLocalFunction< SubFatherBasisFunctionSetType, SubDofVectorTypeFather >
          subLFFather( subFatherBasisFunctionSet, fatherSubDofVector );

          localRestrictProlong_.prolongLocal( subLFFather, subLFSon, geometryInFather, initialize );
        }
      }

      bool needCommunication () const
      {
        return localRestrictProlong_.needCommunication();
      }

    private:
      LocalRestrictProlongType localRestrictProlong_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMBINEDSPACE_POWERLOCALRESTRICPROLONG_HH
