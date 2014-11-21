#ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
#define DUNE_FEM_COMBINEDADAPTMANAGER_HH

#include <dune/common/exceptions.hh>
//- local includes
#include <dune/fem/common/subvector.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>



namespace Dune
{

  namespace Fem
  {

    template< class , class >
    class CombinedDiscreteFunctionSpace;

    template< class SP1, class SP2 >
    class CombinedLocalRestrictProlong
    {
      typedef CombinedLocalRestrictProlong< SP1, SP2> ThisType;

      typedef DefaultLocalRestrictProlong< SP1 > LocalRestrictProlong1Type;
      typedef DefaultLocalRestrictProlong< SP2 > LocalRestrictProlong2Type;

    public:
      typedef SP1 DiscreteFunctionSpaceType1;
      typedef SP2 DiscreteFunctionSpaceType2;

      typedef typename DiscreteFunctionSpaceType1::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType1::RangeFieldType RangeFieldType;
      typedef typename DiscreteFunctionSpaceType1::RangeType RangeType1;
      typedef typename DiscreteFunctionSpaceType2::RangeType RangeType2;

      typedef typename DiscreteFunctionSpaceType1::GridPartType GridPartType;

      CombinedLocalRestrictProlong ( const CombinedDiscreteFunctionSpace< SP1, SP2> &space )
      : rp1_( space.space1() ),
        rp2_( space.space2() )
      {}

      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        rp1_.setFatherChildWeight( weight );
        rp2_.setFatherChildWeight( weight );
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        typedef DenseSubVector< typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;
        typedef DenseSubVector< const typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;

        typedef typename DiscreteFunctionSpaceType1 :: BasisFunctionSetType BasisFunctionSetType1;
        typedef typename DiscreteFunctionSpaceType2 :: BasisFunctionSetType BasisFunctionSetType2;

        const BasisFunctionSetType1 &fatherBasisFunctionSet1 = lfFather.basisFunctionSet().template subBasisFunctionSet< 0 >();
        const BasisFunctionSetType2 &fatherBasisFunctionSet2 = lfFather.basisFunctionSet().template subBasisFunctionSet< 1 >();

        const BasisFunctionSetType1 &sonBasisFunctionSet1 = lfSon.basisFunctionSet().template subBasisFunctionSet< 0 >();
        const BasisFunctionSetType2 &sonBasisFunctionSet2 = lfSon.basisFunctionSet().template subBasisFunctionSet< 1 >();

        const std::size_t fatherSize1 = fatherBasisFunctionSet1.size();
        const std::size_t fatherSize2 = fatherBasisFunctionSet2.size();
        const std::size_t fatherOffset1 = lfFather.basisFunctionSet().template offset<0>();
        const std::size_t fatherOffset2 = lfFather.basisFunctionSet().template offset<1>();

        const std::size_t sonSize1 = sonBasisFunctionSet1.size();
        const std::size_t sonSize2 = sonBasisFunctionSet2.size();
        const std::size_t sonOffset1 = lfSon.basisFunctionSet().template offset<0>();
        const std::size_t sonOffset2 = lfSon.basisFunctionSet().template offset<1>();

        LocalFunction< BasisFunctionSetType1, SubDofVectorTypeFather > lfFather1 (
            fatherBasisFunctionSet1, SubDofVectorTypeFather( lfFather.localDofVector(), fatherSize1, fatherOffset1 ) );
        LocalFunction< BasisFunctionSetType2, SubDofVectorTypeFather > lfFather2 (
            fatherBasisFunctionSet2, SubDofVectorTypeFather( lfFather.localDofVector(), fatherSize2, fatherOffset2 ) );

        BasicConstLocalFunction< BasisFunctionSetType1, SubDofVectorTypeSon > lfSon1 (
            sonBasisFunctionSet1, SubDofVectorTypeSon( lfSon.localDofVector(), sonSize1, sonOffset1 ) );
        BasicConstLocalFunction< BasisFunctionSetType2, SubDofVectorTypeSon > lfSon2 (
            sonBasisFunctionSet2, SubDofVectorTypeSon( lfSon.localDofVector(), sonSize2, sonOffset2 ) );

        rp1_.restrictLocal( lfFather1, lfSon1, geometryInFather, initialize );
        rp2_.restrictLocal( lfFather2, lfSon2, geometryInFather, initialize );
      }


      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        typedef DenseSubVector< const typename LFFather::LocalDofVectorType > SubDofVectorTypeFather;
        typedef DenseSubVector< typename LFSon::LocalDofVectorType > SubDofVectorTypeSon;

        typedef typename DiscreteFunctionSpaceType1 :: BasisFunctionSetType BasisFunctionSetType1;
        typedef typename DiscreteFunctionSpaceType2 :: BasisFunctionSetType BasisFunctionSetType2;

        const BasisFunctionSetType1 &fatherBasisFunctionSet1 = lfFather.basisFunctionSet().template subBasisFunctionSet< 0 >();
        const BasisFunctionSetType2 &fatherBasisFunctionSet2 = lfFather.basisFunctionSet().template subBasisFunctionSet< 1 >();

        const BasisFunctionSetType1 &sonBasisFunctionSet1 = lfSon.basisFunctionSet().template subBasisFunctionSet< 0 >();
        const BasisFunctionSetType2 &sonBasisFunctionSet2 = lfSon.basisFunctionSet().template subBasisFunctionSet< 1 >();

        const std::size_t fatherSize1 = fatherBasisFunctionSet1.size();
        const std::size_t fatherSize2 = fatherBasisFunctionSet2.size();
        const std::size_t fatherOffset1 = lfFather.basisFunctionSet().template offset<0>();
        const std::size_t fatherOffset2 = lfFather.basisFunctionSet().template offset<1>();

        const std::size_t sonSize1 = sonBasisFunctionSet1.size();
        const std::size_t sonSize2 = sonBasisFunctionSet2.size();
        const std::size_t sonOffset1 = lfSon.basisFunctionSet().template offset<0>();
        const std::size_t sonOffset2 = lfSon.basisFunctionSet().template offset<1>();

        LocalFunction< BasisFunctionSetType1, SubDofVectorTypeSon > lfSon1 (
            sonBasisFunctionSet1, SubDofVectorTypeSon( lfSon.localDofVector(), sonSize1, sonOffset1 ) );
        LocalFunction< BasisFunctionSetType2, SubDofVectorTypeSon > lfSon2 (
            sonBasisFunctionSet2, SubDofVectorTypeSon( lfSon.localDofVector(), sonSize2, sonOffset2 ) );

        BasicConstLocalFunction< BasisFunctionSetType1, SubDofVectorTypeFather > lfFather1 (
            fatherBasisFunctionSet1, SubDofVectorTypeFather( lfFather.localDofVector(), fatherSize1, fatherOffset1 ) );
        BasicConstLocalFunction< BasisFunctionSetType2, SubDofVectorTypeFather > lfFather2 (
            fatherBasisFunctionSet2, SubDofVectorTypeFather( lfFather.localDofVector(), fatherSize2, fatherOffset2 ) );

        rp1_.prolongLocal( lfFather1, lfSon1, geometryInFather, initialize );
        rp2_.prolongLocal( lfFather2, lfSon2, geometryInFather, initialize );
      }

      bool needCommunication () const { return  rp1_.needCommunication() || rp2_.needCommunication(); }

    private:
      LocalRestrictProlong1Type rp1_;
      LocalRestrictProlong2Type rp2_;
    };


    template< class SP1, class SP2 >
    class DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > >
    : public CombinedLocalRestrictProlong< SP1, SP2 >
    {
      typedef DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > > ThisType;
      typedef CombinedLocalRestrictProlong< SP1, SP2 > BaseType;
    public:

      DefaultLocalRestrictProlong ( const CombinedDiscreteFunctionSpace< SP1, SP2 > & space )
      : BaseType( space )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif //  #ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
