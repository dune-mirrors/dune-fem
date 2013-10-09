#ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
#define DUNE_FEM_COMBINEDADAPTMANAGER_HH

#include <dune/common/exceptions.hh>
//- local includes  
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/combinedspace/localsubfunction.hh>


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
      template< class FT, class ST, class LocalGeometry >
      void restrictLocal ( LocalFunction< FT > &lfFather, const LocalFunction< ST > &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        LocalSubFunction< FT, 0 >  lfFather1( lfFather );
        ConstLocalSubFunction< ST, 0 >  lfSon1( lfSon );
        LocalSubFunction< FT, 1 >  lfFather2( lfFather );
        ConstLocalSubFunction< ST, 1 >  lfSon2( lfSon );

        rp1_.restrictLocal( lfFather1, lfSon1, geometryInFather, initialize );
        rp2_.restrictLocal( lfFather2, lfSon2, geometryInFather, initialize );
      }


      template< class FT, class ST, class LocalGeometry >
      void prolongLocal ( const LocalFunction< FT > &lfFather, LocalFunction< ST > &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        ConstLocalSubFunction< FT, 0 >  lfFather1( lfFather );
        LocalSubFunction< ST, 0 >  lfSon1( lfSon );
        ConstLocalSubFunction< FT, 1 >  lfFather2( lfFather );
        LocalSubFunction< ST, 1 >  lfSon2( lfSon );

        rp1_.prolongLocal( lfFather1, lfSon1, geometryInFather, initialize );
        rp2_.prolongLocal( lfFather2, lfSon2, geometryInFather, initialize );
      }

      bool needCommunication () const { return  rp1_.needCommunication() || rp2_.needCommunication(); }

    private:
      LocalRestrictProlong1Type rp1_;
      LocalRestrictProlong2Type rp2_;
    };


    template< class SP1, class SP2 >
    struct DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > >
    : public CombinedLocalRestrictProlong< SP1, SP2 >
    {
      typedef DefaultLocalRestrictProlong< CombinedDiscreteFunctionSpace< SP1, SP2 > > ThisType;
      typedef CombinedLocalRestrictProlong< SP1, SP2 > BaseType;

      DefaultLocalRestrictProlong ( const CombinedDiscreteFunctionSpace< SP1, SP2 > & space )
      : BaseType( space )
      {}
    };

  } // namespace Fem

} // namespace Dune 

#endif //  #ifndef DUNE_FEM_COMBINEDADAPTMANAGER_HH
