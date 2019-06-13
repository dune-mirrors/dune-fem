#ifndef DUNE_FEMPY_PY_SPACEADAPT_HH
#define DUNE_FEMPY_PY_SPACEADAPT_HH

#include <dune/fempy/pybind11/pybind11.hh>
#include <dune/fempy/space/adaptation.hh>

namespace Dune
{

  namespace FemPy
  {

    template <class DF, class... options >
    void registerSpaceAdaptation( pybind11::module module, pybind11::class_<SpaceAdaptation<DF>, options... > cls )
    {
      typedef SpaceAdaptation<DF> SpaceAdapt;
      typedef typename SpaceAdapt::DiscreteFunctionSpaceType Space;
      typedef typename Space::EntityType Element;
      using pybind11::operator""_a;
      cls.def(pybind11::init([] (Space &space)
            {return new SpaceAdapt(space);}
            ),pybind11::keep_alive<1,2>(), "space"_a );
      cls.def( "adapt", [] ( SpaceAdapt &self, std::function<int(Element)> marker, const std::list< std::reference_wrapper<DF> > &dfList ) {
          self.adapt( marker, dfList.begin(), dfList.end() );
        } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACEADAPT_HH
