#ifndef DUNE_FEM_GENERIC_HH
#define DUNE_FEM_GENERIC_HH



// dune-fem
#include <dune/fem/space/generic/space.hh>

// dune-localfunctions
#include<dune/localfunctions/lagrange.hh>
#include<dune/localfunctions/utility/localfiniteelement.hh>

// TODO add more headers from localfunctions

namespace Dune
{
namespace Fem
{
	//-----------------------------
	//Dune-spaces
	//-----------------------------
	namespace Generic
	{
		template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
		class LagrangeDiscreteFunctionSpace
			: public Dune::Fem::GenericDiscreteFunctionSpace< GridPart, 
			                                                  LagrangeLocalFiniteElement< EquidistantPointSet,
																															                      FunctionSpace::dimDomain,
																											                            	typename FunctionSpace::DomainFieldType,
																											                            	typename FunctionSpace::RangeFieldType >
																																										
																												polOrder, 
																												Storage >
		{

		};
		
		
		template <class FunctionSpace, class GridPart, int polOrd, template <class> class Storage = CachingStorage >
    class LagrangeDiscontinuousGalerkinSpace
			: public Dune::Fem::GenericDiscreteFunctionSpace< GridPart, 
			                                                  DGLocalFiniteElement< LagrangeLocalFiniteElement< EquidistantPointSet,
																															                      FunctionSpace::dimDomain,
																											                            	typename FunctionSpace::DomainFieldType,
																											                            	typename FunctionSpace::RangeFieldType > >
																																										
																												polOrder, 
																												Storage >
		{

		};

		// TODO add more spaces
	
	}

}
}

#endif //end DUNE_FEM_GENERIC_HH
