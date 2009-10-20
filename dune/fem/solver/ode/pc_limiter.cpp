#include "limiter.hpp"
#include "mesh_base.hpp"
#include "dg.hpp"
#include "blas.hpp"


using namespace pardg;



PCLimiter::PCLimiter(int dim_system, int num_base_polys, 
		     MeshBase &mesh, DG &dg) :
  dim_system(dim_system), num_base_polys(num_base_polys), mesh(mesh), dg(dg)
{
  dg.provide_smoothness_indicator(true);
}



void PCLimiter::operator()(double *U)
{
  const double *smoothness = dg.smoothness();

  for(int i=0; i<mesh.max_number_of_cells(); i++){

    if (smoothness[i] > 1.0){ // not smooth
    //if (true){
      //std::cout << i << " " << smoothness[i] << std::endl;

      double *_Uhi = U + i*dim_system*num_base_polys + dim_system;
      dset( (num_base_polys-1)*dim_system, 0.0, _Uhi, 1 );
    }
  }
}
