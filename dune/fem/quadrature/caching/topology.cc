#include "topology.hh"

namespace Dune
{
  namespace Fem
  {
    //- class FaceTopologyMapping
    template <>
    int FaceTopologyMapping<tetra>::
    twist(int index, int faceTwist) {
      return (faceTwist < 0) ?
        (7 - index + faceTwist)%3 : (faceTwist + index)%3 ;
    }

    template <>
    int FaceTopologyMapping<hexa>::
    twist(int index, int faceTwist) {
      return (faceTwist < 0) ?
        (9 - index + faceTwist)%4 : (faceTwist + index)%4 ;
    }

    template <>
    int FaceTopologyMapping<tetra>::
    invTwist(int index, int faceTwist)
    {
      return (faceTwist < 0) ?
        (7 - index + faceTwist)%3 : (3 + index - faceTwist)%3;
    }

    template <>
    int FaceTopologyMapping<hexa>::
    invTwist(int index, int faceTwist) {
      return (faceTwist < 0) ?
        (9 - index + faceTwist)%4 : (4 + index - faceTwist)%4;
    }

    // alu triangle face are oriented just the other way then dune faces
    // therefore vertex 1 and 2 are swapped because
    // ALUGrid tetra face are oriented just the other way compared to Dune
    // tetra faces, see also gitter_geo.cc of the ALUGrid code
    template <>
    const int FaceTopologyMapping<tetra>::
    dune2aluVertex_[EntityCount<tetra>::numVerticesPerFace] = {0, 2, 1};

    // mapping of twists from alu 2 dune
    template <>
    const int FaceTopologyMapping<tetra>::
    alu2duneTwist_[ 2 * EntityCount<tetra>::numVerticesPerFace] = { -2, -3, -1, 0, 2, 1 };

    // mapping of twists from alu 2 dune
    template <>
    const int FaceTopologyMapping<tetra>::
    aluTwistMap_[ 2 * EntityCount<tetra>::numVerticesPerFace] = { 1, 2, 0, -1, -2, -3 };

    // mapping of twists from alu 2 dune
    template <>
    const int FaceTopologyMapping<hexa>::
    aluTwistMap_[ 2 * EntityCount<hexa>::numVerticesPerFace] = { -2, -3, -4, -1, 0, 1, 2, 3 } ;

    // the mapping of vertices in the reference quad
    // this is used for hexa face during intersection iterator build
    // and to calculate the intersectionSelfLocal and
    // intersectionSelfNeighbor geometries.
    template <>
    const int FaceTopologyMapping<hexa>::
    dune2aluVertex_[EntityCount<hexa>::numVerticesPerFace] = {0, 3, 1, 2};

    // alu triangle face are oriented just the other way then dune faces
    // therefore vertex 1 and 2 are swaped
    template <>
    const int FaceTopologyMapping<tetra>::
    alu2duneVertex_[EntityCount<tetra>::numVerticesPerFace] = {0, 2, 1};

    template <>
    const int FaceTopologyMapping<hexa>::
    alu2duneVertex_[EntityCount<hexa>::numVerticesPerFace] = {0, 2, 3, 1};

    template <>
    const int FaceTopologyMapping<tetra>::
    dune2aluEdge_[EntityCount<tetra>::numEdgesPerFace] = {1, 2, 0};

    template <>
    const int FaceTopologyMapping<hexa>::
    dune2aluEdge_[EntityCount<hexa>::numEdgesPerFace] = {0, 2, 3, 1};

    template <>
    const int FaceTopologyMapping<tetra>::
    alu2duneEdge_[EntityCount<tetra>::numEdgesPerFace] = {2, 0, 1};

    template <>
    const int FaceTopologyMapping<hexa>::
    alu2duneEdge_[EntityCount<hexa>::numEdgesPerFace] = {0, 3, 1, 2};

  } // end namespace Fem
} // end namespace Dune
