#ifndef DUNE_FEM_ALU3DGRIDTOPOLOGY_HH
#define DUNE_FEM_ALU3DGRIDTOPOLOGY_HH

//- system includes
#include <cassert>

namespace Dune
{

  namespace Fem
  {

    // types of the elementes,
    // i.e . tetra or hexa, mixed is not implemeneted
    enum ALU3dGridElementType { tetra = 4, hexa = 7, mixed, error };

    template <ALU3dGridElementType type>
    struct EntityCount {};

    template <>
    struct EntityCount<tetra> {
      enum {numFaces = 4};
      enum {numVertices = 4};
      enum {numEdges = 6};
      enum {numVerticesPerFace = 3};
      enum {numEdgesPerFace = 3};
    };

    template <>
    struct EntityCount<hexa> {
      enum {numFaces = 6};
      enum {numVertices = 8};
      enum {numEdges = 12};
      enum {numVerticesPerFace = 4};
      enum {numEdgesPerFace = 4};
    };


    //! Maps indices of the Dune reference face onto the indices of the
    //! ALU3dGrid reference face and vice-versa.
    template <ALU3dGridElementType type>
    class FaceTopologyMapping {
    public:
      //! Maps vertex index from Dune onto ALU3dGrid reference face
      static int dune2aluVertex(int index);
      //! Maps vertex index from Dune onto ALU3dGrid reference face, where the
      //! face in the ALU3dGrid has the twist <i>twist</i> compared to the orientation
      //! of the respective face in the reference element
      //! \param index local Dune vertex index on the particular face (i.e. the
      //! face which has a twist <i>twist</i> compared to the reference element's face
      //! \param twist twist of the face in consideration
      //! \return local ALU3dGrid vertex index on reference element face
      static int dune2aluVertex(int index, int twist);
      //! Maps vertex index from ALU3dGrid onto Dune reference face
      static int alu2duneVertex(int index);
      //! Maps vertex index from ALU3dGrid onto Dune reference face, where the
      //! face in the ALU3dGrid has the twist <i>twist</i> compared to the orientation
      //! of the respective face in the reference element
      //! \param index local ALU3dGrid vertex index on the particular face (i.e.
      //! the face which has a twist <i>twist</i> compared to the reference element's
      //! face
      //! \param twist twist of the face in consideration
      //! \return local Dune vertex index on reference element face
      static int alu2duneVertex(int index, int twist);
      //! Maps edge index from Dune onto ALU3dGrid reference face
      static int dune2aluEdge(int index);
      //! Maps edge index from ALU3dGrid onto Dune reference face
      static int alu2duneEdge(int index);
      //  private:
      static int twist(int index, int faceTwist);
      static int invTwist(int index, int faceTwist);

      static int twistedDuneIndex( const int idx, const int twist );

      // for each aluTwist apply additional mapping
      static int aluTwistMap(const int aluTwist);
    private:
      const static int dune2aluVertex_[EntityCount<type>::numVerticesPerFace];
      const static int alu2duneVertex_[EntityCount<type>::numVerticesPerFace];

      const static int dune2aluEdge_[EntityCount<type>::numEdgesPerFace];
      const static int alu2duneEdge_[EntityCount<type>::numEdgesPerFace];

      const static int alu2duneTwist_[ 2 * EntityCount<type>::numVerticesPerFace ];
      const static int aluTwistMap_[ 2 * EntityCount<type>::numVerticesPerFace ];
    };

    //- class FaceTopologyMapping
    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::dune2aluVertex(int index) {
      assert(index >= 0 && index < EntityCount<type>::numVerticesPerFace);
      return dune2aluVertex_[index];
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::dune2aluVertex(int index, int twist) {
      assert(index >= 0 && index < EntityCount<type>::numVerticesPerFace);
      return invTwist(dune2aluVertex_[index], twist);
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::alu2duneVertex(int index) {
      assert(index >= 0 && index < EntityCount<type>::numVerticesPerFace);
      return alu2duneVertex_[index];
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::alu2duneVertex(int index, int twist)
    {
      assert(index >= 0 && index < EntityCount<type>::numVerticesPerFace);
      return alu2duneVertex_[invTwist(index, twist)];
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::alu2duneEdge(int index) {
      assert(index >= 0 && index < EntityCount<type>::numEdgesPerFace);
      return alu2duneEdge_[index];
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::
    aluTwistMap(const int aluTwist)
    {
      // this map has been calculated by grid/test/checktwists.cc
      // and the dune-fem twist calculator
      // this should be revised after the release 2.1
      return aluTwistMap_[ aluTwist + ((type == tetra) ? 3 : 4) ];
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::
    twistedDuneIndex(const int duneIdx, const int aluTwist)
    {
      if( type == tetra )
      {
        // apply alu2dune twist mapping (only for tetra)
        const int twist = alu2duneTwist_[ aluTwist + 3 ];
        return alu2duneVertex( dune2aluVertex(duneIdx) , twist );
      }
      else
       return alu2duneVertex( dune2aluVertex(duneIdx) , aluTwist );
    }

    template <ALU3dGridElementType type>
    inline int FaceTopologyMapping<type>::dune2aluEdge(int index) {
      assert(index >= 0 && index < EntityCount<type>::numEdgesPerFace);
      return dune2aluEdge_[index];
    }


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALU3DGRIDTOPOLOGY_HH
