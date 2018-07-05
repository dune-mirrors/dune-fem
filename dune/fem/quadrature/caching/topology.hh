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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALU3DGRIDTOPOLOGY_HH
