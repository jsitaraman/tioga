#ifndef DMESHBLOCK_H
#define DMESHBLOCK_H
#include "TiogaMeshInfo.h"
namespace TIOGA {
class dMeshBlock
{
public:
  int meshtag{0};
  int nnodes{0};
  int nwbc{0};
  int nobc{0};
  int ncells{0};
  int ntypes{0};
  double *x{nullptr};
  int *iblank{nullptr};
  int *iblank_cell{nullptr};
  int *wbcnode{nullptr};
  int *obcnode{nullptr};
  int *nc{nullptr};
  int *nv{nullptr};
  int* vconn_ptrs[TIOGA::MeshBlockInfo::max_vertex_types];
  int **vconn;
  int block_size{512};

  dMeshBlock() {}

  ~dMeshBlock();

  void setData( TIOGA::MeshBlockInfo* m_info_device);
  void resetIblanks();

};

} // namespace TIOGA

#endif 
