#ifndef DMESHBLOCK_H
#define DMESHBLOCK_H
#include "TiogaMeshInfo.h"
#include "ADT.h"
namespace TIOGA {
class dMeshBlock
{
public:
  TIOGA::MeshBlockInfo* m_info_device{nullptr};
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
  int block_size{128};

  dMeshBlock() {}

  ~dMeshBlock() {};

  void setMinfo( TIOGA::MeshBlockInfo* m_info_in);
  void setData( TIOGA::MeshBlockInfo* m_info_device);
  void resetIblanks();
  void search(ADT *adt,int *elementList, double *xsearch, int *donorId, int nsearch);

};

} // namespace TIOGA

#endif 
