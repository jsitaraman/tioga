#ifndef DMESHBLOCK_H
#define DMESHBLOCK_H
#include "TiogaMeshInfo.h"
#include "ADT.h"
namespace TIOGA {
class dMeshBlock
{
public:
  TIOGA::MeshBlockInfo* m_info_device{nullptr};
  int myid{0};
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

  ~dMeshBlock() {};

  void setMinfo( TIOGA::MeshBlockInfo* m_info_in, int myid);
  void setData( TIOGA::MeshBlockInfo* m_info_device);
  void resetIblanks();
  void search(double *x,ADT *adt,int *elementList, double *xsearch, int *donorId, int nsearch, int nnodes, int ntypes, int *nc,int *nv, int **vconn);

};

} // namespace TIOGA

#endif 
