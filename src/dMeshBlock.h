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
  int myid;

  int ninterp,nweights;
  int *interpList_wcft{nullptr};
  int *interpList_inode{nullptr};
  double *interpList_weights{nullptr};
	
  dMeshBlock() {}

  ~dMeshBlock();

  void setMinfo( TIOGA::MeshBlockInfo* m_info_in, int myid_in);
  void resetIblanks();
  void search(ADT *adt,int *elementList, double *xsearch, int *donorId, int nsearch);

  void pushInterpListsToDevice(int ninterp_in, int nweights_in,
			       int *interpList_wcft_in,
			       int *interpList_inode_in,
			       double *interpList_weights_in);

  void getInterpolatedData(double *realData, int nvar);


};

    
} // namespace TIOGA

#endif 
