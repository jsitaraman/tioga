#ifndef DMESHBLOCK_H
#define DMESHBLOCK_H
#include "TiogaMeshInfo.h"
#include "ADT.h"
namespace TIOGA {
class dMeshBlock
{
public:
  TIOGA::MeshBlockInfo* m_info_device{nullptr};

  int block_size{128};
  int myid;
  int ninterp,nweights;
  int *interpList_wcft{nullptr};
  int *interpList_inode{nullptr};
  double *interpList_weights{nullptr};
	
  dMeshBlock() {}

  ~dMeshBlock();

  void setMinfo( TIOGA::MeshBlockInfo* m_info_in, int myid_in);
  void resetIblanks(const int num_nodes);
  void search(ADT *adt,int *elementList, double *xsearch, int *donorId, int nsearch);

  void pushInterpListsToDevice(int ninterp_in, int nweights_in,
			       int *interpList_wcft_in,
			       int *interpList_inode_in,
			       double *interpList_weights_in);

  void getInterpolatedData(double *realData, int nvar, TIOGA::MeshBlockInfo* m_info_in);

  void updateSolution(std::vector<int>& q_fringe_ind, std::vector<double>& q_fringe);
};

    
} // namespace TIOGA

#endif 
