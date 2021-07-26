#include "dMeshBlock.h"
#include "gpu_global_functions.h"
#include<chrono>
#include<thread>
#include <vector>

namespace TIOGA {

  void dMeshBlock::setMinfo(TIOGA::MeshBlockInfo *m_info_in, int myid_in)
  {
    if (m_info_device == nullptr) {
      m_info_device = TIOGA::gpu::allocate_on_device<TIOGA::MeshBlockInfo>(
        sizeof(TIOGA::MeshBlockInfo));
    }
    TIOGA::gpu::copy_to_device(m_info_device, m_info_in, sizeof(TIOGA::MeshBlockInfo));

    myid=myid_in;
  }

  // destructor that deallocates all the
  // the dynamic objects inside
  //
  dMeshBlock::~dMeshBlock()
  {
    if (m_info_device) TIOGA_FREE_DEVICE(m_info_device);
  };

  void dMeshBlock::update_minfo_device(TIOGA::MeshBlockInfo *m_info_in) {
#ifdef TIOGA_HAS_GPU
    TIOGA::gpu::copy_to_device(m_info_device, m_info_in, sizeof(TIOGA::MeshBlockInfo));
#endif
  }

  void dMeshBlock::resetIblanks(const int num_nodes) {
#ifdef TIOGA_HAS_GPU
    int n_blocks = num_nodes/block_size + (num_nodes%block_size == 0 ? 0:1);
    TIOGA_GPU_LAUNCH_FUNC(g_reset_iblanks, n_blocks, block_size, 0, 0, m_info_device);
    TIOGA::gpu::synchronize();
#endif
  }

  void dMeshBlock::search(ADT *adt,int *elementList_host, double *xsearch_host, int *donorId_host, 
                          int nsearch)
  {
#ifdef TIOGA_HAS_GPU
   int nelem=adt->get_nelem();
   int ndim=adt->get_ndim();
   // 
   // create gpu memory
   //
   double *adtExtents =TIOGA::gpu::push_to_device<double>(adt->get_Extents(),sizeof(double)*ndim);
   int    *adtIntegers=TIOGA::gpu::push_to_device<int>(adt->get_Integers(),sizeof(int)*nelem*4);
   double *adtReals=TIOGA::gpu::push_to_device<double>(adt->get_Reals(),sizeof(double)*nelem*ndim);
   double *coord=TIOGA::gpu::push_to_device<double>(adt->get_coord(),sizeof(double)*nelem*ndim);
   int    *elementList=TIOGA::gpu::push_to_device<int>(elementList_host,sizeof(int)*nelem);
   int    *donorId=TIOGA::gpu::push_to_device<int>(donorId_host,sizeof(int)*nsearch);
   double *xsearch=TIOGA::gpu::push_to_device<double>(xsearch_host,sizeof(double)*nsearch*3);
   //
   // perform the gpu based adt search now
   //  
   int n_blocks = nsearch/block_size + (nsearch%block_size == 0 ? 0:1);
   TIOGA_GPU_LAUNCH_FUNC(g_adt_search,n_blocks,block_size,0,0,m_info_device,
                         coord,adtExtents,adtIntegers,adtReals,
                         elementList,donorId,xsearch,ndim,nelem,nsearch,myid);

   //g_adt_search(double *x, int **vconn,int *nc, int *nv, int ntypes, 
   //             double *coord, double *adtExtents, double *adtIntegers, double *adtReals,
   //             int *elementList, double *donorId, double *xsearch, int ndim, int nelem, int nsearch);
   TIOGA::gpu::synchronize();
   TIOGA::gpu::pull_from_device<int>(donorId_host,donorId,sizeof(int)*nsearch);    
   //
   //release all gpu memory
   // 
   TIOGA_FREE_DEVICE(adtExtents);
   TIOGA_FREE_DEVICE(adtIntegers);
   TIOGA_FREE_DEVICE(adtReals);
   TIOGA_FREE_DEVICE(coord);
   TIOGA_FREE_DEVICE(elementList); 
   TIOGA_FREE_DEVICE(donorId);
   TIOGA_FREE_DEVICE(xsearch);
#endif
  }

  void dMeshBlock::pushInterpListsToDevice(int ninterp_in, int nweights_in,
					   int *interpList_wcft_in,
					   int *interpList_inode_in,
					   double *interpList_weights_in)
  {
    
    if (interpList_wcft) TIOGA_FREE_DEVICE(interpList_wcft);
    if (interpList_inode) TIOGA_FREE_DEVICE(interpList_inode);
    if (interpList_weights) TIOGA_FREE_DEVICE(interpList_weights);
    
    ninterp=ninterp_in;
    nweights=nweights_in;

    interpList_wcft=TIOGA::gpu::push_to_device<int>(interpList_wcft_in,
						    sizeof(int)*(ninterp+1));
    interpList_weights=TIOGA::gpu::push_to_device<double>(interpList_weights_in,
							  sizeof(double)*nweights);
    interpList_inode=TIOGA::gpu::push_to_device<int>(interpList_inode_in,
						     sizeof(int)*nweights);
  }

  void dMeshBlock::getInterpolatedData(double *realData,  int nvar, TIOGA::MeshBlockInfo *m_info_in)
  {
#ifdef TIOGA_HAS_GPU
    int n_blocks = nvar*ninterp/block_size 
      + ((nvar*ninterp)%block_size == 0 ? 0:1);
    double *realData_d=TIOGA::gpu::allocate_on_device<double>(sizeof(double)*ninterp*nvar);
    TIOGA_GPU_LAUNCH_FUNC(g_interp_data, n_blocks, block_size, 0, 0, 
                          interpList_wcft,
                          interpList_weights,
                          interpList_inode,
                          ninterp,
                          nvar,
                          realData_d,
                          m_info_device);
    TIOGA::gpu::synchronize();
    TIOGA::gpu::pull_from_device<double>(realData,realData_d,sizeof(double)*ninterp*nvar);    
    TIOGA_FREE_DEVICE(realData_d);
#endif    
  }

  void dMeshBlock::updateSolution(
    std::vector<int>& q_fringe_ind,
    std::vector<double>& q_fringe,
    TIOGA::MeshBlockInfo* m_info_in)
  {
#ifdef TIOGA_HAS_GPU
    int num_updates = q_fringe_ind.size();
    int n_blocks = num_updates/block_size + (num_updates%block_size == 0 ? 0:1);

    // create gpu memory
    int* fringe_ind_d =
        TIOGA::gpu::push_to_device<int>(q_fringe_ind.data(), sizeof(int)*num_updates);
    double* fringe_val_d =
        TIOGA::gpu::push_to_device<double>(q_fringe.data(), sizeof(double)*num_updates);

    TIOGA_GPU_LAUNCH_FUNC(g_update_sol, n_blocks, block_size, 0, 0,
                          num_updates,
                          fringe_ind_d,
                          fringe_val_d,
                          m_info_device);

    TIOGA::gpu::synchronize();

    TIOGA_FREE_DEVICE(fringe_ind_d);
    TIOGA_FREE_DEVICE(fringe_val_d);
#endif
  }
} //namespace TIOGA
