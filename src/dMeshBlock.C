#include "dMeshBlock.h"
#include "gpu_global_functions.h"
#include<chrono>
#include<thread>
\
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

  void dMeshBlock::resetIblanks() {
#ifdef TIOGA_HAS_GPU
    int n_blocks = nnodes/block_size + (nnodes%block_size == 0 ? 0:1);
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
  
} //namespace TIOGA
