#include "dMeshBlock.h"
#include "gpu_global_functions.h"
#include<chrono>
#include<thread>
\
namespace TIOGA {

  void dMeshBlock::setMinfo(TIOGA::MeshBlockInfo *m_info_in, int myid_in)
  {
    myid=myid_in;
    m_info_device=m_info_in;
  }

  void dMeshBlock::setData(TIOGA::MeshBlockInfo* m_info)
  {
    // Populate legacy data
#if defined(TIOGA_HAS_GPU) 
 #if defined(TIOGA_FAKE_GPU)
    meshtag = m_info->meshtag;
    nnodes = m_info->num_nodes;
    x = m_info->xyz.dptr;
    iblank = m_info->iblank_node.dptr;
    iblank_cell = m_info->iblank_cell.dptr;
    nwbc = m_info->wall_ids.sz;
    nobc = m_info->overset_ids.sz;
    wbcnode = m_info->wall_ids.dptr;
    obcnode = m_info->overset_ids.dptr;
    
    ntypes = m_info->num_vert_per_elem.sz;
    nv = m_info->num_vert_per_elem.dptr;
    nc = m_info->num_cells_per_elem.dptr;
    
    ncells=0;
    for(int i=0;i<ntypes;i++) ncells+=nc[i];
    
    for (int i=0; i < TIOGA::MeshBlockInfo::max_vertex_types; i++) {
      vconn_ptrs[i] = m_info->vertex_conn[i].dptr;
    }
    vconn=&vconn_ptrs[0];
 #else
    meshtag = m_info->meshtag;
    nnodes = m_info->num_nodes;
    x = m_info->xyz.hptr;
    iblank = m_info->iblank_node.dptr;
    iblank_cell = m_info->iblank_cell.dptr;
    nwbc = m_info->wall_ids.sz;
    nobc = m_info->overset_ids.sz;
    wbcnode = m_info->wall_ids.dptr;
    obcnode = m_info->overset_ids.dptr;
    
    ntypes = m_info->num_vert_per_elem.sz;
    nv = m_info->num_vert_per_elem.dptr;
    nc = m_info->num_cells_per_elem.dptr;
//    printf("ntypes=%d\n",ntypes);
//    for(int i=0;i<ntypes;i++) printf("nv[i],nc[i]=%d %d\n",nv[i],nc[i]);  
    ncells=0;
    for(int i=0;i<ntypes;i++) ncells+=m_info->num_cells_per_elem.hptr[i];
    
    for (int i=0; i < TIOGA::MeshBlockInfo::max_vertex_types; i++) {
      vconn_ptrs[i] = m_info->vertex_conn[i].dptr;
      //if (i==0) ndc4= m_info->vertex_conn[i].hptr;
      //if (i==1) ndc5= m_info->vertex_conn[i].hptr;
      //if (i==2) ndc6= m_info->vertex_conn[i].hptr;
      //if (i==0) ndc8= m_info->vertex_conn[i].hptr;
    }
    vconn=&vconn_ptrs[0];
  #endif
#endif
  }

  void dMeshBlock::resetIblanks() {
#ifdef TIOGA_HAS_GPU
    int n_blocks = nnodes/block_size + (nnodes%block_size == 0 ? 0:1);
//    TIOGA_GPU_LAUNCH_FUNC(g_reset_iblanks, n_blocks, block_size, 0, 0, iblank, nnodes);
    TIOGA_GPU_LAUNCH_FUNC(g_reset_iblanks, n_blocks, block_size, 0, 0, m_info_device);
    TIOGA::gpu::synchronize();
    //g_reset_iblanks<<<n_blocks,block_size>>>(iblank,nnodes);
#endif
  }

  void dMeshBlock::search(double *xyz, ADT *adt,int *elementList_host, double *xsearch_host, int *donorId_host, 
                          int nsearch, int nnodes, int ntypes,int *nc_cpu, 
                          int *nv_cpu, int **vconn_cpu)
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
   //int    *donorId=TIOGA::gpu::push_to_device<int>(donorId_host,sizeof(int)*nsearch);
   int    *donorId=TIOGA::gpu::allocate_on_device<int>(sizeof(int)*nsearch);
   double *xsearch=TIOGA::gpu::push_to_device<double>(xsearch_host,sizeof(double)*nsearch*3);
   double *x=TIOGA::gpu::push_to_device<double>(xyz,sizeof(double)*nnodes*3);
   //int *vconn[4];
   //for(int i=0;i<ntypes;i++)
   //  vconn[i]=TIOGA::gpu::push_to_device<int>(vconn_cpu[i],sizeof(int)*nc_cpu[i]*nv_cpu[i]);
   int *vconn=TIOGA::gpu::push_to_device<int>(vconn_cpu[0],sizeof(int)*nc_cpu[0]*nv_cpu[0]);
   // 
   int *nc=TIOGA::gpu::push_to_device<int>(nc_cpu,sizeof(int)*ntypes);
   int *nv=TIOGA::gpu::push_to_device<int>(nv_cpu,sizeof(int)*ntypes);
   //
   // perform the gpu based adt search now
   //  
   int n_blocks = nsearch/block_size + (nsearch%block_size == 0 ? 0:1);
   TIOGA_GPU_LAUNCH_FUNC(g_adt_search,n_blocks,block_size,0,0,x,m_info_device,ntypes,nc,nv,vconn,
                         coord,adtExtents,adtIntegers,adtReals,
                         elementList,donorId,xsearch,ndim,nelem,nsearch);

   //g_adt_search(double *x, int **vconn,int *nc, int *nv, int ntypes, 
   //             double *coord, double *adtExtents, double *adtIntegers, double *adtReals,
   //             int *elementList, double *donorId, double *xsearch, int ndim, int nelem, int nsearch);
   TIOGA::gpu::synchronize();
   TIOGA::gpu::pull_from_device<int>(donorId_host,donorId,sizeof(int)*nsearch);    
   printf("myid/donorId_host[0]=%d %d\n",myid,donorId_host[0]);
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
   TIOGA_FREE_DEVICE(x);
   TIOGA_FREE_DEVICE(vconn);
   TIOGA_FREE_DEVICE(nc);
   TIOGA_FREE_DEVICE(nv);
#endif
  }
  
} //namespace TIOGA
