#include "dMeshBlock.h"
#include "gpu_globals.h"

namespace TIOGA {

  void dMeshBlock::setData(TIOGA::MeshBlockInfo* m_info)
  {
    // Populate legacy data
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
    
    for (int i=0; i < TIOGA::MeshBlockInfo::max_vertex_types; ++i) {
      vconn_ptrs[i] = m_info->vertex_conn[i].dptr;
    }
    vconn=&vconn_ptrs[0];
  }

  void dMeshBlock::resetIblanks() {
    int n_blocks = nnodes/block_size + (nnodes%block_size == 0 ? 0:1);
    TIOGA_GPU_LAUNCH_FUNC(g_reset_iblanks, n_blocks, block_size, 0, 0, iblank, nnodes);
    //g_reset_iblanks<<<n_blocks,block_size>>>(iblank,nnodes);
  }
} //namespace TIOGA
