#ifndef TIOGAMESHINFO_H
#define TIOGAMESHINFO_H

#include <cstdlib>
#include <cstdint>

namespace TIOGA {

template <typename T>
struct TiogaView
{
    T* hptr{nullptr};
    T* dptr{nullptr};

    size_t sz{0};
    bool host_modified{false};
    bool device_modified{false};
};

/** Representation of an unstructured mesh
 */
struct MeshBlockInfo
{
    static constexpr int max_vertex_types = 4;
    using GlobalIDType = uint64_t;

    enum QVarType {
        ROW = 0,
        COL
    };

    TiogaView<int> wall_ids;
    TiogaView<int> overset_ids;
    TiogaView<int> num_vert_per_elem;
    TiogaView<int> num_cells_per_elem;
    TiogaView<int> vertex_conn[max_vertex_types];

    TiogaView<double> xyz;
    TiogaView<int> iblank_node;
    TiogaView<int> iblank_cell;
    TiogaView<double> qnode;

    TiogaView<GlobalIDType> cell_gid;
    TiogaView<GlobalIDType> node_gid;

    TiogaView<double> node_res;
    TiogaView<double> cell_res;

    int meshtag{0};
    int num_nodes{0};
    int num_vars{0};

    QVarType qtype{ROW};
};

/** Representation of an AMReX mesh
 */
struct AMRMeshInfo
{
    // Patch info common to all MPI ranks [ngrids_global]
    TiogaView<int> level;
    TiogaView<int> mpi_rank;
    TiogaView<int> local_id;
    TiogaView<int> ilow;
    TiogaView<int> ihigh;
    TiogaView<double> xlo;
    TiogaView<double> dx;

    // Patch data for all patches belonging to this MPI rank
    // [ngrids_local]
    TiogaView<int> global_idmap;
    TiogaView<int> iblank_node;
    TiogaView<int> iblank_cell;
    TiogaView<double> qcell;
    TiogaView<double> qnode;

    int num_ghost{0};
};

} // namespace TIOGA

#endif /* TIOGAMESHINFO_H */
