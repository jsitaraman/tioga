TEMPLATE = library
CONFIG -= app_bundle
CONFIG -= qt
DEFINES = _GPU __CUDACC__

INCLUDEPATH += include/

SOURCES += src/tioga.C \
    src/ADT.C \
    src/bookKeeping.C \
    src/CartBlock.C \
    src/CartGrid.C \
    src/cartOps.C \
    src/checkContainment.C \
    src/dataUpdate.C \
    src/exchangeAMRDonors.C \
    src/exchangeBoxes.C \
    src/exchangeDonors.C \
    src/exchangeSearchData.C \
    src/buildADTrecursion.c \
    src/computeCellVolume.c \
    src/funcs.cpp \
    src/getCartReceptors.C \
    src/highOrder.C \
    src/holeMap.C \
    src/linklist.C \
    src/math.C \
    src/MeshBlock.C \
    src/parallelComm.C \
    src/points.cpp \
    src/search.C \
    src/searchADTrecursion.C \
    src/tiogaInterface.C \
    src/utils.C \
    src/writeOutput.C \
    src/get_amr_index_xyz.c \
    src/test_amr_functions.c \
    src/cellVolume.f90 \
    src/highOrder_kernels.cu \
    src/kaiser.f \
    src/median.F90 \
    src/writeqnode.f90 \
    src/dADT.cu \
    src/dMeshBlock.cu \
    src/superMesh.cpp

HEADERS += \
    include/ADT.h \
    include/CartBlock.h \
    include/CartGrid.h \
    include/codetypes.h \
    include/cuda_funcs.h \
    include/error.hpp \
    include/funcs.hpp \
    include/globals.h \
    include/highOrder_kernels.h \
    include/MeshBlock.h \
    include/parallelComm.h \
    include/points.hpp \
    include/tioga.h \
    include/tioga.i \
    include/tiogaInterface.h \
    include/dADT.h \
    include/dMeshBlock.h \
    include/superMesh.hpp \
    include/convert.i

DISTFILES += \
    run/tiogaInterface.py
