// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

/* Fallback for builds configured without TIOGA_ENABLE_CUDA. Keeps the
   MeshBlock interface identical so callers can probe at run time. */

#include "codetypes.h"
#include "MeshBlock.h"

int MeshBlock::search_gpu(void)
{
  static int warned = 0;
  if (!warned) {
    fprintf(stderr,
            "#tioga: search_gpu() unavailable, "
            "rebuild with -DTIOGA_ENABLE_CUDA=ON\n");
    warned = 1;
  }
  return -1;
}

void MeshBlock::freeGpuSearchData(void)
{
}

int MeshBlock::search_gpu_batch(MeshBlock ** /*blocks*/, int /*nblocks*/,
                                SEARCHTIMERS * /*timers*/)
{
  fprintf(stderr,
          "#tioga: search_gpu_batch() unavailable, "
          "rebuild with -DTIOGA_ENABLE_CUDA=ON\n");
  return -1;
}

void MeshBlock::freeGpuBatchData(void)
{
}
