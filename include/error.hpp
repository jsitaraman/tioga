#pragma once

#include <cstdlib>
#include <execinfo.h>
#include <stdio.h>
#include <iostream>
#include <unistd.h>

/* Nice NVTX macro from Parallel Forall blog */
#if defined(_NVTX)
#include "cuda_runtime.h"
#include "nvToolsExt.h"

const uint32_t colors[] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define PUSH_NVTX_RANGE(name,cid) { \
      int color_id = cid; \
      color_id = color_id%num_colors;\
      nvtxEventAttributes_t eventAttrib = {0}; \
      eventAttrib.version = NVTX_VERSION; \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
      eventAttrib.colorType = NVTX_COLOR_ARGB; \
      eventAttrib.color = colors[color_id]; \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
      eventAttrib.message.ascii = name; \
      nvtxRangePushEx(&eventAttrib); \
}
#define POP_NVTX_RANGE {nvtxRangePop();}
#else
#define PUSH_NVTX_RANGE(name,cid)
#define POP_NVTX_RANGE
#endif

#ifdef _GPU
#define check_error() \
{ \
  cudaError_t err = cudaGetLastError(); \
  if (err != cudaSuccess) \
  { \
    std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": " << std::endl; \
    FatalErrorST(cudaGetErrorString(err)); \
  } \
}
#else
#define check_error()
#endif

//! Prints the error message, the source file and line number, and exits
#define FatalError(s) {                                             \
  printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        \
  exit(1); }

//! Prints the error message, the source file ane line number, the full stack trace, and exits
#define FatalErrorST(s) {                                           \
  void* array[10];                                                  \
  size_t size;                                                      \
  size = backtrace(array,10);                                       \
  printf("Fatal error '%s' at %s:%d\n\n",s,__FILE__,__LINE__);      \
  backtrace_symbols_fd(array, size, STDERR_FILENO);                 \
  exit(1); }


#define _(x) cout << #x << ": " << x << endl;
#define _print(x,y) cout << #x << ": " << x << ", " << #y << ": " << y << endl;
