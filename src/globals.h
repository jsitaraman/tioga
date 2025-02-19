// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

/**
 global definitions for
 the interface code
*/

#ifndef GLOBALS_H
#define GLOBALS_H

#include "tioga.h"
using namespace TIOGA;
#define MAXBLOCKS 100
tioga *tg;
/*
** pointer storage for connectivity arrays that
** comes from external solver
*/
typedef struct inpdata
{
int **vconn;
int *nc;
int *nv;
} inpdata;
inpdata idata[MAXBLOCKS];

#endif /* GLOBALS_H */
