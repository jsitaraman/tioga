// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include<stdio.h>
#include<stdlib.h>
#define REAL double
#include "inputBlk.h"

void MeshBlock::writeOutput(int bid)
{
  char fname[80];
  char hash,c;
  int i,j;
  int bodytag;
  FILE *fp;

  sprintf(fname,"cell_dcf%d.plt",bid);
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"DCF output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",blk[bid].nnodes,
	  blk[bid].ncells);
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");
  for(i=0;i<blk[bid].nnodes;i++)
    fprintf(fp,"%.14e\n",blk[bid].x[3*i]);
  for(i=0;i<blk[bid].nnodes;i++)
    fprintf(fp,"%.14e\n",blk[bid].x[3*i+1]);
  for(i=0;i<blk[bid].nnodes;i++)
    fprintf(fp,"%.14e\n",blk[bid].x[3*i+2]);
  for(i=0;i<blk[bid].ncells;i++)
    fprintf(fp,"%d\n",blk[bid].iblank[i]);

  if (blk[bid].ntetra > 0) {
    for(i=0;i<blk[bid].ntetra;i++)
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
	      blk[bid].ndc4[4*i]+1,
	      blk[bid].ndc4[4*i+1]+1,
	      blk[bid].ndc4[4*i+2]+1,
	      blk[bid].ndc4[4*i+2]+1,
	      blk[bid].ndc4[4*i+3]+1,
	      blk[bid].ndc4[4*i+3]+1,
	      blk[bid].ndc4[4*i+3]+1,
	      blk[bid].ndc4[4*i+3]+1);
  }

  if (blk[bid].npyra > 0) {
    for(i=0;i<blk[bid].ntetra;i++)
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
	      blk[bid].ndc5[5*i]+1,
	      blk[bid].ndc5[5*i+1]+1,
	      blk[bid].ndc5[5*i+2]+1,
	      blk[bid].ndc5[5*i+3]+1,
	      blk[bid].ndc4[5*i+4]+1,
	      blk[bid].ndc4[5*i+4]+1,
	      blk[bid].ndc4[5*i+4]+1,
	      blk[bid].ndc4[5*i+4]+1);
  }

  if (blk[bid].nprizm > 0) {
    for(i=0;i<blk[bid].nprizm;i++)
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
	      blk[bid].ndc6[6*i]+1,
	      blk[bid].ndc6[6*i+1]+1,
	      blk[bid].ndc6[6*i+2]+1,
	      blk[bid].ndc6[6*i+2]+1,
	      blk[bid].ndc6[6*i+3]+1,
	      blk[bid].ndc6[6*i+4]+1,
	      blk[bid].ndc6[6*i+5]+1,
	      blk[bid].ndc6[6*i+5]+1);
  }
  
	      
  if (blk[bid].nhexa > 0) {
    for(i=0;i<blk[bid].nhexa;i++)
      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
	      blk[bid].ndc8[8*i]+1,
	      blk[bid].ndc8[8*i+1]+1,
	      blk[bid].ndc8[8*i+2]+1,
	      blk[bid].ndc8[8*i+3]+1,
	      blk[bid].ndc8[8*i+4]+1,
	      blk[bid].ndc8[8*i+5]+1,
	      blk[bid].ndc8[8*i+6]+1,
	      blk[bid].ndc8[8*i+7]+1);
  }
	      
  fclose(fp);

  return;
}
  
  
  
  
