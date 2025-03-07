// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#include "codetypes.h"

void deallocateLinkList(DONORLIST *temp)
{
  if (temp!=NULL) 
    {
      deallocateLinkList(temp->next);
      TIOGA_FREE(temp);
    }
}

void deallocateLinkList2(INTEGERLIST *temp)
{
  if (temp!=NULL) 
    {
      deallocateLinkList2(temp->next);
      TIOGA_FREE(temp);
    }
}

void deallocateLinkList3(INTEGERLIST2 *temp)
{
  if (temp!=NULL)
    {
      if (temp->intData) TIOGA_FREE(temp->intData);
      if (temp->realData) TIOGA_FREE(temp->realData);
      deallocateLinkList3(temp->next);
      TIOGA_FREE(temp);
    }
}

void deallocateLinkList4(INTERPLIST2 *temp)
{
  if (temp!=NULL)
    {
      if (temp->inode) TIOGA_FREE(temp->inode);
      if (temp->weights) TIOGA_FREE(temp->weights);
      deallocateLinkList4(temp->next);
      TIOGA_FREE(temp);
    }
}

void insertInList(DONORLIST **donorList,DONORLIST *temp1)
{
  DONORLIST *temp;
  DONORLIST *ptemp;
  int inserted;  
  temp=*donorList;
  inserted=0;
  ptemp=NULL;
  while(temp!=NULL && !inserted)
    {
      if (fabs(temp->donorRes) > temp1->donorRes) 
	{
	  temp1->next=temp;
	  if (ptemp!=NULL) 
	    {
	      ptemp->next=temp1;
	    }
	  else
	    {
	      *donorList=temp1;
	    }	  
	  inserted=1;
	}
      else
	{
	  ptemp=temp;
	  temp=temp->next;
	}
    }
  if (!inserted) 
    {
     if (*donorList) 
      {
       temp1->next=NULL;
       ptemp->next=temp1;
      }
     else
       {
        temp1->next=NULL;
        *donorList=temp1;
       }
    }
}

