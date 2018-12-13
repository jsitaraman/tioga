//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef PARALLELCOMM_H
#define PARALLELCOMM_H
#include "codetypes.h"
#include <cstdlib>
#include "mpi.h"

struct PACKET;

/**
* Parallel communication methods
* the MPI calls are abstracted into
* these methods and they provide the
* ability to transfer packets across
* processors */
class parallelComm
{
 private:
  int nsend;
  int nrecv;
  int *sndMap;
  int *rcvMap;

 public :
  int myid;
  int numprocs;
  MPI_Comm scomm;
  
  parallelComm() { sndMap=NULL; rcvMap=NULL;}
  
 ~parallelComm() { if (sndMap) free(sndMap);
                   if (rcvMap) free(rcvMap);}

  void sendRecvPacketsAll(PACKET *sndPack,PACKET *rcvPack);
  
  void sendRecvPackets(PACKET *sndPack,PACKET *rcvPack);

  void sendRecvPacketsCheck(PACKET *sndPack,PACKET *rcvPack);

  void setMap(int ns, int nr, int *snd,int *rcv);

  void getMap(int *ns,int *nr, int **snd, int **rcv);

  void initPackets(PACKET *sndPack, PACKET *rcvPack);

  void clearPackets(PACKET *sndPack, PACKET *rcvPack);
  void clearPackets2(PACKET *sndPack, PACKET *rcvPack);
  
};
  



#endif /* PARALLELCOMM_H */
