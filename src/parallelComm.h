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
#include "mpi.h"
#include <vector>

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

  /// TESTING SEPARATE SEND/RECV FOR OVERLAP
  std::vector<int> scount, rcount;
  std::vector<MPI_Request> reqs;
  std::vector<MPI_Status> stats;
  int nwait;
  
  parallelComm() { sndMap=NULL; rcvMap=NULL;}
  
 ~parallelComm() { if (sndMap) free(sndMap);
                   if (rcvMap) free(rcvMap);}

  void sendRecvPacketsAll(PACKET *sndPack,PACKET *rcvPack);
  
  void sendRecvPackets(PACKET *sndPack,PACKET *rcvPack);

  void sendRecvPacketsV(std::vector<VPACKET> &sndPack, std::vector<VPACKET> &rcvPack);

  void sendPacketsV(std::vector<VPACKET> &sndPack, std::vector<VPACKET> &rcvPack);
  void sendPacketsV2(std::vector<PACKET> &sndPack, std::vector<VPACKET> &rcvPack);
  void recvPacketsV(void);

  void sendRecvPacketsCheck(PACKET *sndPack,PACKET *rcvPack);

  void setMap(int ns, int nr, int *snd,int *rcv);

  void getMap(int *ns,int *nr, int **snd, int **rcv);

  void initPackets(PACKET *sndPack, PACKET *rcvPack);
  void initPacketsV(std::vector<VPACKET>& sndPack, std::vector<VPACKET>& rcvPack);
  void initPacketsV2(std::vector<PACKET>& sndPack, std::vector<VPACKET>& rcvPack);

  void clearPackets(PACKET *sndPack, PACKET *rcvPack);
  void clearPackets2(PACKET *sndPack, PACKET *rcvPack);
  void clearPacketsV(std::vector<VPACKET> &sndPack, std::vector<VPACKET> &rcvPack);
  void clearPacketsV2(std::vector<PACKET> &sndPack, std::vector<VPACKET> &rcvPack);
  
};
  
