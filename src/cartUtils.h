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
#ifndef CARTUTILS_H
#define CARTUTILS_H

namespace cart_utils
{
  //Q[nq,nZ+2*nf,nY+2*nf,nX+2*nf]--> C++ Cell storage
  inline int get_cell_index(int nX,int nY,int nf,int i,int j,int k)
  {
    return (nY+2*nf)*(nX+2*nf)*(k+nf) + (nX+2*nf)*(j+nf) + (i+nf);
  }

  //Q[nq,nZ+1+2*nf,nY+1+2*nf,nX+1+2*nf]--> C++ node storage
  inline int get_node_index(int nX,int nY,int nf,int i,int j,int k)
  {
    return (nY+1+2*nf)*(nX+1+2*nf)*(k+nf) + (nX+1+2*nf)*(j+nf) + (i+nf);
  }

  //Q[nq,nZ+1+2*nf,nY+1+2*nf,nX+1+2*nf]--> C++ node storage
  // for arrays with both node and cell indices,
  // node index is assumed to follow cell index
  inline int get_concatenated_node_index(int nX,int nY,int nZ,int nf,int i,int j,int k)
  {
    return (nY+1+2*nf)*(nX+1+2*nf)*(k+nf) + (nX+1+2*nf)*(j+nf) + (i+nf)
        + (nX+2*nf)*(nY+2*nf)*(nZ+2*nf);
  }
} // namespacee cart_utils

#endif /* CARTUTILS_H */
