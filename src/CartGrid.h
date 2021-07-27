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

#ifndef CARTGRID_H
#define CARTGRID_H

#include <cstdlib>

namespace TIOGA
{
struct AMRMeshInfo;
}

class CartGrid
{
 private:
  double xlosup[3];
  double *dxlvl{nullptr};
  int *lcount{nullptr};
  int maxlevel;

  bool own_data_ptrs{true};
  bool own_amr_mesh_info{false};

 public :
  TIOGA::AMRMeshInfo* m_info{nullptr};

  int *global_id{nullptr};
  int *level_num{nullptr};
  int *proc_id{nullptr};
  int *local_id{nullptr};
  int *ilo{nullptr};
  int *ihi{nullptr};
  int *dims{nullptr};
  int myid{0};
  int nf{0};
  double *xlo{nullptr};
  double *dx{nullptr};
  int ngrids{0};
  void (*donor_frac) (int *,double *,int *,double *) = nullptr;

  CartGrid() = default;
  ~CartGrid();

  void registerData(TIOGA::AMRMeshInfo* m_info);
  void registerData(int nf,int *idata,double *rdata,int ngridsin);
  void preprocess(void);     
  void search(double *x,int *donorid,int nsearch);
  void setcallback(void (*f1)(int *,double *,int *,double *)) 
  {   donor_frac=f1;  }

  void create_mesh_info();
};

#endif /* CARTGRID_H */
