#define NQ 5
void get_amr_index_xyz(int i,int j,int k,
		   int pBasis,
		   int nX,int nY,int nZ,
		   int nf,
		   double *xlo,double *dx,
		   double *qnodes,
		   int* index, double* xyz)
{
  int ind = 0;
  int base;
  int two_nf = 2*nf;
  int starting_entry;
  int qp_1d = pBasis + 1;  
        
  int grid_stride = NQ*(nX + two_nf)*(nY + two_nf)*(nZ + two_nf);
  int z_stride = grid_stride*(qp_1d*qp_1d);
  int y_stride = grid_stride*(qp_1d);
        
  double cell_xlo = xlo[0] + i*dx[0];
  double cell_ylo = xlo[1] + j*dx[1];
  double cell_zlo = xlo[2] + k*dx[2];
        
  double half_dx = 0.5*dx[0];
  double half_dy = 0.5*dx[1];
  double half_dz = 0.5*dx[2];
  
  int z,y,x;
  
      
  //Q[nX+2*nf,nY+2*nf,nZ+2*nf,nq,p+1,p+1,p+1]
  starting_entry = (nY + two_nf)*(nX + two_nf)*(k+nf)
    + (nX + two_nf)*(j+nf)
    + (i + nf);
        
  for (z = 0;z < qp_1d;z++){
    for (y = 0;y < qp_1d;y++){
      for (x = 0;x < qp_1d;x++){                    
                    index[ind] = starting_entry
		      + z_stride*(z)
		      + y_stride*(y)
		      + grid_stride*(x);
                                                            
                    base = 3*ind;
                    xyz[base]   = cell_xlo + half_dx*(1.0 + qnodes[x]);
                    xyz[base+1] = cell_ylo + half_dy*(1.0 + qnodes[y]);
                    xyz[base+2] = cell_zlo + half_dz*(1.0 + qnodes[z]);
                    
                    ind++;
      }
    }
  }
}
