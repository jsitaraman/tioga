void get_amr_index_xyz( int nq,int i,int j,int k,
			int pBasis,
			int nX,int nY,int nZ,
			int nf,
			double *xlo,double *dx,
			double qnodes[],
			int* index, double* xyz){
        
  int base;
  int starting_entry;
        
  int ind = 0;               
  int two_nf = 2*nf;
  int qp_1d = pBasis + 1;  
        
  int grid_stride = nq*(nX + two_nf)*(nY + two_nf)*(nZ + two_nf);
  int z_stride = grid_stride*(qp_1d*qp_1d);
  int y_stride = grid_stride*(qp_1d);
        
  double cell_xlo = xlo[0] + i*dx[0];
  double cell_ylo = xlo[1] + j*dx[1];
  double cell_zlo = xlo[2] + k*dx[2];
        
  double half_dx = 0.5*dx[0];
  double half_dy = 0.5*dx[1];
  double half_dz = 0.5*dx[2];
        
  double yco = 0.0;
  double zco = 0.0;
  int x,y,z;
                
  //Q[p+1,p+1,p+1,nq,nZ+2*nf,nY+2*nf,nX+2*nf]--> C++ Storage
  starting_entry = (nY+two_nf)*(nX+two_nf)*(k+nf)
    + (nX+two_nf)*(j+nf)
    + (i+nf);
        
  for (z = 0;z < qp_1d;z++){
    zco = cell_zlo + half_dz*(1.0 + qnodes[z]);
    for (y = 0;y < qp_1d;y++){
      yco = cell_ylo + half_dy*(1.0 + qnodes[y]);
      for ( x = 0;x < qp_1d;x++){                    
                    index[ind] = starting_entry
		      + z_stride*(z)
		      + y_stride*(y)
		      + grid_stride*(x);
                                                            
                    base = 3*ind;
                    xyz[base]   = cell_xlo + half_dx*(1.0 + qnodes[x]);
                    xyz[base+1] = yco;
                    xyz[base+2] = zco;
                    
                    ind++;
      }
    }
  }
}

void amr_index_to_ijklmn(int pBasis,int nX,int nY,int nZ, int nf, int nq,
			 int index, int* ijklmn)
{
        
  //Q[p+1,p+1,p+1,nq,nZ+2*nf,nY+2*nf,nX+2*nf]--> C++ Storage

  int remainder;

  int two_nf = 2*nf;
  int qp_1d = pBasis + 1;  

  int x_size = (nX + two_nf);
  int y_size = (nY + two_nf);
  int z_size = (nZ + two_nf);

  int xy_stride = x_size*y_size;
  int xyz_stride = xy_stride*z_size;

  int grid_stride = nq*x_size*y_size*z_size;
  int x_stride = grid_stride;
  int y_stride = grid_stride*(qp_1d);
  int z_stride = grid_stride*(qp_1d*qp_1d);

  remainder = index;
  int n = index / z_stride;
  remainder = index % z_stride;

  int m = remainder / y_stride;
  remainder %= y_stride;                

  int l = remainder / x_stride;
  remainder %= x_stride;

  //Skip fields
  remainder %= xyz_stride;

  int k = remainder / xy_stride;
  remainder %= xy_stride;

  int j = remainder / x_size;
  remainder %= x_size;

  int i = remainder;

  ijklmn[0] = i-nf;
  ijklmn[1] = j-nf;
  ijklmn[2] = k-nf;
  ijklmn[3] = l;
  ijklmn[4] = m;
  ijklmn[5] = n;    
}
