#include "dMeshBlock.h"
#include "funcs.hpp"

/* ------ dMeshBlock Member Functions ------ */

void dMeshBlock::dataToDevice(int ndims, int nnodes, int ncells, int ncells_adt,
    int nsearch, int* nv, int* nc, int* eleList, double* eleBBox, int* isearch,
    double* xsearch)
{
  this->nnodes = nnodes;
  this->ncells = ncells;
  this->nc_adt = ncells_adt;

  this->nv = nv;
  this->nc = nc;

  nvert = nv[0];

  this->eleBBox.assign(eleBBox, ncells_adt*ndims*2);
  this->eleList.assign(eleList, ncells_adt);

  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);

  auto ijk2gmsh_h = tg_funcs::structured_to_gmsh_hex(nvert);
  ijk2gmsh.assign(ijk2gmsh_h.data(), ijk2gmsh_h.size());

  int nSide = std::cbrt(nvert);
  std::vector<double> xlist_h(nSide);
  double dxi = 2./(nSide-1);

  for (int i = 0; i < nSide; i++)
    xlist_h[i] = -1. + i*dxi;

  xlist.assign(xlist_h.data(), xlist_h.size());
}

void dMeshBlock::extraDataToDevice(int* vconn)
{
//  c2v.assign(vconn, nvert*ncells);
}

void dMeshBlock::updateSearchPoints(int nsearch, int *isearch, double *xsearch)
{
  this->nsearch = nsearch;
  this->isearch.assign(isearch, nsearch);
  this->xsearch.assign(xsearch, nsearch*nDims);
  rst.resize(nsearch*nDims);
  donorId.resize(nsearch);
}

void dMeshBlock::setDeviceData(double* vx, double* ex, int* ibc, int* ibf)
{
  x = vx;
  iblank_cell = ibc;
  iblank_face = ibf;
  coord = ex;
}

void dMeshBlock::setTransform(double* mat, double* off, int ndim)
{
  if (ndim != nDims)
    ThrowException("dMeshBlock::set_transform: input ndim != nDims");

  rrot = true;
  Rmat.assign(mat, ndim*ndim);
  offset.assign(off, ndim);
}

/* ---------------------------- Direct Cut Method Functions --------------------------- */

static
__device__
double lineSegmentDistance(double *p1, double *p2, double *p3, double *p4, double *dx)
{
  // Get the line equations
  double U[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
  double V[3] = {p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2]};
  double W[3] = {p1[0]-p3[0], p1[1]-p3[1], p1[2]-p3[2]};
  double uu = U[0]*U[0] + U[1]*U[1] + U[2]*U[2];
  double vv = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
  double uv = U[0]*V[0] + U[1]*V[1] + U[2]*V[2];

  double uw = U[0]*W[0] + U[1]*W[1] + U[2]*W[2];
  double vw = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];

  double den = uu*vv - uv*uv;

  // NOTE: not finding exact minimum distance between the line segments in all
  // cases; plenty close enough for our purposes
  // (see http://geomalgorithms.com/a07-_distance.html for full algo)

  // Calculate line parameters (if nearly parallel, set one & calculate other)
  double s = (den < 1e-10) ? 0 : (uv*vw - vv*uw) / den;
  double t = (den < 1e-10) ? uw / uv: (uu*vw - uv*uw) / den;

  s = fmin(max(s, 0.), 1.);
  t = fmin(max(t, 0.), 1.);

  // vec = closest distance from segment 1 to segment 2
  for (int i = 0; i < 3; i++)
    dx[i] = (p3[i] + t*V[i]) - (p1[i] + s*U[i]);

  double dist = 0;
  for (int i = 0; i < 3; i++)
    dist += dx[i]*dx[i];

  return sqrt(dist);
}

/*! Modified Moller triangle-triangle intersection algorithm
 *  Determines if triangles intersect, or returns an approximate minimum
 *  distance between them otherwise */
static
__device__
double triTriDistance(double* T1, double* T2, double* minVec, double tol)
{
  double dist = 1e15;
  double vec[3];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      int i2 = (i+1) % 3;
      int j2 = (j+1) % 3;
      double D = lineSegmentDistance(&T1[3*i], &T1[3*i2], &T2[3*j], &T2[3*j2], vec);

      if (D < dist)
      {
        for (int d = 0; d < 3; d++)
          minVec[d] = vec[d];
        dist = D;
      }
    }
  }

  if (dist < tol) return 0.;

  /// TODO: optimize (no use of point class)
  // Compute pi2 - N2 * X + d2 = 0
  dPoint V01 = dPoint(T1);
  dPoint V11 = dPoint(T1+3);
  dPoint V21 = dPoint(T1+6);

  dPoint V02 = dPoint(T2);
  dPoint V12 = dPoint(T2+3);
  dPoint V22 = dPoint(T2+6);

  dPoint N2 = (V12-V02).cross(V22-V02);
  N2 /= N2.norm();
  double d2 = -(N2*V02);

  dPoint N1 = (V11-V01).cross(V21-V01);
  N1 /= N1.norm();
  double d1 = -(N1*V01);

  // Signed distances of T1 points to T2's plane
  double d01 = N2*V01 + d2;
  double d11 = N2*V11 + d2;
  double d21 = N2*V21 + d2;

  // Signed distances of T2 points to T1's plane
  double d02 = N1*V02 + d1;
  double d12 = N1*V12 + d1;
  double d22 = N1*V22 + d1;

  // Round values near 0 to 0
  d01 = (fabs(d01) < 1e-10) ? 0 : d01;
  d11 = (fabs(d11) < 1e-10) ? 0 : d11;
  d21 = (fabs(d21) < 1e-10) ? 0 : d21;

  d02 = (fabs(d02) < 1e-10) ? 0 : d02;
  d12 = (fabs(d12) < 1e-10) ? 0 : d12;
  d22 = (fabs(d22) < 1e-10) ? 0 : d22;

  if (fabs(d01) + fabs(d11) + fabs(d21) < 3*tol ||
      fabs(d02) + fabs(d12) + fabs(d22) < 3*tol)
  {
    /* Approximately coplanar; check if one triangle is inside the other */

    // Check if a point in T1 is inside T2
    bool inside = true;
    inside = inside && N2*((V12-V02).cross(V01-V02)) > 0;
    inside = inside && N2*((V02-V22).cross(V01-V22)) > 0;
    inside = inside && N2*((V22-V12).cross(V01-V12)) > 0;

    if (inside) return 0.;

    // Check if a point in T2 is inside T1
    inside = true;
    inside = inside && N1*((V11-V01).cross(V02-V01)) > 0;
    inside = inside && N1*((V01-V21).cross(V02-V21)) > 0;
    inside = inside && N1*((V21-V11).cross(V02-V11)) > 0;

    if (inside) return 0.;
  }

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d01) == sgn(d11) && sgn(d01) == sgn(d21)) // && fabs(d01) > tol)
  {
    // No intersection; return result from edge distances
    return dist;
  }

  // Check for intersection with plane - one point should have opposite sign
  if (sgn(d02) == sgn(d12) && sgn(d02) == sgn(d22)) // && fabs(d02) > tol)
    return dist;

  // Compute intersection line
  dPoint L = N1.cross(N2);
  L /= L.norm();

  double p0 = L*V01;
  double p1 = L*V11;
  double p2 = L*V21;

  double q0 = L*V02;
  double q1 = L*V12;
  double q2 = L*V22;

  // Figure out which point of each triangle is opposite the other two
  int npt1 = (sgn(d01) != sgn(d11)) ? ( (sgn(d11) == sgn(d21)) ? 0 : 1 ) : 2;
  int npt2 = (sgn(d02) != sgn(d12)) ? ( (sgn(d12) == sgn(d22)) ? 0 : 1 ) : 2;

  double s1, s2;
  switch (npt1)
  {
    case 0:
      s1 = p1 + (p0-p1) * (d11 / (d11-d01));
      s2 = p2 + (p0-p2) * (d21 / (d21-d01));
      break;
    case 1:
      s1 = p0 + (p1-p0) * (d01 / (d01-d11));
      s2 = p2 + (p1-p2) * (d21 / (d21-d11));
      break;
    case 2:
      s1 = p0 + (p2-p0) * (d01 / (d01-d21));
      s2 = p1 + (p2-p1) * (d11 / (d11-d21));
      break;
  }

  double t1, t2;
  switch (npt2)
  {
    case 0:
      t1 = q1 + (q0-q1) * (d12 / (d12-d02));
      t2 = q2 + (q0-q2) * (d22 / (d22-d02));
      break;
    case 1:
      t1 = q0 + (q1-q0) * (d02 / (d02-d12));
      t2 = q2 + (q1-q2) * (d22 / (d22-d12));
      break;
    case 2:
      t1 = q0 + (q2-q0) * (d02 / (d02-d22));
      t2 = q1 + (q2-q1) * (d12 / (d12-d22));
      break;
  }

  if (s1 > s2)
    swap(s1,s2);

  if (t1 > t2)
    swap(t1,t2);

  if (s2 < t1 || t2 < s1)
  {
    // No overlap; return min of dt*L and minDist
    double dt = fmin(fabs(t1-s2), fabs(s1-t2));
    double dl = (L*dt).norm();

    if (dl < dist)
    {
      dist = dl;
      for (int i = 0; i < 3; i++)
        minVec[i] = sgn(t1-s2)*dt*L[i]; // Ensure vec is T1 -> T2
    }

    return dist;
  }

  return 0.;
}

static
__device__
dPoint faceNormal(const double* xv)
{
  /* Assuming nodes of face ordered CCW such that right-hand rule gives
     * outward normal */

  // Triangle #1
  dPoint pt0 = dPoint(&xv[0]);
  dPoint pt1 = dPoint(&xv[3]);
  dPoint pt2 = dPoint(&xv[6]);
  dPoint a = pt1 - pt0;
  dPoint b = pt2 - pt0;
  dPoint norm1 = a.cross(b);           // Face normal vector

  // Triangle #2
  pt1 = dPoint(&xv[9]);
  a = pt1 - pt0;
  dPoint norm2 = b.cross(a);

  // Average the two triangle's normals
  dPoint norm = 0.5*(norm1+norm2);
  norm /= norm.norm();

  return norm;
}

__device__
dPoint intersectionCheck(dMeshBlock &mb, double *fxv, int nfv, double *exv, int nev)
{
  /* --- Prerequisites --- */

  // NOTE: Structured ordering  |  btm,top,left,right,front,back
  short TriPts[12][3] = {{0,1,3},{0,3,2},{4,7,5},{4,6,7},{0,2,6},{0,6,4},
                       {1,3,7},{1,7,5},{0,4,5},{0,5,1},{2,3,7},{2,6,7}};

  double tol = 1e-9;

  int nsideC = cbrt((float)nev);
  int sorderC = nsideC-1;

  int nsideF = sqrt((float)nfv);
  int sorderF = nsideF-1;

  double TC[9], TF[9];
  double minDist = BIG_DOUBLE;
  double minVec[3] = {minDist,minDist,minDist};
  double bboxC[6], bboxF[6];

  cuda_funcs::getBoundingBox<3>(fxv, nfv, bboxF);
  cuda_funcs::getBoundingBox<3>(exv, nev, bboxC);

  /* Only 3 cases possible:
   * 1) Face entirely contained within element
   * 2) Face intersects with element's boundary
   * 3) Face and element do not intersect
   */

  // 1) In case of face entirely inside element, check if a pt is inside ele
  bool inside = false;
  double rst[3];
  switch (nsideC)
  {
    case 2:
      inside = mb.getRefLoc<2>(exv, bboxC, fxv, rst);
      break;
    case 3:
      inside = mb.getRefLoc<3>(exv, bboxC, fxv, rst);
      break;
    case 4:
      inside = mb.getRefLoc<4>(exv, bboxC, fxv, rst);
      break;
    case 5:
      inside = mb.getRefLoc<5>(exv, bboxC, fxv, rst);
      break;
    default:
      printf("nsideC case not implemented for checkPtInEle");
  }

  if (inside)
    return dPoint(0., 0., 0.);

  // 2) Check outer faces of element for intersection with face

  for (int f = 0; f < 6; f++)
  {
    for (int g = 0; g < sorderC*sorderC; g++)
    {
      int I, J, K;
      switch (f)
      {
        case 0: // Bottom
          I = g / sorderC;
          J = g % sorderC;
          K = 0;
          break;
        case 1: // Top
          I = g / sorderC;
          J = g % sorderC;
          K = sorderC - 1;
          break;
        case 2: // Left
          I = 0;
          J = g / sorderC;
          K = g % sorderC;
          break;
        case 3: // Right
          I = sorderC - 1;
          J = g / sorderC;
          K = g % sorderC;
          break;
        case 4: // Front
          I = g / sorderC;
          J = 0;
          K = g % sorderC;
          break;
        case 5: // Back
          I = g / sorderC;
          J = sorderC - 1;
          K = g % sorderC;
          break;
      }

      int i0 = I+nsideC*(J+nsideC*K);
      int j0 = i0 + nsideC*nsideC;
      int lin2curv[8] = {i0, i0+1, i0+nsideC, i0+nsideC+1, j0, j0+1, j0+nsideC, j0+nsideC+1};
      for (int i = 0; i < 8; i++)
        lin2curv[i] = mb.ijk2gmsh[lin2curv[i]]; /// TODO

      // Get triangles for the sub-hex of the larger curved hex
      for (int i = f; i < f+2; i++)
      {
        for (int p = 0; p < 3; p++)
        {
          int ipt = lin2curv[TriPts[i][p]];
          for (int d = 0; d < 3; d++)
            TC[3*p+d] = exv[3*ipt+d];
        }

        cuda_funcs::getBoundingBox<3>(TC, 3, bboxC);
        //double btol = .05*(bboxC[3]-bboxC[0]+bboxC[4]-bboxC[1]+bboxC[5]-bboxC[2]);
        double btol = (bboxC[3]-bboxC[0]+bboxC[4]-bboxC[1]+bboxC[5]-bboxC[2]);
        btol = fmin(btol, minDist);
        if (!cuda_funcs::boundingBoxCheck<3>(bboxC,bboxF,btol)) continue;

        for (int M = 0; M < sorderF; M++)
        {
          for (int N = 0; N < sorderF; N++)
          {
            int m0 = M + nsideF*N;
            int TriPtsF[2][3] = {{m0, m0+1, m0+nsideF+1}, {m0, m0+nsideF+1, m0+nsideF}};

            // Intersection check between element face tris & cutting-face tris
            for (int j = 0; j < 2; j++)
            {
              for (int p = 0; p < 3; p++)
              {
                int ipt = mb.ijk2gmsh_quad[TriPtsF[j][p]];
                for (int d = 0; d < 3; d++)
                  TF[3*p+d] = fxv[3*ipt+d];
              }

              double vec[3];
              double dist = triTriDistance(TF, TC, vec, tol);

              if (dist < tol)
                return dPoint(0.,0.,0);

              if (dist < minDist)
              {
                for (int d = 0; d < 3; d++)
                  minVec[d] = vec[d];
                minDist = dist;
              }
            }
          }
        }
      }
    }
  }

  if (minDist == BIG_DOUBLE) // Definitely no intersection; use centroids to get vector
  {
    cuda_funcs::getCentroid<3>(exv,nev,minVec);
    dPoint vec(minVec);
    cuda_funcs::getCentroid<3>(fxv,nfv,minVec);
    return (vec - dPoint(minVec));
  }

  return dPoint(minVec);
}

template<int nDims, int nvert>
__global__
void fillCutMap(dMeshBlock mb, dvec<double> cutFaces, int nCut, int nvertf,
                dCutMap cutMap, int cutType)
{
  int ic = blockIdx.x * blockDim.x + threadIdx.x;

//  if (ic >= mb.ncells) return;
  if (ic != 0) return;

  double xv[nDims*nvert];
  extern __shared__ double fxv[];

  // Load up the cell nodes into an array
  for (int i = 0; i < nvert; i++)
    for (int d = 0; d < nDims; d++)
      xv[nDims*i+d] = mb.coord[nDims*(i+nvert*ic)+d];

  double bboxC[2*nDims], bboxF[2*nDims];

  cuda_funcs::getBoundingBox<nDims>(xv, nvert, bboxC);

  // btol == 10 times the average side length of the cell's bounding box
  double btol = (bboxC[3]-bboxC[0]+bboxC[4]-bboxC[1]+bboxC[5]-bboxC[2]) * 3.3;

  int stride = nDims*nvertf;

  for (int ff = 0; ff < nCut; ff++)
  {
    for (int i = threadIdx.x; i < nvertf*nDims; i += blockDim.x)
    {
      fxv[i] = cutFaces[ff*stride+i];
    }

    __syncthreads();

    /// TODO: read cutface into shared memory (cuz why not? They're all in sync!)
    /*if (mb.rrot)
      getBoundingBox(&cutFaces[ff*stride], nvertf, nDims, bbox, Rmat.data());
    else*/
    cuda_funcs::getBoundingBox<nDims>(fxv, nvertf, bboxF);

    if (cuda_funcs::boundingBoxCheck<nDims>(bboxC, bboxF, btol))
    {
      // Find distance from face to cell
      dPoint vec = intersectionCheck(mb, fxv, nvertf, xv, nvert);
      double dist = vec.norm();
      vec /= dist;
      printf("(%d, %d): dist = %f\n",ic,ff,dist);
      continue;
      if (dist == 0.) // They intersect
      {
        cutMap.flag[ic] = DC_CUT;
        cutMap.dist[ic] = 0.;
      }
      else if (cutMap.flag[ic] == DC_UNASSIGNED) // || dist < (cutMap.dist[ic] - tol))
      {
        // Unflagged cell, or have a closer face to use
        dPoint norm = faceNormal(fxv);
        if (cutType == 0) norm *= -1;

        double dot = norm*vec;

        cutMap.dist[ic] = dist;
        for (int d = 0; d < 3; d++)
          cutMap.norm[3*ic+d] = norm[d];
        cutMap.dot[ic] = dot;
        cutMap.nMin[ic] = 1;

        if (dot < 0) /// TODO: decide on standard orientation
          cutMap.flag[ic] = DC_HOLE; // outwards normal = inside hole
        else
          cutMap.flag[ic] = DC_NORMAL;
      }
      else if (dist < 1.25*cutMap.dist[ic]) //if (fabs(dist-cutMap.dist[ic]) < 2e-3)
      {
        // Approx. same dist. to two faces; avg. their normals to decide
        dPoint norm = faceNormal(fxv);
        if (cutType == 0) norm *= -1;

        // First, check to see which normal is more direct
        // i.e. try to ignore faces around corners from the ele
        double dot = norm*vec;
        double adot = fabs(dot);

        double dot0 = cutMap.dot[ic];
        double adot0 = fabs(dot0);

        if (dist < .9*cutMap.dist[ic] && adot > adot0)
        {
          // Clearly better; use it
          cutMap.dist[ic] = dist;
          for (int d = 0; d < 3; d++)
            cutMap.norm[3*ic+d] = norm[d];
          cutMap.dot[ic] = dot;
          cutMap.nMin[ic] = 1;

          if (dot < 0) /// TODO: decide on standard orientation
            cutMap.flag[ic] = DC_HOLE; // outwards normal = inside hole
          else
            cutMap.flag[ic] = DC_NORMAL;
        }
        else if (fabs(dist = cutMap.dist[ic]) < 1e-6)
        {
          // They're not that far apart; average them
          cutMap.dist[ic] = fmin(cutMap.dist[ic], dist);
          int N = cutMap.nMin[ic];
          for (int d = 0; d < nDims; d++)
            norm[d] = (N*cutMap.norm[3*ic+d] + norm[d]) / (N + 1.);
          norm /= norm.norm();

          for (int d = 0; d < 3; d++)
            cutMap.norm[3*ic+d] = norm[d];
          cutMap.dot[ic] = (adot0 > adot) ? dot0 : dot;
          cutMap.nMin[ic]++;

          if (cutMap.dot[ic] < 0)
            cutMap.flag[ic] = DC_HOLE;
          else
            cutMap.flag[ic] = DC_NORMAL;
        }
        else if (.866*adot > adot0)
        {
          // This one's at least 30deg closer; use it instead
          for (int d = 0; d < 3; d++)
            cutMap.norm[3*ic+d] = norm[d];
          cutMap.dist[ic] = dist;
          cutMap.dot[ic] = dot;

          if (dot < 0)
            cutMap.flag[ic] = DC_HOLE;
          else
            cutMap.flag[ic] = DC_NORMAL;
        }
        else if (adot >= .866*adot0)
        {
          // They're not that far apart; average them
          int N = cutMap.nMin[ic];
          for (int i = 0; i < 3; i++)
            norm[i] = cutMap.norm[3*ic+i]*adot0*N + norm[i]*adot;
//          for (int d = 0; d < nDims; d++)
//            norm[d] = (N*cutMap.norm[ic][d] + norm[d]) / (N + 1.);
          norm /= norm.norm();

          cutMap.dist[ic] = fmin(cutMap.dist[ic], dist);
          for (int d = 0; d < 3; d++)
            cutMap.norm[3*ic+d] = norm[d];
          cutMap.dot[ic] = (adot0 > adot) ? dot0 : dot;
          cutMap.nMin[ic]++;

          if (cutMap.dot[ic] < 0)
            cutMap.flag[ic] = DC_HOLE;
          else
            cutMap.flag[ic] = DC_NORMAL;
        }
        else
        {
          // Face at least 30deg out of alignment; ignore it (probably perpendicular)
        }
      }
    }
  }
}

//void dMeshBlock::directCut(dvec<double> &cutFaces, int nCut, int nvertf, dCutMap &cutMap, int cutType)
void dMeshBlock::directCut(double* cutFaces_h, int nCut, int nvertf, int* cutFlag, int cutType)
{
  // Setup cutMap TODO: create initialization elsewhere?
  dCutMap cutMap;
  cutMap.flag.resize(ncells);
  cutMap.dist.resize(ncells);
  cutMap.nMin.resize(ncells);
  cutMap.norm.resize(nDims*ncells);
  cutMap.dot.resize(ncells);

  dvec<double> cutFaces;
  cutFaces.assign(cutFaces_h, nCut*nvertf*nDims);

  int threads = 32;
  int blocks = (ncells + threads - 1) / threads;

  int nbShare = sizeof(double)*nvertf*nDims;

  if (ijk2gmsh_quad.size() != nvertf)
  {
    auto ijk2gmsh_quad_h = tg_funcs::structured_to_gmsh_quad(nvert);
    ijk2gmsh_quad.assign(ijk2gmsh_quad_h.data(), nvertf);
  }

  switch(nvert)
  {
    case 8:
      fillCutMap<3,8><<<blocks, threads, nbShare>>>(*this, cutFaces, nCut, nvertf, cutMap, cutType);
      break;
//    case 27:
//      fillCutMap<3,27><<<blocks, threads, nbShare>>>(*this, cutFaces, nCut, nvertf, cutMap, cutType);
//      break;
//    case 64:
//      fillCutMap<3,64><<<blocks, threads, nbShare>>>(*this, cutFaces, nCut, nvertf, cutMap, cutType);
//      break;
//    case 125:
//      fillCutMap<3,125><<<blocks, threads, nbShare>>>(*this, cutFaces, nCut, nvertf, cutMap, cutType);
//      break;
    default:
      printf("nvert = %d\n",nvert);
      ThrowException("nvert case not implemented for directCut on device");
  }

  cudaDeviceSynchronize();
  check_error();

  cuda_copy_d2h(cutMap.flag.data(), cutFlag, ncells);
}
