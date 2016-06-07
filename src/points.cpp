/*!
 * \file points.cpp
 * \brief Functions related to quadrature point locations and weights
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * This file is part of the Tioga software library
 *
 * Tioga  is a tool for overset grid assembly on parallel distributed systems
 * Copyright (C) 2015-2016 Jay Sitaraman
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
#include <sstream>

#include "points.hpp"

std::vector<point> getLocCGPts(int order, int nDims)
{
  int nPts = (order+1) * (order+1);
  if (nDims == 3) nPts *= (order+1);

  std::vector<point> outPts(nPts);

  double dxi = 2. / order;

  if (nDims == 2)
  {
    for (int i = 0; i < order+1; i++) {
      for (int j = 0; j < order+1; j++) {
        outPts[j+i*(order+1)].x = -1 + dxi*j;
        outPts[j+i*(order+1)].y = -1 + dxi*i;
      }
    }
  }
  else
  {
    for (int k = 0; k < order+1; k++) {
      for (int j = 0; j < order+1; j++) {
        for (int i = 0; i < order+1; i++) {
          outPts[i+(order+1)*(j+(order+1)*k)].x = -1 + dxi*i;
          outPts[i+(order+1)*(j+(order+1)*k)].y = -1 + dxi*j;
          outPts[i+(order+1)*(j+(order+1)*k)].z = -1 + dxi*k;
        }
      }
    }
  }

  return outPts;
}

std::vector<point> getLocSpts(int eType, int order, std::string sptsType)
{
  std::vector<point> outPts;

  if (eType == QUAD) {
    // Tensor-product element
    std::vector<double> spts1D = getPts1D(sptsType,order);
    outPts.resize((order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        outPts[j+i*(order+1)].x = spts1D[j];
        outPts[j+i*(order+1)].y = spts1D[i];
      }
    }
  }
  else if (eType == HEX) {
    // Tensor-product element
    std::vector<double> spts1D = getPts1D(sptsType,order);
    outPts.resize((order+1)*(order+1)*(order+1));
    for (int k=0; k<order+1; k++) {
      for (int j=0; j<order+1; j++) {
        for (int i=0; i<order+1; i++) {
          outPts[i+(order+1)*(j+(order+1)*k)].x = spts1D[i];
          outPts[i+(order+1)*(j+(order+1)*k)].y = spts1D[j];
          outPts[i+(order+1)*(j+(order+1)*k)].z = spts1D[k];
        }
      }
    }
  }

  return outPts;
}

std::vector<point> getLocFpts(int eType, int order, std::string sptsType)
{
  std::vector<point> outPts;
  std::vector<double> pts1D;

  if (eType == QUAD) {
    outPts.resize(4*(order+1));
    pts1D = getPts1D(sptsType,order);
    for (int i=0; i<order+1; i++) {
      // Face 0
      outPts[i].x = pts1D[i];
      outPts[i].y = -1.;
      // Face 1
      outPts[i+order+1].x = 1.;
      outPts[i+order+1].y = pts1D[i];
      // Face 2
      outPts[i+2*(order+1)].x = pts1D[order-i];
      outPts[i+2*(order+1)].y = 1.;
      // Face 3
      outPts[i+3*(order+1)].x = -1.;
      outPts[i+3*(order+1)].y = pts1D[order-i];
    }
  }
  else if (eType == HEX) {
    int P12 = (order+1)*(order+1);
    outPts.resize(6*P12);
    pts1D = getPts1D(sptsType,order);
    // Flux points are ordered such that, as seen from inside the
    // element, the id's increase btm-left->top-right fashion on
    // each face, starting with lowest dimension first ('x' or 'y')
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        int ind = i+j*(order+1);
        // Face 0 - bottom
        outPts[ind].x = pts1D[i];
        outPts[ind].y = pts1D[j];
        outPts[ind].z = -1.;
        // Face 1 - top
        outPts[ind+P12].x = pts1D[order-i];
        outPts[ind+P12].y = pts1D[j];
        outPts[ind+P12].z = 1.;
        // Face 2 - left
        outPts[ind+2*P12].x = -1;
        outPts[ind+2*P12].y = pts1D[i];
        outPts[ind+2*P12].z = pts1D[j];
        // Face 3 - right
        outPts[ind+3*P12].x = 1;
        outPts[ind+3*P12].y = pts1D[order-i];
        outPts[ind+3*P12].z = pts1D[j];
        // Face 4 - front
        outPts[ind+4*P12].x = pts1D[order-i];
        outPts[ind+4*P12].y = -1;
        outPts[ind+4*P12].z = pts1D[j];
        // Face 5 - back
        outPts[ind+5*P12].x = pts1D[i];
        outPts[ind+5*P12].y = 1;
        outPts[ind+5*P12].z = pts1D[j];
      }
    }
  }

  return outPts;
}

std::vector<point> getLocPpts(int eType, int order, std::string sptsType)
{
  std::vector<point> outPts;
  std::vector<double> pts1D;

  int nPts1D = order + 3;
  pts1D = getPts1D(sptsType,order);
  pts1D.insert(pts1D.begin(),1,-1.);
  pts1D.push_back(1);

  if (eType == QUAD) {
    outPts.resize(nPts1D*nPts1D);
    for (int ppt=0; ppt<outPts.size(); ppt++) {
      int i = ppt%(order+3);
      int j = ppt/(order+3);
      outPts[ppt].x = pts1D[i];
      outPts[ppt].y = pts1D[j];
    }
  }
  else if (eType == HEX) {
    outPts.resize(nPts1D*nPts1D*nPts1D);
    for (int ppt = 0; ppt < outPts.size(); ppt++) {
      int k = ppt/(nPts1D*nPts1D);
      int j = (ppt-nPts1D*nPts1D*k)/nPts1D;
      int i = ppt - nPts1D*j - nPts1D*nPts1D*k;
      outPts[ppt].x = pts1D[i];
      outPts[ppt].y = pts1D[j];
      outPts[ppt].z = pts1D[k];
    }
  }

  return outPts;
}

std::vector<double> getPts1D(std::string ptsType, int order)
{
  std::vector<double> outPts(order+1);

  if (!ptsType.compare("Equidistant")) {
    double dxi = 2. / order;
    for (int i = 0; i < order+1; i++)
      outPts[i] = -1 + dxi*i;
  }
  else if (!ptsType.compare("Legendre")) { // Gauss-Legendre
    if (order == 0) {
      outPts[0] =  0.0;
    }
    else if(order == 1) {
      outPts[0] = -0.577350269189626;
      outPts[1] =  0.577350269189626;
    }
    else if(order == 2) {
      outPts[0] = -0.774596669241483;
      outPts[1] =  0.000000000000000;
      outPts[2] =  0.774596669241483;
    }
    else if(order == 3) {
      outPts[0] = -0.861136311594053;
      outPts[1] = -0.339981043584856;
      outPts[2] =  0.339981043584856;
      outPts[3] =  0.861136311594053;
    }
    else if(order == 4) {
      outPts[0] = -0.906179845938664;
      outPts[1] = -0.538469310105683;
      outPts[2] =  0.000000000000000;
      outPts[3] =  0.538469310105683;
      outPts[4] =  0.906179845938664;
    }
    else if(order == 5) {
      outPts[0] = -0.932469514203152;
      outPts[1] = -0.661209386466264;
      outPts[2] = -0.238619186083197;
      outPts[3] =  0.238619186083197;
      outPts[4] =  0.661209386466265;
      outPts[5] =  0.932469514203152;
    }
    else if(order == 6) {
      outPts[0] = -0.949107912342758;
      outPts[1] = -0.741531185599395;
      outPts[2] = -0.405845151377397;
      outPts[3] =  0.000000000000000;
      outPts[4] =  0.405845151377397;
      outPts[5] =  0.741531185599394;
      outPts[6] =  0.949107912342758;
    }
    else if(order == 7) {
      outPts[0] = -0.960289856497536;
      outPts[1] = -0.796666477413627;
      outPts[2] = -0.525532409916329;
      outPts[3] = -0.183434642495650;
      outPts[4] =  0.183434642495650;
      outPts[5] =  0.525532409916329;
      outPts[6] =  0.796666477413627;
      outPts[7] =  0.960289856497536;
    }
    else if(order == 8) {
      outPts[0] = -0.968160239507626;
      outPts[1] = -0.836031107326636;
      outPts[2] = -0.613371432700591;
      outPts[3] = -0.324253423403809;
      outPts[4] =  0.000000000000000;
      outPts[5] =  0.324253423403809;
      outPts[6] =  0.613371432700591;
      outPts[7] =  0.836031107326636;
      outPts[8] =  0.968160239507626;
    }
    else if(order == 9) {
      outPts[0] = -0.973906528517172;
      outPts[1] = -0.865063366688985;
      outPts[2] = -0.679409568299024;
      outPts[3] = -0.433395394129247;
      outPts[4] = -0.148874338981631;
      outPts[5] =  0.148874338981632;
      outPts[6] =  0.433395394129247;
      outPts[7] =  0.679409568299024;
      outPts[8] =  0.865063366688984;
      outPts[9] =  0.973906528517171;
    }
    else if(order == 10) {
      outPts[0] = -0.978228658146057;
      outPts[1] = -0.887062599768095;
      outPts[2] = -0.730152005574049;
      outPts[3] = -0.519096129206812;
      outPts[4] = -0.269543155952345;
      outPts[5] =  0.000000000000000;
      outPts[6] =  0.269543155952345;
      outPts[7] =  0.519096129206812;
      outPts[8] =  0.730152005574049;
      outPts[9] =  0.887062599768095;
      outPts[10]=  0.978228658146057;
    }
    else {
      std::stringstream ss; ss << order;
      std::string errMsg = "Gauss-Ledgendre point locations for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }
  else if (!ptsType.compare("Lobatto")) { // Gauss-Lobatto
    if (order == 0) {
        outPts[0] =  0.0;
    }
    else if(order == 1) {
      outPts[0] = -1.000000000000000;
      outPts[1] =  1.000000000000000;
    }
    else if(order == 2) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.000000000000000;
      outPts[2] =  1.000000000000000;
    }
    else if(order == 3) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.447213595499958;
      outPts[2] =  0.447213595499958;
      outPts[3] =  1.000000000000000;
    }
    else if(order == 4) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.654653670707977;
      outPts[2] = -0.000000000000000;
      outPts[3] =  0.654653670707977;
      outPts[4] =  1.000000000000000;
    }
    else if(order == 5) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.765055323929465;
      outPts[2] = -0.285231516480645;
      outPts[3] =  0.285231516480645;
      outPts[4] =  0.765055323929465;
      outPts[5] =  1.000000000000000;
    }
    else if(order == 6) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.830223896278567;
      outPts[2] = -0.468848793470714;
      outPts[3] =  0.000000000000000;
      outPts[4] =  0.468848793470714;
      outPts[5] =  0.830223896278567;
      outPts[6] =  1.000000000000000;
    }
    else if(order == 7) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.871740148509607;
      outPts[2] = -0.591700181433142;
      outPts[3] = -0.209299217902479;
      outPts[4] =  0.209299217902479;
      outPts[5] =  0.591700181433142;
      outPts[6] =  0.871740148509607;
      outPts[7] =  1.000000000000000;
    }
    else if(order == 8) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.899757995411460;
      outPts[2] = -0.677186279510738;
      outPts[3] = -0.363117463826178;
      outPts[4] = -0.000000000000000;
      outPts[5] =  0.363117463826178;
      outPts[6] =  0.677186279510738;
      outPts[7] =  0.899757995411460;
      outPts[8] =  1.000000000000000;
    }
    else if(order == 9) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.919533908166459;
      outPts[2] = -0.738773865105505;
      outPts[3] = -0.477924949810445;
      outPts[4] = -0.165278957666387;
      outPts[5] =  0.165278957666387;
      outPts[6] =  0.477924949810444;
      outPts[7] =  0.738773865105505;
      outPts[8] =  0.919533908166459;
      outPts[9] =  1.000000000000000;
    }
    else if(order == 10) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.934001430408059;
      outPts[2] = -0.784483473663144;
      outPts[3] = -0.565235326996205;
      outPts[4] = -0.295758135586939;
      outPts[5] =  0.000000000000000;
      outPts[6] =  0.295758135586939;
      outPts[7] =  0.565235326996205;
      outPts[8] =  0.784483473663144;
      outPts[9] =  0.934001430408059;
      outPts[10] =  1.000000000000000;
    }
    else {
      std::stringstream ss;  ss << order;
      std::string errMsg = "Gauss-Lobatto point locations for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }
  else {
    std::string str = "Points type not recognized: " + ptsType;
    FatalError(str.c_str());
  }

  return outPts;
}

std::vector<double> getQptWeights(int order, int nDims)
{
  // Tensor-product elements
  std::vector<double> qwts1D = getQptWeights1D(order);
  if (nDims == 1) return qwts1D;

  std::vector<double> outWts;
  if (nDims == 2) {
    outWts.resize((order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        outWts[j+i*(order+1)] = qwts1D[i]*qwts1D[j];
      }
    }
  }
  else if (nDims == 3) {
    outWts.resize((order+1)*(order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        for (int k=0; k<order+1; k++) {
          outWts[k+(order+1)*(j+i*(order+1))] = qwts1D[i]*qwts1D[j]*qwts1D[k];
        }
      }
    }
  }

  return outWts;
}


std::vector<double> getQptWeights1D(int order)
{
  // Order here refers to the order of a polynomial fit through
  // the Gauss points, not the order of accuracy of integration
  // using the same number of points

  std::vector<double> outWts(order+1);

  if (order == 0) {
    outWts[0] =  2.0;
  }
  else if(order == 1) {
    outWts[0] = 1.0;
    outWts[1] = 1.0;
  }
  else if(order == 2) {
    outWts[0] = 0.5555555555555556;
    outWts[1] = 0.8888888888888888;
    outWts[2] = 0.5555555555555556;
  }
  else if(order == 3) {
    outWts[0] = 0.3478548451374538;
    outWts[1] = 0.6521451548625461;
    outWts[2] = 0.6521451548625461;
    outWts[3] = 0.3478548451374538;
  }
  else if(order == 4) {
    outWts[0] = 0.2369268850561891;
    outWts[1] = 0.4786286704993665;
    outWts[2] = 0.5688888888888889;
    outWts[3] = 0.4786286704993665;
    outWts[4] = 0.2369268850561891;
  }
  else if (order == 5) {
    outWts[0] = 0.1713244923791704;
    outWts[1] = 0.3607615730481386;
    outWts[2] = 0.4679139345726910;
    outWts[3] = 0.4679139345726910;
    outWts[4] = 0.3607615730481386;
    outWts[5] = 0.1713244923791704;
  }
  else if (order == 6) {
    outWts[0] = 0.1294849661688697;
    outWts[1] = 0.2797053914892766;
    outWts[2] = 0.3818300505051189;
    outWts[3] = 0.4179591836734694;
    outWts[4] = 0.3818300505051189;
    outWts[5] = 0.2797053914892766;
    outWts[6] = 0.1294849661688697;
  }
  else if (order == 7) {
    outWts[0] = 0.1012285362903763;
    outWts[1] = 0.2223810344533745;
    outWts[2] = 0.3137066458778873;
    outWts[3] = 0.3626837833783620;
    outWts[4] = 0.3626837833783620;
    outWts[5] = 0.3137066458778873;
    outWts[6] = 0.2223810344533745;
    outWts[7] = 0.1012285362903763;
  }
  else if (order == 8) {
    outWts[0] = 0.0812743883615744;
    outWts[1] = 0.1806481606948574;
    outWts[2] = 0.2606106964029354;
    outWts[3] = 0.3123470770400029;
    outWts[4] = 0.3302393550012598;
    outWts[5] = 0.3123470770400029;
    outWts[6] = 0.2606106964029354;
    outWts[7] = 0.1806481606948574;
    outWts[8] = 0.0812743883615744;
  }
  else if (order == 9) {
    outWts[0] = 0.0666713443086881;
    outWts[1] = 0.1494513491505806;
    outWts[2] = 0.2190863625159820;
    outWts[3] = 0.2692667193099963;
    outWts[4] = 0.2955242247147529;
    outWts[5] = 0.2955242247147529;
    outWts[6] = 0.2692667193099963;
    outWts[7] = 0.2190863625159820;
    outWts[8] = 0.1494513491505806;
    outWts[9] = 0.0666713443086881;
  }
  else if (order == 10) {
    outWts[0 ] = 0.0556685671161737;
    outWts[1 ] = 0.1255803694649046;
    outWts[2 ] = 0.1862902109277343;
    outWts[3 ] = 0.2331937645919905;
    outWts[4 ] = 0.2628045445102467;
    outWts[5 ] = 0.2729250867779006;
    outWts[6 ] = 0.2628045445102467;
    outWts[7 ] = 0.2331937645919905;
    outWts[8 ] = 0.1862902109277343;
    outWts[9 ] = 0.1255803694649046;
    outWts[10] = 0.0556685671161737;
  }
  else {
    std::stringstream ss; ss << order;
    std::string errMsg = "Gauss quadrature weights for order " + ss.str() + " not implemented.";
    FatalError(errMsg.c_str());
  }

  return outWts;
}

void getQuadRuleTet(int order, std::vector<point> &locQpts, std::vector<double> &weights)
{
  /* Linbo Zhang, Tau Cui and Hui Liu, "A Set of Symmetric Quadrature Rules on
   * Triangles and Tetrahedra." J. Comp. Math., 2009.
   * See also lsec.cc.ac.cn/phg to obtain source files containing explicit
   * quadrature rules. */
  switch (order) {
    case 1:
      locQpts.resize(1);
      locQpts[0].x = .25; locQpts[0].y = .25; locQpts[0].z = .25;
      weights = {1};
      break;

    case 2: {
      locQpts.resize(4);
      double loc1 = .138196601125011;
      double loc2 = .585410196624969;
      locQpts[0].x = loc1; locQpts[0].y = loc1; locQpts[0].z = loc1;
      locQpts[1].x = loc2; locQpts[1].y = loc1; locQpts[1].z = loc1;
      locQpts[2].x = loc1; locQpts[2].y = loc2; locQpts[2].z = loc1;
      locQpts[3].x = loc1; locQpts[3].y = loc1; locQpts[3].z = loc2;
      weights = {.25,.25,.25,.25};
      break;
    }
    case 3: {
      locQpts.resize(8);
      double loc1 = .328054696711427;
      double loc2 = 1-3*loc1;
      locQpts[0].x = loc1; locQpts[0].y = loc1; locQpts[0].z = loc1;
      locQpts[1].x = loc2; locQpts[1].y = loc1; locQpts[1].z = loc1;
      locQpts[2].x = loc1; locQpts[2].y = loc2; locQpts[2].z = loc1;
      locQpts[3].x = loc1; locQpts[3].y = loc1; locQpts[3].z = loc2;
      loc1 = .328054696711427;
      loc2 = 1-3*loc1;
      locQpts[4].x = loc1; locQpts[4].y = loc1; locQpts[4].z = loc1;
      locQpts[5].x = loc2; locQpts[5].y = loc1; locQpts[5].z = loc1;
      locQpts[6].x = loc1; locQpts[6].y = loc2; locQpts[6].z = loc1;
      locQpts[7].x = loc1; locQpts[7].y = loc1; locQpts[7].z = loc2;
      double wt1 = .138527966511862;
      double wt2 = .111472033488138;
      weights = {wt1,wt1,wt1,wt1,wt2,wt2,wt2,wt2};
      break;
    }
    case 4: {
      locQpts.resize(14);
      double loc1 = .0927352503108912;
      double loc2 = 1-3*loc1;
      locQpts[0].x = loc1; locQpts[0].y = loc1; locQpts[0].z = loc1;
      locQpts[1].x = loc2; locQpts[1].y = loc1; locQpts[1].z = loc1;
      locQpts[2].x = loc1; locQpts[2].y = loc2; locQpts[2].z = loc1;
      locQpts[3].x = loc1; locQpts[3].y = loc1; locQpts[3].z = loc2;
      loc1 = .3108859192633006;
      loc2 = 1-3*loc1;
      locQpts[4].x = loc1; locQpts[4].y = loc1; locQpts[4].z = loc1;
      locQpts[5].x = loc2; locQpts[5].y = loc1; locQpts[5].z = loc1;
      locQpts[6].x = loc1; locQpts[6].y = loc2; locQpts[6].z = loc1;
      locQpts[7].x = loc1; locQpts[7].y = loc1; locQpts[7].z = loc2;
      loc1 = .0455037041256497;
      loc2 = 0.5-loc1;
      locQpts[8].x  = loc1; locQpts[8].y  = loc1; locQpts[8].z = loc2;
      locQpts[9].x  = loc1; locQpts[9].y  = loc2; locQpts[9].z = loc1;
      locQpts[10].x = loc1; locQpts[10].y = loc2; locQpts[10].z = loc2;
      locQpts[11].x = loc2; locQpts[11].y = loc1; locQpts[11].z = loc2;
      locQpts[12].x = loc2; locQpts[12].y = loc1; locQpts[12].z = loc1;
      locQpts[13].x = loc2; locQpts[13].y = loc2; locQpts[13].z = loc1;
      double wt1 = .073493043116362;
      double wt2 = .112687925718016;
      double wt3 = .042546020777082;
      weights = {wt1,wt1,wt1,wt1,wt2,wt2,wt2,wt2,wt3,wt3,wt3,wt3,wt3,wt3};
      break;
    }
    case 6: {
      locQpts.resize(24);
      double loc1, loc2, loc3;
      int N = 0;
      loc1 = .21460287125915202928883921938628499;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      loc1 = .04067395853461135311557944895641006;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      loc1 = .32233789014227551034399447076249213;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      loc1 = .06366100187501752529923552760572698;
      loc2 = .60300566479164914136743113906093969;
      loc3 = 1. - 2.*loc1 - loc2;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      double wt1 = .03992275025816749209969062755747998;
      double wt2 = .01007721105532064294801323744593686;
      double wt3 = .05535718154365472209515327785372602;
      double wt4 = 27./560.;
      weights = {wt1,wt1,wt1,wt1,
                 wt2,wt2,wt2,wt2,
                 wt3,wt3,wt3,wt3,
                 wt4,wt4,wt4,wt4,wt4,wt4,wt4,wt4,wt4,wt4,wt4,wt4};
      break;
    }
    case 8: {
      locQpts.resize(46);
      double loc1, loc2, loc3;
      int N = 0;
      loc1 = .03967542307038990126507132953938949;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      N += 4;
      loc1 = .31448780069809631378416056269714830;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      N += 4;
      loc1 = .1019866930627033;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      N += 4;
      loc1 = .18420369694919151227594641734890918;
      loc2 = 1-3*loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc1;
      locQpts[N+1].x = loc2; locQpts[N+1].y = loc1; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc1;
      locQpts[N+3].x = loc1; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      N += 4;
      loc1 = .06343628775453989240514123870189827;
      loc2 = 0.5-loc1;
      locQpts[N+0].x = loc1; locQpts[N+0].y = loc1; locQpts[N+0].z = loc2;
      locQpts[N+1].x = loc1; locQpts[N+1].y = loc2; locQpts[N+1].z = loc1;
      locQpts[N+2].x = loc1; locQpts[N+2].y = loc2; locQpts[N+2].z = loc2;
      locQpts[N+3].x = loc2; locQpts[N+3].y = loc1; locQpts[N+3].z = loc2;
      locQpts[N+4].x = loc2; locQpts[N+4].y = loc1; locQpts[N+4].z = loc1;
      locQpts[N+5].x = loc2; locQpts[N+5].y = loc2; locQpts[N+5].z = loc1;
      N += 6;
      loc1 = .02169016206772800480266248262493018;
      loc2 = .71993192203946593588943495335273478;
      loc3 = 1. - 2.*loc1 - loc2;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      N+=12;
      loc1 = .20448008063679571424133557487274534;
      loc2 = .58057719012880922417539817139062041;
      loc3 = 1. - 2.*loc1 - loc2;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      double wt1 = .00639714777990232132145142033517302;
      double wt2 = .04019044802096617248816115847981783;
      double wt3 = .02430797550477032117486910877192260;
      double wt4 = .05485889241369744046692412399039144;
      double wt5 = .03571961223409918246495096899661762;
      double wt6 = .00718319069785253940945110521980376;
      double wt7 = .01637218194531911754093813975611913;
      weights = {wt1,wt1,wt1,wt1,
                 wt2,wt2,wt2,wt2,
                 wt3,wt3,wt3,wt3,
                 wt4,wt4,wt4,wt4,
                 wt5,wt5,wt5,wt5,wt5,wt5,
                 wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,
                 wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7};
      break;
    }
    case 10: {
      locQpts.resize(79);
      double loc1 = .25;
      locQpts[0].x = loc1; locQpts[0].y = loc1; locQpts[0].z = loc1;
      loc1 = .11425191803006935688146412277598412;
      double loc2 = 1-3*loc1;
      locQpts[1].x = loc1; locQpts[1].y = loc1; locQpts[1].z = loc1;
      locQpts[2].x = loc2; locQpts[2].y = loc1; locQpts[2].z = loc1;
      locQpts[3].x = loc1; locQpts[3].y = loc2; locQpts[3].z = loc1;
      locQpts[4].x = loc1; locQpts[4].y = loc1; locQpts[4].z = loc2;
      loc1 = .01063790234539248531264164411274776;
      loc2 = 1-3*loc1;
      locQpts[5].x = loc1; locQpts[5].y = loc1; locQpts[5].z = loc1;
      locQpts[6].x = loc2; locQpts[6].y = loc1; locQpts[6].z = loc1;
      locQpts[7].x = loc1; locQpts[7].y = loc2; locQpts[7].z = loc1;
      locQpts[8].x = loc1; locQpts[8].y = loc1; locQpts[8].z = loc2;
      loc1 = .31274070833535645859816704980806110;
      loc2 = 1-3*loc1;
      locQpts[9].x = loc1;  locQpts[9].y = loc1;  locQpts[9].z = loc1;
      locQpts[10].x = loc2; locQpts[10].y = loc1; locQpts[10].z = loc1;
      locQpts[11].x = loc1; locQpts[11].y = loc2; locQpts[11].z = loc1;
      locQpts[12].x = loc1; locQpts[12].y = loc1; locQpts[12].z = loc2;
      loc1 = .01631296303281644;
      loc2 = 0.5-loc1;
      locQpts[13].x = loc1; locQpts[13].y = loc1; locQpts[13].z = loc2;
      locQpts[14].x = loc1; locQpts[14].y = loc2; locQpts[14].z = loc1;
      locQpts[15].x = loc1; locQpts[15].y = loc2; locQpts[15].z = loc2;
      locQpts[16].x = loc2; locQpts[16].y = loc1; locQpts[16].z = loc2;
      locQpts[17].x = loc2; locQpts[17].y = loc1; locQpts[17].z = loc1;
      locQpts[18].x = loc2; locQpts[18].y = loc2; locQpts[18].z = loc1;
      int N = 19;
      loc1 = .03430622963180452385835196582344460;
      loc2 = .59830121060139461905983787517050400;
      double loc3 = 1. - 2.*loc1 - loc2;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      N += 12;
      loc1 = .12346418534551115945916818783743644;
      loc2 = .47120066204746310257913700590727081;
      loc3 = 1. - 2.*loc1 - loc2;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      loc1 = .40991962933181117418479812480531207;
      loc2 = .16546413290740130923509687990363569;
      loc3 = 1. - 2.*loc1 - loc2;
      N += 12;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      loc1 = .17397243903011716743177479785668929;
      loc2 = .62916375300275643773181882027844514;
      loc3 = 1. - 2.*loc1 - loc2;
      N += 12;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      loc1 = .03002157005631784150255786784038011;
      loc2 = .81213056814351208262160080755918730;
      loc3 = 1. - 2.*loc1 - loc2;
      N += 12;
      locQpts[N+0].x  = loc1; locQpts[N+0].y  = loc1; locQpts[N+0].z  = loc2;
      locQpts[N+1].x  = loc1; locQpts[N+1].y  = loc1; locQpts[N+1].z  = loc3;
      locQpts[N+2].x  = loc1; locQpts[N+2].y  = loc2; locQpts[N+2].z  = loc1;
      locQpts[N+3].x  = loc1; locQpts[N+3].y  = loc2; locQpts[N+3].z  = loc3;
      locQpts[N+4].x  = loc1; locQpts[N+4].y  = loc3; locQpts[N+4].z  = loc1;
      locQpts[N+5].x  = loc1; locQpts[N+5].y  = loc3; locQpts[N+5].z  = loc2;
      locQpts[N+6].x  = loc2; locQpts[N+6].y  = loc1; locQpts[N+6].z  = loc1;
      locQpts[N+7].x  = loc2; locQpts[N+7].y  = loc1; locQpts[N+7].z  = loc3;
      locQpts[N+8].x  = loc2; locQpts[N+8].y  = loc3; locQpts[N+8].z  = loc1;
      locQpts[N+9].x  = loc3; locQpts[N+9].y  = loc1; locQpts[N+9].z  = loc1;
      locQpts[N+10].x = loc3; locQpts[N+10].y = loc1; locQpts[N+10].z = loc2;
      locQpts[N+11].x = loc3; locQpts[N+11].y = loc2; locQpts[N+11].z = loc1;
      double wt1 = .04574189830483037077884770618329337;
      double wt2 = .01092727610912416907498417206565671;
      double wt3 = .00055352334192264689534558564012282;
      double wt4 = .02569337913913269580782688316792080;
      double wt5 = .00055387649657283109312967562590035;
      double wt6 = .01044842402938294329072628200105773;
      double wt7 = .02513844602651287118280517785487423;
      double wt8 = .01178620679249594711782155323755017;
      double wt9 = .01332022473886650471019828463616468;
      double wt10 = .00615987577565961666092767531756180;
      weights = {wt1,
                 wt2,wt2,wt2,wt2,
                 wt3,wt3,wt3,wt3,
                 wt4,wt4,wt4,wt4,
                 wt5,wt5,wt5,wt5,wt5,wt5,
                 wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,wt6,
                 wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,wt7,
                 wt8,wt8,wt8,wt8,wt8,wt8,wt8,wt8,wt8,wt8,wt8,wt8,
                 wt9,wt9,wt9,wt9,wt9,wt9,wt9,wt9,wt9,wt9,wt9,wt9,
                 wt10,wt10,wt10,wt10,wt10,wt10,wt10,wt10,wt10,wt10,wt10,wt10};
      break;
    }
    default:  {
      std::stringstream ss; ss << order;
      std::string errMsg = "Tetrahedron quadrature rules for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }
}

void getQuadRuleTri(int order, std::vector<point> &locQpts, std::vector<double> &weights)
{
  /* Linbo Zhang, Tau Cui and Hui Liu, "A Set of Symmetric Quadrature Rules on
   * Triangles and Tetrahedra." J. Comp. Math., 2009.
   * See also lsec.cc.ac.cn/phg to obtain source files containing explicit
   * quadrature rules. */
  switch(order) {
    case 1: {
      locQpts = {point({1./3.,1./3.,0.})};
      weights = {1.};
      break;
    }
    case 2: {
      locQpts.resize(3);
      double loc1 = 1./6;
      double loc2 = 1.-2.*loc1;
      locQpts[0] = point({loc2,loc1,0});
      locQpts[1] = point({loc1,loc2,0});
      locQpts[2] = point({loc1,loc1,0});
      weights = {1./3.,1./3.,1./3.};
      break;
    }
    case 3: {
      locQpts.resize(6);
      double loc1 = .1628828503958919;
      double loc2 = 1-2.*loc1;
      locQpts[0] = point({loc2, loc1, 0.});
      locQpts[1] = point({loc1, loc2, 0.});
      locQpts[2] = point({loc1, loc1, 0.});
      loc1 = .4779198835675637;
      loc2 = 1-2.*loc1;
      locQpts[3] = point({loc2, loc1, 0.});
      locQpts[4] = point({loc1, loc2, 0.});
      locQpts[5] = point({loc1, loc1, 0.});
      double wt1 = .2811498024409796;
      double wt2 = .0521835308923537;
      weights = {wt1,wt1,wt1,wt2,wt2,wt2};
      break;
    }
    case 4: {
      locQpts.resize(6);
      double loc1 = .4459484909159649; // = (8-sqrt(10)+sqrt(38-44*sqrt(2/5))) / 18
      double loc2 = 1-2.*loc1;
      locQpts[0] = point({loc2, loc1, 0.});
      locQpts[1] = point({loc1, loc2, 0.});
      locQpts[2] = point({loc1, loc1, 0.});
      loc1 = .0915762135097707;        // = (8-sqrt(10)-sqrt(38-44*sqrt(2/5))) / 18
      loc2 = 1-2.*loc1;
      locQpts[3] = point({loc2, loc1, 0.});
      locQpts[4] = point({loc1, loc2, 0.});
      locQpts[5] = point({loc1, loc1, 0.});
      double wt1 = .2233815896780115; // = (620 + sqrt(213125 - 53320 * sqrt(10))) / 3720
      double wt2 = .1099517436553219; // = (620 - sqrt(213125 - 53320 * sqrt(10))) / 3720
      weights = {wt1,wt1,wt1,wt2,wt2,wt2};
      break;
    }
    case 5: {
      locQpts.resize(7);
      double loc1 = .1012865073234563;
      double loc2 = 1-2.*loc1;
      locQpts[0] = point({loc2, loc1, 0.});
      locQpts[1] = point({loc1, loc2, 0.});
      locQpts[2] = point({loc1, loc1, 0.});
      loc1 = .4701420641051151;
      loc2 = 1-2.*loc1;
      locQpts[3] = point({loc2, loc1, 0.});
      locQpts[4] = point({loc1, loc2, 0.});
      locQpts[5] = point({loc1, loc1, 0.});
      loc1 = 1./3.;
      locQpts[6] = point({loc1, loc1, 0.});
      double wt1 = .1259391805448272; // = (155 - sqrt(15)) / 1200
      double wt2 = .1323941527885062; // = (155 + sqrt(15)) / 1200
      double wt3 = 9./40.;
      weights = {wt1,wt1,wt1,wt2,wt2,wt2,wt3};
      break;
    }
    case 6: {
      locQpts.resize(12);
      double loc1 = .0630890144915022;
      double loc2 = 1-2.*loc1;
      locQpts[0] = point({loc2, loc1, 0.});
      locQpts[1] = point({loc1, loc2, 0.});
      locQpts[2] = point({loc1, loc1, 0.});
      loc1 = .2492867451709104;
      loc2 = 1-2.*loc1;
      locQpts[3] = point({loc2, loc1, 0.});
      locQpts[4] = point({loc1, loc2, 0.});
      locQpts[5] = point({loc1, loc1, 0.});
      loc1 = .0531450498448169;
      loc2 = .3103524510337844;
      double loc3 = 1-loc1-loc2;
      locQpts[6]  = point({loc1, loc2, 0.});
      locQpts[7]  = point({loc1, loc3, 0.});
      locQpts[8]  = point({loc2, loc1, 0.});
      locQpts[9]  = point({loc2, loc3, 0.});
      locQpts[10] = point({loc3, loc1, 0.});
      locQpts[11] = point({loc3, loc2, 0.});
      double wt1 = .0508449063702068;
      double wt2 = .1167862757263794;
      double wt3 = .0828510756183736;
      weights = {wt1,wt1,wt1,wt2,wt2,wt2,wt3,wt3,wt3,wt3,wt3,wt3};
      break;
    }
    case 7: {
      locQpts.resize(15);
      double loc1 = .0282639241560763;
      double loc2 = 1-2.*loc1;
      locQpts[0] = point({loc2, loc1, 0.});
      locQpts[1] = point({loc1, loc2, 0.});
      locQpts[2] = point({loc1, loc1, 0.});
      loc1 = .4743113232672226;
      loc2 = 1-2.*loc1;
      locQpts[3] = point({loc2, loc1, 0.});
      locQpts[4] = point({loc1, loc2, 0.});
      locQpts[5] = point({loc1, loc1, 0.});
      loc1 = .2411433258498488;
      loc2 = 1-2.*loc1;
      locQpts[6] = point({loc2, loc1, 0.});
      locQpts[7] = point({loc1, loc2, 0.});
      locQpts[8] = point({loc1, loc1, 0.});
      loc1 = .7612227480245238;
      loc2 = .0462708777988089;
      double loc3 = 1-loc1-loc2;
      locQpts[9]  = point({loc1, loc2, 0.});
      locQpts[10] = point({loc1, loc3, 0.});
      locQpts[11] = point({loc2, loc1, 0.});
      locQpts[12] = point({loc2, loc3, 0.});
      locQpts[13] = point({loc3, loc1, 0.});
      locQpts[14] = point({loc3, loc2, 0.});
      double wt1 = .0135338625156656;
      double wt2 = .0789512544320110;
      double wt3 = .1286079278189061;
      double wt4 = .0561201442833754;
      weights = {wt1,wt1,wt1,wt2,wt2,wt2,wt3,wt3,wt3,wt4,wt4,wt4,wt4,wt4,wt4};
      break;
    }
    case 8: {
      locQpts.resize(16);
      locQpts[0] = point({1./3., 1./3., 0.});
      double loc1 = .1705693077517602;
      double loc2 = 1-2.*loc1;
      locQpts[1] = point({loc2, loc1, 0.});
      locQpts[2] = point({loc1, loc2, 0.});
      locQpts[3] = point({loc1, loc1, 0.});
      loc1 = .0505472283170310;
      loc2 = 1-2.*loc1;
      locQpts[4] = point({loc2, loc1, 0.});
      locQpts[5] = point({loc1, loc2, 0.});
      locQpts[6] = point({loc1, loc1, 0.});
      loc1 = .4592925882927232;
      loc2 = 1-2.*loc1;
      locQpts[7] = point({loc2, loc1, 0.});
      locQpts[8] = point({loc1, loc2, 0.});
      locQpts[9] = point({loc1, loc1, 0.});
      loc1 = .2631128296346381;
      loc2 = .00839477740995761;
      double loc3 = 1-loc1-loc2;
      locQpts[10] = point({loc1, loc2, 0.});
      locQpts[11] = point({loc1, loc3, 0.});
      locQpts[12] = point({loc2, loc1, 0.});
      locQpts[13] = point({loc2, loc3, 0.});
      locQpts[14] = point({loc3, loc1, 0.});
      locQpts[15] = point({loc3, loc2, 0.});
      double wt0 = .1443156076777872;
      double wt1 = .1032173705347183;
      double wt2 = .0324584976231981;
      double wt3 = .0950916342672846;
      double wt4 = .0272303141744350;
      weights = {wt0,wt1,wt1,wt1,wt2,wt2,wt2,wt3,wt3,wt3,wt4,wt4,wt4,wt4,wt4,wt4};
      break;
    }
    case 10: {
      locQpts.resize(25);
      locQpts[0] = point({1./3., 1./3., 0.});
      double loc1 = .4272731788467755;
      double loc2 = 1-2.*loc1;
      locQpts[1] = point({loc2, loc1, 0.});
      locQpts[2] = point({loc1, loc2, 0.});
      locQpts[3] = point({loc1, loc1, 0.});
      loc1 = .1830992224486750;
      loc2 = 1-2.*loc1;
      locQpts[4] = point({loc2, loc1, 0.});
      locQpts[5] = point({loc1, loc2, 0.});
      locQpts[6] = point({loc1, loc1, 0.});
      loc1 = .4904340197011305;
      loc2 = 1-2.*loc1;
      locQpts[7] = point({loc2, loc1, 0.});
      locQpts[8] = point({loc1, loc2, 0.});
      locQpts[9] = point({loc1, loc1, 0.});
      loc1 = .0125724455515805;
      loc2 = 1-2.*loc1;
      locQpts[10] = point({loc2, loc1, 0.});
      locQpts[11] = point({loc1, loc2, 0.});
      locQpts[12] = point({loc1, loc1, 0.});
      loc1 = .65426866792006614;
      loc2 = .30804600168524770;
      double loc3 = 1-loc1-loc2;
      locQpts[13] = point({loc1, loc2, 0.});
      locQpts[14] = point({loc1, loc3, 0.});
      locQpts[15] = point({loc2, loc1, 0.});
      locQpts[16] = point({loc2, loc3, 0.});
      locQpts[17] = point({loc3, loc1, 0.});
      locQpts[18] = point({loc3, loc2, 0.});
      loc1 = .12280457706855927;
      loc2 = 3.3371833739304786e-2;
      loc3 = 1-loc1-loc2;
      locQpts[19] = point({loc1, loc2, 0.});
      locQpts[20] = point({loc1, loc3, 0.});
      locQpts[21] = point({loc2, loc1, 0.});
      locQpts[22] = point({loc2, loc3, 0.});
      locQpts[23] = point({loc3, loc1, 0.});
      locQpts[24] = point({loc3, loc2, 0.});
      double wt0 = 8.093742879762288e-2;
      double wt1 = 7.729858800296312e-2;
      double wt2 = 7.845763861237173e-2;
      double wt3 = 1.746916799592949e-2;
      double wt4 = 4.292374184832828e-3;
      double wt5 = 3.746885821046764e-2;
      double wt6 = 2.694935259187996e-2;
      weights = {wt0,wt1,wt1,wt1,wt2,wt2,wt2,wt3,wt3,wt3,wt4,wt4,wt4,
                     wt5,wt5,wt5,wt5,wt5,wt5,wt6,wt6,wt6,wt6,wt6,wt6};
      break;
    }
    default:  {
      std::stringstream ss; ss << order;
      std::string errMsg = "Triangle quadrature rules for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }
}
