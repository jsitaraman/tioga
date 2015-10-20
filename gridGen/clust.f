c
c Tioga is a library for overset grid assembly on parallel distributed systems
c Copyright (C) 2015 Jay Sitaraman
c
c=======================================================================
      subroutine distrib(ncluster,cntrlpt,arclen,ds,out)
c
c  distribute points into out based on input using vinokur stretching
c  input - ncluster is number of points at which clustering is specified
c          cntrlpt are the indices
c          arclen are the distances  
c          ds are the spacings
c
c=======================================================================
      implicit real*8(a-h,o-z)
	parameter(ndim=1800)
c
      dimension cntrlpt(*),arclen(*),ds(*),out(*)
      dimension arcdist(ndim)
c
      do 10 n=1,ncluster-1
         node0 = cntrlpt(n)
         node1 = cntrlpt(n+1)
         alength = arclen(n+1)-arclen(n)
         ds0 = ds(n)/alength
         ds1 = ds(n+1)/alength
         call clust(arcdist,ds0,ds1,node0,node1)
         do 20 i=node0,node1
            out(i) = arclen(n) + alength*arcdist(i)
 20      continue
 10   continue
c
      return
      end

c=======================================================================
      subroutine clust(y,ds0,ds1,jmin,jmax)
c
c..one step corrector for clustering
c
c=======================================================================
      implicit real*8(a-h,o-z)
      dimension y(*)
c
      del = 1./(jmax-jmin)
      s0=del/ds0
      s1=del/ds1
      call clust2(y,s0,s1,jmin,jmax)
c
      alpha0 = (y(jmin+1)-y(jmin))/ds0
      alpha1 = (y(jmax)-y(jmax-1))/ds1
      s0=alpha0*del/ds0
      s1=alpha1*del/ds1
      call clust2(y,s0,s1,jmin,jmax)
c
      return
      end

c=======================================================================
      subroutine clust2(y,s0,s1,jmin,jmax)
c
c  vinokur's end point clustering function (should be run double 
c  precision if not on cray).  y(j) is stretched function between 0 & 1
c  (i.e. it is normalized and must be rescaled).
c
c=======================================================================
cdouble      implicit real*8(a-h,o-z)
      implicit real*8(a-h,o-z)
      dimension y(*)
      data eps/.001/
c
      jm1=jmax-1
      jp1=jmin+1
      del=1./(jmax-jmin)
      b=sqrt(s0*s1)
      a=b/s1
      y(jmax)=1.0
      y(jmin)=0.0
c
      if(b.gt.1.+eps) then
        dz=asinhf(b)
        coshdz=cosh(dz)
        omacdz=1.-a*coshdz
        asindz=a*sqrt(coshdz*coshdz-1.)
        do 11 j=jp1,jm1
          u=tanh(dz*(j-jmin)*del)
          y(j)=u/(asindz+omacdz*u)
 11     continue
      else
        if(b.gt.1.-eps) then
          twobmo=2.*(b-1.)
          onema=1.-a
          do 21 j=jp1,jm1
            x=(j-jmin)*del
            u=x*(1.+twobmo*(x-.5)*(1.-x))
            y(j)=u/(a+onema*u)
 21       continue
        else
          dz=asinf(b)
          cosdz=cos(dz)
          cscdz=1./sqrt(1.-cosdz*cosdz)
          w0=cscdz*(cosdz-1./a)
          z0=atan(w0)
          dw=cscdz*(a-cosdz)-w0
          odw=1./dw
          do  1 j=jp1,jm1
            y(j)=odw*(tan(dz*(j-jmin)*del+z0)-w0)
 1        continue
        endif
      endif
c
      return
      end


c=======================================================================
      function asinhf(u)
c=======================================================================
cdouble      implicit real*8(a-h,o-z)
      implicit real*8(a-h,o-z)
      data a1,a2,a3,a4,a5/-.15,.0573214285714,-.024907294878,.0077424460
     1899,-.0010794122691/
      data b0,b1,b2,b3,b4/-.0204176930892,.2490272170591,1.9496443322775
     1,-2.629454725241,8.5679591096315/
      data u1,u2/2.7829681178603,35.0539798452776/
c
      if(u.le.u1) then
        ub=u-1.
        asinhf=sqrt(6.*ub)*(((((a5*ub+a4)*ub+a3)*ub+a2)*ub+a1)*ub+1.)
      else
cdouble        v=dlog(u)
        v=log(u)
        w=1./u-1./u2
cdouble        asinhf=v+dlog(2.*v)*(1.+1./v)+(((b4*w+b3)*w+b2)*w+b1)*w+b0
        asinhf=v+log(2.*v)*(1.+1./v)+(((b4*w+b3)*w+b2)*w+b1)*w+b0
      endif
c
      return
      end

c=======================================================================
      function asinf(u)
c=======================================================================
cdouble      implicit real*8(a-h,o-z)
      implicit real*8(a-h,o-z)
      data a1,a2,a3,a4,a5/.15,.0573214285714,.0489742834696,
     1-.053337753213,.0758451335824/
      data b3,b4,b5,b6/-2.6449340668482,6.7947319658321,
     1-13.2055008110734,11.7260952338351/
      data u1,pi/.2693897165164,3.14159265358981/
c
      if(u.le.u1) then
        asinf=pi*((((((b6*u+b5)*u+b4)*u+b3)*u+1.)*u-1.)*u+1.)
      else
        ub=1.-u
cdouble        asinf=dsqrt(6.*ub)*(((((a5*ub+a4)*ub+a3)*ub+a2)*ub+a1)*ub+1.)
        asinf=sqrt(6.*ub)*(((((a5*ub+a4)*ub+a3)*ub+a2)*ub+a1)*ub+1.)
      endif
c
      return
      end



