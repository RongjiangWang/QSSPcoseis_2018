      subroutine gccmp(expl,clvd,ss12,ss11,ds31,ds23,csa,ssa,cs2a,ss2a,
     &                 irg,izrg1,izrg2,nzrg,y)
      use gcalloc
      implicit none
c
      integer*4 irg,izrg1,izrg2,nzrg
      real*8 expl,clvd,ss12,ss11,ds31,ds23,csa,ssa,cs2a,ss2a
      real*8 y(22,nzrg)
c
      integer*4 izrg
      real*8 dzz,dzs,dze,dsz,dss,dse,dez,des,dee
c
c     azimuth (phi) dependence:
c     explosion & clvd: 1.0
c     strike-slip m(1,2)=m(2,1):
c                 ur,ut,err,ert,etr,ett,epp,gr = sin(2*phi), up,erp,etp,epr,ept = cos(2*phi)
c     dip-slip m(1,3)=m(3,1)!=0:
c                 ur,ut,err,ert,etr,ett,epp,gr = cos(phi), up,erp,etp,epr,ept = sin(phi)
c     strike-slip m(1,1)=-m(2,2) => strike-slip m(1,2)=m(2,1) through phi -> phi+45°:
c                 ur,ut,err,ert,etr,ett,epp,gr = cos(2*phi), up,erp,etp,epr,ept = -sin(2*phi)
c     dip-slip m(2,3)=m(3,2)!=0 => dip-slip m(1,3)=m(3,1) through phi -> phi-90°:
c                 ur,ut,err,ert,etr,ett,epp,gr = sin(phi), up,erp,etp,epr,ept = -cos(phi)
c
      do izrg=izrg1,izrg2
c
c       1. Uz, 2. Ut, 3. Up
c
        y(1,izrg)=expl*ur(izrg,irg,1)+clvd*ur(izrg,irg,4)
     &          +(ss12*ss2a+ss11*cs2a)*ur(izrg,irg,2)
     &          +(ds31*csa+ds23*ssa)*ur(izrg,irg,3)
        y(2,izrg)=expl*ut(izrg,irg,1)+clvd*ut(izrg,irg,4)
     &          +(ss12*ss2a+ss11*cs2a)*ut(izrg,irg,2)
     &          +(ds31*csa+ds23*ssa)*ut(izrg,irg,3)
        y(3,izrg)=(ss12*cs2a-ss11*ss2a)*up(izrg,irg,2)
     &          +(ds31*ssa-ds23*csa)*up(izrg,irg,3)
c
        dzz=expl*err(izrg,irg,1)+clvd*err(izrg,irg,4)
     &     +(ss12*ss2a+ss11*cs2a)*err(izrg,irg,2)
     &     +(ds31*csa+ds23*ssa)*err(izrg,irg,3)
        dsz=expl*etr(izrg,irg,1)+clvd*etr(izrg,irg,4)
     &     +(ss12*ss2a+ss11*cs2a)*etr(izrg,irg,2)
     &     +(ds31*csa+ds23*ssa)*etr(izrg,irg,3)
        dez=(ss12*cs2a-ss11*ss2a)*epr(izrg,irg,2)
     &     +(ds31*ssa-ds23*csa)*epr(izrg,irg,3)
c
        dzs=expl*ert(izrg,irg,1)+clvd*ert(izrg,irg,4)
     &     +(ss12*ss2a+ss11*cs2a)*ert(izrg,irg,2)
     &     +(ds31*csa+ds23*ssa)*ert(izrg,irg,3)
        dss=expl*ett(izrg,irg,1)+clvd*ett(izrg,irg,4)
     &     +(ss12*ss2a+ss11*cs2a)*ett(izrg,irg,2)
     &     +(ds31*csa+ds23*ssa)*ett(izrg,irg,3)
        des=(ss12*cs2a-ss11*ss2a)*ept(izrg,irg,2)
     &     +(ds31*ssa-ds23*csa)*ept(izrg,irg,3)
c
        dze=(ss12*cs2a-ss11*ss2a)*erp(izrg,irg,2)
     &     +(ds31*ssa-ds23*csa)*erp(izrg,irg,3)
        dse=(ss12*cs2a-ss11*ss2a)*etp(izrg,irg,2)
     &     +(ds31*ssa-ds23*csa)*etp(izrg,irg,3)
        dee=expl*epp(izrg,irg,1)+clvd*epp(izrg,irg,4)
     &     +(ss12*ss2a+ss11*cs2a)*epp(izrg,irg,2)
     &     +(ds31*csa+ds23*ssa)*epp(izrg,irg,3)
c
c       4. Ezz, 5. Ezt, 6. Ezp, 7. Ett, 8. Etp, 9. Epp, 10. Roz, 11. Rot, 12. Rop
c
        y(4,izrg)=dzz
        y(5,izrg)=0.5d0*(dsz+dzs)
        y(6,izrg)=0.5d0*(dez+dze)
        y(7,izrg)=dss
        y(8,izrg)=0.5d0*(dse+des)
        y(9,izrg)=dee
c
        y(10,izrg)=0.5d0*(des-dse)
        y(11,izrg)=0.5d0*(dez-dze)
        y(12,izrg)=0.5d0*(dzs-dsz)
c
c       13-15. gravity (space fixed without free-air gradient effect)
c       16. geopotential
c
        y(13,izrg)=expl*gr(izrg,irg,1)+clvd*gr(izrg,irg,4)
     &      +(ss12*ss2a+ss11*cs2a)*gr(izrg,irg,2)
     &      +(ds31*csa+ds23*ssa)*gr(izrg,irg,3)
        y(14,izrg)=expl*gt(izrg,irg,1)+clvd*gt(izrg,irg,4)
     &          +(ss12*ss2a+ss11*cs2a)*gt(izrg,irg,2)
     &          +(ds31*csa+ds23*ssa)*gt(izrg,irg,3)
        y(15,izrg)=(ss12*cs2a-ss11*ss2a)*gp(izrg,irg,2)
     &          +(ds31*ssa-ds23*csa)*gp(izrg,irg,3)
        y(16,izrg)=expl*po(izrg,irg,1)+clvd*po(izrg,irg,4)
     &      +(ss12*ss2a+ss11*cs2a)*po(izrg,irg,2)
     &      +(ds31*csa+ds23*ssa)*po(izrg,irg,3)
      enddo
      return
      end
