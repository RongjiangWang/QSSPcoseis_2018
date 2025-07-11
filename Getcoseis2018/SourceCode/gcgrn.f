      subroutine gcgrn(spatgrndir,ns,lats,lons,deps,spa,
     &                 nr,latr,lonr,depr,obs)
c
c     related subroutines:
c     disazi.f, gcalloc.f, gccmp.f, gcread.f
c
      use gcalloc
      implicit none
c
c     input parameters:
c     spatgrndir = directory of Green's function database
c     ns: no of sources
c     lats,lons,deps: source locations
c     vol=spa(0): anelastic strain volume
c     mee,men,mez,mnn,mnz,mzz=spa(1-6): moment tensor if vol = 0, else anelastic strain tensor.
c     nr: no of observation locations
c     latr,lonr,depr: receiver locations
c
      integer*4 ns,nr
      integer*4 nrg,nzrg,ng
      real*8 lats(ns),lons(ns),deps(ns)
      real*8 spa(0:6,ns)
      real*8 latr(nr),lonr(nr),depr(nr)
      character*80 spatgrndir
c
c     parameters to be returned:
c     1-3: Ue, Un, Uz (displacement);
c     4-9: Eee,Een,Eez,Enn,Enz,Ezz (strain);
c     10-12: Roe,Ron,Roz (rotation);
c     13-18: Eee,Een,Eez,Enn,Enz,Ezz (stress);
c     19-21: ge,gn,gz (space-based gravity change, gz downward positive)
c     22: geoid
c
      real*8 obs(22,nr)
c
c     local variables
c
      integer*4 i,j,k,ks,ierr,is,ir,iz,ig,flen,unit
      integer*4 igsel,offset,nrg1,nrg2
      integer*4 izr,izrg1,izrg2,irg,irg1,irg2
      integer*4 izrg,nzrg1,nzrg2,iglast
      real*8 rn,re,z1,z2,delta,dismax,dismin,deprmin,deprmax
      real*8 expl,clvd,ss12,ss11,ds31,ds23,depsmin,depsmax
      real*8 csa,ssa,cs2a,ss2a
      real*8 ssb,csb,yr,yt,yz,eii,lamr,muer,lams,mues
      real*8 mtt,mpp,mrr,mtp,mpr,mrt
      real*8 wr1,wr2,wz1,wz2
      real*8 rtz(3,3),enz(3,3),rot(3,3),swp(3,3)
c
      real*8, allocatable:: sum1(:,:),sum2(:,:),sumg(:,:)
c
      do flen=80,1,-1
        if(spatgrndir(flen:flen).ne.' ')goto 100
      enddo
      do i=1,flen
        if(spatgrndir(i:i).eq.'/')spatgrndir(i:i)='\'
      enddo
100   if(spatgrndir(flen:flen).ne.'\')then
        spatgrndir=spatgrndir(1:flen)//'\'
        flen=flen+1
        if(flen.gt.80)then
          stop ' Error in gcgrn: too long name of GF database!'
        endif 
      endif
c
      unit=10
      open(unit,file=spatgrndir(1:flen)//'GreenInfo.dat',status='old')
      call skipdoc(unit)
      read(unit,*)nrg
c
      allocate(rrg(nrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rrg not allocated!'
c
      call skipdoc(unit)
      read(unit,*)(rrg(i),i=1,nrg)
      do irg=1,nrg
        rrg(irg)=rrg(irg)*km2m
      enddo
c
      call skipdoc(unit)
      read(unit,*)nzrg
c
      allocate(zrg(nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: zrg not allocated!'
      allocate(vprg(nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vprg not allocated!'
      allocate(vsrg(nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vsrg not allocated!'
      allocate(rhorg(nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rhorg not allocated!'
c
      call skipdoc(unit)
      read(unit,*)((zrg(izrg),vprg(izrg),vsrg(izrg),rhorg(izrg)),
     &             izrg=1,nzrg)
      do izrg=1,nzrg
        zrg(izrg)=zrg(izrg)*km2m
        vprg(izrg)=vprg(izrg)*km2m
        vsrg(izrg)=vsrg(izrg)*km2m
        rhorg(izrg)=rhorg(izrg)*km2m
      enddo
c
      call skipdoc(unit)
      read(unit,*)ng
c
      allocate(zg(ng),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: zg not allocated!'
      allocate(vpg(ng),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vpg not allocated!'
      allocate(vsg(ng),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vsg not allocated!'
      allocate(rhog(ng),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rhog not allocated!'
      allocate(grnfile(ng),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: grnfile not allocated!'
c
      allocate(ur(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ur not allocated!'
      allocate(ut(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ut not allocated!'
      allocate(up(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: up not allocated!'
c
      allocate(gr(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gr not allocated!'
      allocate(gt(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gt not allocated!'
      allocate(gp(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gp not allocated!'
      allocate(po(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: po not allocated!'
c
      allocate(err(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: err not allocated!'
      allocate(ert(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ert not allocated!'
      allocate(erp(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: erp not allocated!'
      allocate(etr(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: etr not allocated!'
      allocate(ett(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ett not allocated!'
      allocate(etp(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: etp not allocated!'
      allocate(epr(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: epr not allocated!'
      allocate(ept(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ept not allocated!'
      allocate(epp(nzrg,nrg,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: epp not allocated!'
c
      call skipdoc(unit)
      read(unit,*)((zg(ig),vpg(ig),vsg(ig),rhog(ig),
     &             grnfile(ig)),ig=1,ng)
      do ig=1,ng
        zg(ig)=zg(ig)*km2m
        vpg(ig)=vpg(ig)*km2m
        vsg(ig)=vsg(ig)*km2m
        rhog(ig)=rhog(ig)*km2m
      enddo
c
      close(unit)
c
      z1=depr(1)
      z2=depr(1)
      do ir=2,nr
        z1=dmin1(z1,depr(ir))
        z2=dmax1(z2,depr(ir))
      enddo
c
      nzrg1=1
      do izr=1,nzrg
        if(z1.ge.zrg(izr))nzrg1=izr
      enddo
c
      nzrg2=nzrg
      do izr=nzrg,nzrg1,-1
        if(z2.le.zrg(izr))nzrg2=izr
      enddo
c
      allocate(sum1(22,nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sum1 not allocated!'
      allocate(sum2(22,nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sum2 not allocated!'
      allocate(sumg(22,nzrg),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sumg not allocated!'
      do izr=nzrg1,nzrg2
        do i=1,22
          sum1(i,izr)=0.d0
          sum2(i,izr)=0.d0
          sumg(i,izr)=0.d0
        enddo
      enddo
c
      allocate(dis(nr,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: dis not allocated!'
      allocate(azi(nr,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: azi not allocated!'
      allocate(bazi(nr,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: bazi not allocated!'
      allocate(next(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: next not allocated!'
c
      dismax=0.d0
      dismin=rearth*pi
      depsmax=deps(1)
      do is=1,ns
        next(is)=.true.
        depsmax=dmax1(depsmax,deps(is))
        do ir=1,nr
          call disazi(rearth,lats(is),lons(is),latr(ir),lonr(ir),rn,re)
c
c         dis = epicentral distance
c         azi = azimuth (from south to east) of vector source to receiver
c
          dis(ir,is)=dsqrt(rn*rn+re*re)
c
          dismax=dmax1(dismax,dis(ir,is))
          dismin=dmin1(dismin,dis(ir,is))
c
          if(dis(ir,is).gt.0.d0)then
            azi(ir,is)=datan2(re,-rn)
          else
c
c           assume southern receiver in case of 0 distance
c
            azi(ir,is)=datan2(0.d0,1.d0)
          endif
c
          call disazi(rearth,latr(ir),lonr(ir),lats(is),lons(is),rn,re)
c
          if(dis(ir,is).gt.0.d0)then
c
c           azimuth (from south to east) of opposite vector receiver to source
c
            bazi(ir,is)=datan2(-re,rn)
          else
            bazi(ir,is)=datan2(0.d0,1.d0)
          endif
        enddo
      enddo
c
      deprmax=depr(1)
      deprmin=depr(1)
      do ir=2,nr
        deprmax=dmax1(deprmax,depr(ir))
        deprmin=dmin1(deprmin,depr(ir))
      enddo
c
      if(dismin.lt.rrg(1).or.dismax.gt.rrg(nrg))then
        print *,' Warning in gcgrn: '
     &     //'insufficient distance coverage of Green functions!'
      endif
      if(deps(1).lt.zg(1).or.deps(ns).gt.zg(ng))then
        print *,' Warning in gcgrn: '
     &     //'insufficient source depth coverage of Green functions!'
      endif
      if(deprmin.lt.zrg(1).or.deprmax.gt.zrg(nzrg))then
        print *,' Warning in gcgrn: '
     &     //'insufficient receiver depth coverage of Green functions!'
      endif
c
      nrg1=1
      do irg=2,nrg
        if(dismin.ge.rrg(irg))nrg1=irg
      enddo
c
      nrg2=nrg
      do irg=nrg-1,1,-1
        if(dismax.le.rrg(irg))nrg2=irg
      enddo
c
      iglast=0
c
      do ir=1,nr
        do i=1,22
          obs(i,ir)=0.d0
        enddo
      enddo
c
      do ks=1,ns
        depsmin=depsmax
        k=0
        do is=1,ns
          if(next(is).and.depsmin.ge.deps(is))then
            k=is
            depsmin=deps(is)
          endif
        enddo
        is=k
        next(is)=.false.
c
        ig=1
        delta=dabs(deps(is)-zg(1))
        do i=2,ng
          if(delta.gt.dabs(deps(is)-zg(i)))then
            delta=dabs(deps(is)-zg(i))
            ig=i
          endif
        enddo
c
        if(dabs(spa(0,is)).gt.0.d0)then
          lams=rhog(ig)*(vpg(ig)**2-2.d0*vsg(ig)**2)
          mues=rhog(ig)*vsg(ig)**2
          eii=spa(1,is)+spa(4,is)+spa(6,is)
          spa(1,is)=spa(0,is)*(lams*eii+2.d0*mues*spa(1,is))
          spa(2,is)=spa(0,is)*2.d0*mues*spa(2,is)
          spa(3,is)=spa(0,is)*2.d0*mues*spa(3,is)
          spa(4,is)=spa(0,is)*(lams*eii+2.d0*mues*spa(4,is))
          spa(5,is)=spa(0,is)*2.d0*mues*spa(5,is)
          spa(6,is)=spa(0,is)*(lams*eii+2.d0*mues*spa(6,is))
        endif
c
        mpp=spa(1,is)
        mtp=-spa(2,is)
        mpr=spa(3,is)
        mtt=spa(4,is)
        mrt=-spa(5,is)
        mrr=spa(6,is)
c
        expl=(mtt+mpp+mrr)/3.d0
        clvd=mrr-expl
        ss12=mtp
        ss11=(mtt-mpp)/2.d0
        ds31=mrt
        ds23=mpr
c
        if(ig.gt.iglast)then
          iglast=ig
          unit=20+ig
          open(unit,file=grnfile(ig),form='unformatted',status='old')
c
c         each of nrg distances has 52 records
c         each record includes nzrg depths (= nzrg*8 byes) plus 8 byes alignment
c
          offset=(nrg1-1)*52*(nzrg+1)*8
          call fseek(unit,offset,0,i)
          call gcread(unit,nrg1,nrg2,nzrg)
          close(unit)
        endif
c
        do ir=1,nr
          ssa=dsin(azi(ir,is))
          csa=dcos(azi(ir,is))
          ss2a=dsin(2.d0*azi(ir,is))
          cs2a=dcos(2.d0*azi(ir,is))
c
          ssb=dsin(bazi(ir,is))
          csb=dcos(bazi(ir,is))
c
          rot(1,1)=ssb
          rot(1,2)=csb
          rot(1,3)=0.d0
          rot(2,1)=-csb
          rot(2,2)=ssb
          rot(2,3)=0.d0
          rot(3,1)=0.d0
          rot(3,2)=0.d0
          rot(3,3)=1.d0
c
          irg1=0
          do irg=1,nrg
            if(dis(ir,is).ge.rrg(irg))irg1=irg
          enddo
          if(irg1.eq.0)then
            irg1=1
            irg2=1
          else if(irg1.eq.nrg)then
            irg2=irg1
          else
            irg2=min0(irg1+1,nrg)
          endif
c
          if(irg2.gt.irg1)then
            wr1=(rrg(irg2)-dis(ir,is))/(rrg(irg2)-rrg(irg1))
            wr2=1.d0-wr1
          else
            wr1=1.d0
            wr2=0.d0
          endif
c
          izrg1=nzrg1
          do izrg=nzrg1,nzrg2
            if(depr(ir).ge.zrg(izrg))izrg1=izrg
          enddo
          izrg2=min0(izrg1+1,nzrg2)
          if(izrg2.gt.izrg1)then
            wz1=(zrg(izrg2)-depr(ir))/(zrg(izrg2)-zrg(izrg1))
            wz2=1.d0-wz1
          else
            wz1=1.d0
            wz2=0.d0
          endif
c
          if(dabs(wr1).gt.0.d0)then
            call gccmp(expl,clvd,ss12,ss11,ds31,ds23,csa,ssa,cs2a,ss2a,
     &                 irg1,izrg1,izrg2,nzrg,sum1)
          endif
c
          if(dabs(wr2).gt.0.d0)then
            call gccmp(expl,clvd,ss12,ss11,ds31,ds23,csa,ssa,cs2a,ss2a,
     &                 irg2,izrg1,izrg2,nzrg,sum2)
          endif
c
          do izrg=izrg1,izrg2
            do i=1,16
              sumg(i,izrg)=wr1*sum1(i,izrg)+wr2*sum2(i,izrg)
            enddo
c
            sumg(19,izrg)=sumg(13,izrg)
            sumg(20,izrg)=sumg(14,izrg)
            sumg(21,izrg)=sumg(15,izrg)
            sumg(22,izrg)=sumg(16,izrg)
c
c           transform of displacement vector
c
            yz=sumg(1,izrg)
            yr=sumg(2,izrg)
            yt=sumg(3,izrg)
c
c           1/2/3 = E/N/Z
c
            sumg(1,izrg)= yr*ssb+yt*csb
            sumg(2,izrg)=-yr*csb+yt*ssb
            sumg(3,izrg)=yz
c
c           transform of strain tensor
c
            rtz(1,1)=sumg(7,izrg)
            rtz(1,2)=sumg(8,izrg)
            rtz(1,3)=sumg(5,izrg)
            rtz(2,1)=sumg(8,izrg)
            rtz(2,2)=sumg(9,izrg)
            rtz(2,3)=sumg(6,izrg)
            rtz(3,1)=sumg(5,izrg)
            rtz(3,2)=sumg(6,izrg)
            rtz(3,3)=sumg(4,izrg)
c
            do i=1,3
              do j=1,3
                swp(i,j)=0.d0
                do k=1,3
                  swp(i,j)=swp(i,j)+rot(i,k)*rtz(k,j)
                enddo
              enddo
            enddo
            do i=1,3
              do j=1,3
                enz(i,j)=0.d0
                do k=1,3
                  enz(i,j)=enz(i,j)+swp(i,k)*rot(j,k)
                enddo
              enddo
            enddo
            sumg(4,izrg)=enz(1,1)
            sumg(5,izrg)=enz(1,2)
            sumg(6,izrg)=enz(1,3)
            sumg(7,izrg)=enz(2,2)
            sumg(8,izrg)=enz(2,3)
            sumg(9,izrg)=enz(3,3)
c
c           transform of rotation vector
c
            yz=sumg(10,izrg)
            yr=sumg(11,izrg)
            yt=sumg(12,izrg)
c
            sumg(10,izrg)= yr*ssb+yt*csb
            sumg(11,izrg)=-yr*csb+yt*ssb
            sumg(12,izrg)=yz
c
c           calculate stress from strain
c
            lamr=rhorg(izrg)*(vprg(izrg)**2-2.d0*vsrg(izrg)**2)
            muer=rhorg(izrg)*vsrg(izrg)**2
c
            eii=sumg(4,izrg)+sumg(7,izrg)+sumg(9,izrg)
            sumg(13,izrg)=lamr*eii+2.d0*muer*sumg(4,izrg)
            sumg(14,izrg)=2.d0*muer*sumg(5,izrg)
            sumg(15,izrg)=2.d0*muer*sumg(6,izrg)
            sumg(16,izrg)=lamr*eii+2.d0*muer*sumg(7,izrg)
            sumg(17,izrg)=2.d0*muer*sumg(8,izrg)
            sumg(18,izrg)=lamr*eii+2.d0*muer*sumg(9,izrg)
c
c           transform of gravity vector
c           (note that yz = vertical component, downwards positive!)
c
            yz=sumg(19,izrg)
            yr=sumg(20,izrg)
            yt=sumg(21,izrg)
c
            sumg(19,izrg)= yr*ssb+yt*csb
            sumg(20,izrg)=-yr*csb+yt*ssb
            sumg(21,izrg)=yz
          enddo
c
          if(depr(ir).ge.0.d0)then
            do i=1,22
              obs(i,ir)=obs(i,ir)+wz1*sumg(i,izrg1)+wz2*sumg(i,izrg2)
            enddo
          else
            do i=19,22
              obs(i,ir)=obs(i,ir)+wz1*sumg(i,izrg1)+wz2*sumg(i,izrg2)
            enddo
          endif
        enddo
      enddo
c
      deallocate(rrg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rrg not deallocated!'
c
      deallocate(zrg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: zrg not deallocated!'
      deallocate(vprg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vprg not deallocated!'
      deallocate(vsrg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vsrg not deallocated!'
      deallocate(rhorg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rhorg not deallocated!'
c
      deallocate(zg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: zg not deallocated!'
      deallocate(vpg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vpg not deallocated!'
      deallocate(vsg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: vsg not deallocated!'
      deallocate(rhog,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: rhog not deallocated!'
      deallocate(grnfile,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: grnfile not deallocated!'
c
      deallocate(ur,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ur not deallocated!'
      deallocate(ut,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ut not deallocated!'
      deallocate(up,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: up not deallocated!'
c
      deallocate(gr,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gr not deallocated!'
      deallocate(gt,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gt not deallocated!'
      deallocate(gp,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: gp not deallocated!'
      deallocate(po,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: po not deallocated!'
c
      deallocate(err,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: err not deallocated!'
      deallocate(ert,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ert not deallocated!'
      deallocate(erp,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: erp not deallocated!'
      deallocate(etr,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: etr not deallocated!'
      deallocate(ett,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ett not deallocated!'
      deallocate(etp,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: etp not deallocated!'
      deallocate(epr,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: epr not deallocated!'
      deallocate(ept,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: ept not deallocated!'
      deallocate(epp,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: epp not deallocated!'
c
      deallocate(sum1,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sum1 not deallocated!'
      deallocate(sum2,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sum2 not deallocated!'
      deallocate(sumg,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: sumg not deallocated!'
c
      deallocate(dis,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: dis not deallocated!'
      deallocate(azi,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: azi not deallocated!'
      deallocate(bazi,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: bazi not deallocated!'
      deallocate(next,stat=ierr)
      if(ierr.ne.0)stop ' Error in gcgrn: next not deallocated!'
      return
      end
