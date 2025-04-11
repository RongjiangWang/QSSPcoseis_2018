      subroutine qpcgetinp(unit)
      use qpcalloc
      implicit none
      integer*4 unit
c
c     work space
c
      integer*4 i,j,k,l,flen,ig,istp,ierr
      real*8 wup,wlw,deprsmax,ddep,ddr,delta,vp0,vs0
      real*8 dswap(9)
      character*80 fswap
c
c     read folder of spectral Green's functions
c
      call skipdoc(unit)
      read(unit,*)specgrndir
c
c     read critical and cutoff harmonic degrees
c
      call skipdoc(unit)
      read(unit,*)ldeggr,ldegcut
      nogravity=ldeggr.le.0
c
      ldegmax=ldegcut+1+ndmax
      allocate(tap(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: tap not allocated!'
      allocate(disk(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: disk not allocated!'
c
      allocate(lyuppsv(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lyuppsv not allocated!'
      allocate(lylwpsv(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lylwpsv not allocated!'
      allocate(lyupt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lyupt not allocated!'
      allocate(lylwt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lylwt not allocated!'
c
c     read receiver depth sampling parameters
c
      call skipdoc(unit)
      read(unit,*)zr1,zr2,dzr
      if(zr1.gt.zr2.or.dzr.le.0.d0)then
        stop ' Error in qpcgetinp: bad receiver depth parameters!'
      endif
c
      zr1=zr1*KM2M
      zr2=zr2*KM2M
      dzr=dzr*KM2M
      nz=1+idnint((zr2-zr1)/dzr)
c
      allocate(lyr(nz),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lyr not allocated!'
      allocate(depr(nz),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: depr not allocated!'
c
      if(nz.eq.1)then
        depr(1)=zr1
      else if(nz.eq.2)then
        depr(1)=zr1
        depr(2)=zr2
      else
        dzr=(zr2-zr1)/dble(nz-1)
        depr(1)=zr1
        do i=2,nz-1
          depr(i)=zr1+dble(i-1)*dzr
        enddo
        depr(nz)=zr2
      endif
c
c     read number of Green's functions
c
      call skipdoc(unit)
      read(unit,*)ngrn
c
      allocate(grndep(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: grndep not allocated!'
      allocate(gdds(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: gdds not allocated!'
      allocate(grrs(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: grrs not allocated!'
      allocate(grnsel(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: grnsel not allocated!'
      allocate(spatgrnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: '
     &                //'spatgrnfile not allocated!'
      allocate(grnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: grnfile not allocated!'
      allocate(uspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: uspecfile not allocated!'
      allocate(vspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: vspecfile not allocated!'
      allocate(wspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: wspecfile not allocated!'
      allocate(especfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: especfile not allocated!'
      allocate(fspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: fspecfile not allocated!'
      allocate(gspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: gspecfile not allocated!'
      allocate(pspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: pspecfile not allocated!'
      allocate(qspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: qspecfile not allocated!'
      allocate(lygrn(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lygrn not allocated!'
      allocate(lygrn1(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lygrn1 not allocated!'
      allocate(lygrn2(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: lygrn2 not allocated!'
c
      do ig=1,ngrn
        call skipdoc(unit)
        read(unit,*)grndep(ig),gdds(ig),grrs(ig),grnfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' Error in qpcgetinp: bad selection of spectral GF!'
        endif
        grndep(ig)=grndep(ig)*KM2M
        gdds(ig)=dmax1(0.d0,(4.d0/3.d0)*gdds(ig)*KM2M)
        grrs(ig)=dmax1(0.d0,dsqrt(3.d0)*grrs(ig)*KM2M)
      enddo
c
c     sort spectral Green's function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            dswap(1)=grndep(i)
            dswap(2)=gdds(i)
            dswap(3)=grrs(i)
            fswap=grnfile(i)
            k=grnsel(i)
c
            grndep(i)=grndep(j)
            gdds(i)=gdds(j)
            grrs(i)=grrs(j)
            grnfile(i)=grnfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            gdds(j)=dswap(2)
            grrs(j)=dswap(3)
            grnfile(j)=fswap
            grnsel(j)=k
          endif
        enddo
      enddo
c
c     complete file names of spectral Green's functions
c
      do flen=80,1,-1
        if(specgrndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do i=1,flen
        if(specgrndir(i:i).eq.'/')specgrndir(i:i)='\'
      enddo
      do ig=1,ngrn
        uspecfile(ig)=specgrndir(1:flen)//'U_'//grnfile(ig)
        vspecfile(ig)=specgrndir(1:flen)//'V_'//grnfile(ig)
        wspecfile(ig)=specgrndir(1:flen)//'W_'//grnfile(ig)
        especfile(ig)=specgrndir(1:flen)//'E_'//grnfile(ig)
        fspecfile(ig)=specgrndir(1:flen)//'F_'//grnfile(ig)
        gspecfile(ig)=specgrndir(1:flen)//'G_'//grnfile(ig)
        pspecfile(ig)=specgrndir(1:flen)//'P_'//grnfile(ig)
        qspecfile(ig)=specgrndir(1:flen)//'Q_'//grnfile(ig)
      enddo
      specgrntmp=specgrndir(1:flen)//'specgrn.tmp'
c
c     read folder of spatial Green's functions
c
      call skipdoc(unit)
      read(unit,*)spatgrndir
c
c     complete file names of spatial Green's functions
c
      do flen=80,1,-1
        if(spatgrndir(flen:flen).ne.' ')goto 300
      enddo
300   continue
      do i=1,flen
        if(spatgrndir(i:i).eq.'/')spatgrndir(i:i)='\'
      enddo
      do ig=1,ngrn
        spatgrnfile(ig)=spatgrndir(1:flen)//grnfile(ig)
      enddo
c
c     read receiver locations relative to source
c
      call skipdoc(unit)
      read(unit,*)dr1,dr2,ddr1,ddr2
      if(dr1.lt.0.d0.or.dr1.gt.dr2.or.ddr1.le.0.d0.or.ddr1.gt.ddr2)then
        stop ' Error in qpcgetinp: bad receiver distance parameters!'
      endif
      dr1=dr1*KM2M
      dr2=dr2*KM2M
      ddr1=ddr1*KM2M
      ddr2=ddr2*KM2M
      nr=1+idnint(2.d0*(dr2-dr1)/(ddr1+ddr2))
c
      allocate(idr(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: idr not allocated!'
      allocate(dism(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: dism not allocated!'
      allocate(disrad(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: disrad not allocated!'
      allocate(ssd(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ssd not allocated!'
      allocate(csd(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: csd not allocated!'
      allocate(ssf(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ssf not allocated!'
c
      allocate(plm(0:ldegmax,0:2,nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: plm not allocated!'
c
      allocate(ul0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ul0 not allocated!'
      allocate(vl0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: vl0 not allocated!'
      allocate(wl0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: wl0 not allocated!'
      allocate(el0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: el0 not allocated!'
      allocate(fl0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: fl0 not allocated!'
      allocate(gl0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: gl0 not allocated!'
      allocate(pl0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: pl0 not allocated!'
      allocate(ql0(0:ldegmax,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ql0 not allocated!'
c
      allocate(urlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: urlm not allocated!'
      allocate(utlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: utlm not allocated!'
      allocate(uplm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: uplm not allocated!'
c
      allocate(grlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: grlm not allocated!'
      allocate(gtlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: gtlm not allocated!'
      allocate(gplm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: gplm not allocated!'
      allocate(polm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: polm not allocated!'
c
      allocate(errlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: errlm not allocated!'
      allocate(ertlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ertlm not allocated!'
      allocate(erplm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: erplm not allocated!'
      allocate(etrlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: etrlm not allocated!'
      allocate(ett0lm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ett0lm not allocated!'
      allocate(ettalm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ettalm not allocated!'
      allocate(ettblm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ettblm not allocated!'
      allocate(etp0lm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: etp0lm not allocated!'
      allocate(etpalm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: etpalm not allocated!'
      allocate(etpblm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: etpblm not allocated!'
      allocate(eprlm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: eprlm not allocated!'
      allocate(ept0lm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ept0lm not allocated!'
      allocate(eptalm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: eptalm not allocated!'
      allocate(eptblm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: eptblm not allocated!'
      allocate(epp0lm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: epp0lm not allocated!'
      allocate(eppalm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: eppalm not allocated!'
      allocate(eppblm(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: eppblm not allocated!'
c
      allocate(ypsv(6,4,nz),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ypsv not allocated!'
      allocate(ypsvg(6,4,nz),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ypsvg not allocated!'
      allocate(ysh(2,2,nz),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: ysh not allocated!'
c
      allocate(spatur(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatur not allocated!'
      allocate(spatut(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatut not allocated!'
      allocate(spatup(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatup not allocated!'
c
      allocate(spatgr(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatgr not allocated!'
      allocate(spatgt(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatgt not allocated!'
      allocate(spatgp(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatgp not allocated!'
      allocate(spatpo(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatpo not allocated!'
c
      allocate(spaterr(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spaterr not allocated!'
      allocate(spatert(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatert not allocated!'
      allocate(spaterp(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spaterp not allocated!'
      allocate(spatetr(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatetr not allocated!'
      allocate(spatett(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatett not allocated!'
      allocate(spatetp(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatetp not allocated!'
      allocate(spatepr(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatepr not allocated!'
      allocate(spatept(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatept not allocated!'
      allocate(spatepp(nz,nr,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: spatepp not allocated!'
c
      if(nr.eq.1)then
        dism(1)=dr1
      else if(nr.eq.2)then
        dism(1)=dr1
        dism(2)=dr2
      else
        ddr=(dr2-dr1)/dble(nr-1)
        ddr1=ddr1*2.d0*ddr/(ddr1+ddr2)
        ddr2=ddr2*2.d0*ddr/(ddr1+ddr2)
        dism(1)=dr1
        do i=2,nr-1
          ddr=ddr1+(ddr2-ddr1)*dble(i-2)/dble(nr-2)
          dism(i)=dism(i-1)+ddr
        enddo
        dism(nr)=dr2
      endif
c
c     read information file for spatial Green's functions
c
      call skipdoc(unit)
      read(unit,*)infofile
      infofile=spatgrndir(1:flen)//infofile
c
c     multilayered model parameters
c     =============================
c
      lyadd=nz+ngrn+1
      call skipdoc(unit)
      read(unit,*)l,rearth,gravity
      rearth=rearth*KM2M
c
      l0=l+1
      allocate(dp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: dp0 not allocated!'
      allocate(mue0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: mue0 not allocated!'
      allocate(kap0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: kap0 not allocated!'
      allocate(rho0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: rho0 not allocated!'
c
      allocate(dp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: dp0up not allocated!'
      allocate(mue0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: mue0up not allocated!'
      allocate(kap0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: kap0up not allocated!'
      allocate(rho0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: rho0up not allocated!'
c
      allocate(dp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: dp0lw not allocated!'
      allocate(mue0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: mue0lw not allocated!'
      allocate(kap0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: kap0lw not allocated!'
      allocate(rho0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcgetinp: rho0lw not allocated!'
c
      do i=1,nr
        disrad(i)=dism(i)/rearth
      enddo
c
      do i=1,l
        call skipdoc(unit)
        read(unit,*)j,dp0(i),vp0,vs0,rho0(i)
        if(vp0.le.0.d0.or.vs0.lt.0.d0.or.rho0(i).le.0.d0)then
          stop ' Error in qpcgetinp: bad seismic parameter!'
        endif
c
c       input units:    -,km,  km/s, km/s, g/cm^3
c
        dp0(i)=KM2M*dp0(i)
        vp0=KM2M*vp0
        vs0=KM2M*vs0
        rho0(i)=KM2M*rho0(i)
        mue0(i)=rho0(i)*vs0**2
        kap0(i)=rho0(i)*vp0**2-mue0(i)*4.d0/3.d0
      enddo
c
      if(dp0(1).lt.0.d0.or.dp0(1).gt.0.d0)then
        stop ' Error in qpcgetinp: wrong first interface depth!'
      else if(dp0(l).gt.rearth)then
        stop ' Error in qpcgetinp: too large interface depth!'
      else if(dp0(l).lt.rearth)then
        l=l+1
        dp0(l)=rearth
        kap0(l)=kap0(l-1)
        mue0(l)=mue0(l-1)
        rho0(l)=rho0(l-1)
      endif
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          kap0up(l0)=kap0(i-1)
          mue0up(l0)=mue0(i-1)
          rho0up(l0)=rho0(i-1)
c
          dp0lw(l0)=dp0(i)
          kap0lw(l0)=kap0(i)
          mue0lw(l0)=mue0(i)
          rho0lw(l0)=rho0(i)
        endif
      enddo
c
      deprsmax=0.d0
      do i=1,ngrn
        deprsmax=dmax1(deprsmax,grndep(i))
      enddo
      do i=1,nz
        deprsmax=dmax1(deprsmax,depr(i))
      enddo
c
      ddep=dp0lw(l0)-dp0up(l0)
      if(dp0up(l0).le.deprsmax+THICKMAX)then
        l0=l0+1
        if(l0.ge.lymax-lyadd)then
          stop ' Error in qpsgetinp: lymax defined too small!'
        endif
        if(deprsmax.ge.dp0up(l0-1))then
          delta=deprsmax-dp0up(l0-1)
     &         +dmin1(THICKMAX,0.5d0*(dp0up(l0-1)-deprsmax))
        else
          delta=dmin1(THICKMAX,0.5d0*ddep)
        endif
c
        wlw=delta/ddep
        wup=1.d0-wlw
c
        dp0lw(l0)=dp0lw(l0-1)
        kap0lw(l0)=kap0lw(l0-1)
        mue0lw(l0)=mue0lw(l0-1)
        rho0lw(l0)=rho0lw(l0-1)
c
        dp0up(l0)=dp0up(l0-1)+delta
        kap0up(l0)=wup*kap0up(l0-1)+wlw*kap0lw(l0-1)
        mue0up(l0)=wup*mue0up(l0-1)+wlw*mue0lw(l0-1)
        rho0up(l0)=wup*rho0up(l0-1)+wlw*rho0lw(l0-1)
c
        dp0lw(l0-1)=dp0up(l0)
        kap0lw(l0-1)=kap0up(l0)
        mue0lw(l0-1)=mue0up(l0)
        rho0lw(l0-1)=rho0up(l0)
      endif
c
c     end of inputs
c     =============
c
      return
      end
