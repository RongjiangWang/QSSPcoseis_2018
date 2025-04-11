	subroutine qpcsublayer(ierr)
      use qpcalloc
	implicit none
      integer*4 ierr
c
c	work space
c
	integer*4 i,j,i0,l,ly,lyz,ig,iz
      real*8 f,h,dh,z,zz,slw,wvlcut,up,lw,uplw4
	real*8 xup,xlw,rrs,rrr,dkap,dmue,rho1,drho,deta1,deta2
      real*8 mass,gr0,gr1,gr2,rrcm
      real*8 rr1,rr2,swapup(4),swaplw(4)
      logical jump
c
      lymax=0
      do l=1,l0-1
	  h=dp0lw(l)-dp0up(l)
	  dkap=2.d0*dabs(kap0lw(l)-kap0up(l))/(kap0lw(l)+kap0up(l))
        if(mue0lw(l)+mue0up(l).gt.0.d0)then
	    dmue=2.d0*dabs(mue0lw(l)-mue0up(l))/(mue0lw(l)+mue0up(l))
        else
          dmue=0.d0
        endif
        drho=2.d0*dabs(rho0lw(l)-rho0up(l))/(rho0lw(l)+rho0up(l))
        i0=1+idint(dmax1(dkap/RESOLUT,dmue/RESOLUT,drho/RESOLUT))
        lymax=lymax+i0
      enddo
      lymax=lymax+1+nz+3*ngrn
c
      allocate(izrly(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: izrly not allocated!'
      allocate(ldegpsv(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ldegpsv not allocated!'
      allocate(ldegsh(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ldegsh not allocated!'
      allocate(ndruku(0:ldegmax,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ndruku not allocated!'
      allocate(rrup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: rrup not allocated!'
      allocate(rrlw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: rrlw not allocated!'
      allocate(kapup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: kapup not allocated!'
      allocate(kaplw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: kaplw not allocated!'
      allocate(mueup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: mueup not allocated!'
      allocate(muelw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: muelw not allocated!'
      allocate(rhoup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: rhoup not allocated!'
      allocate(rholw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: rholw not allocated!'
      allocate(grup(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: grup not allocated!'
      allocate(grlw(lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: grlw not allocated!'
c
      allocate(mat2x2(2,2,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: mat2x2 not allocated!'
      allocate(mat2x2inv(2,2,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: '
     &                //'mat2x2inv not allocated!'
      allocate(mas3x3(3,3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: mas3x3 not allocated!'
      allocate(mas3x3inv(3,3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: '
     &                //'mas3x3inv not allocated!'
      allocate(mas4x4(4,4,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: mas4x4 not allocated!'
      allocate(mas4x4inv(4,4,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: '
     &                //'mas4x4inv not allocated!'
      allocate(mas6x6(6,6,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: mas6x6 not allocated!'
      allocate(mas6x6inv(6,6,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: '
     &                //'mas6x6inv not allocated!'
c
      allocate(cypnorm(6,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: cynorm not allocated!'
c
      allocate(yup2(2,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: yup2 not allocated!'
      allocate(ylw2(2,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ylw2 not allocated!'
      allocate(yup3(3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: yup3 not allocated!'
      allocate(ylw3(3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ylw3 not allocated!'
      allocate(yup6(6,3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: yup6 not allocated!'
      allocate(ylw6(6,3,lymax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpcsublayer: ylw2 not allocated!'
c
	ly0=0
      zz=0.d0
      gr2=gravity
	do l=1,l0-1
c
c       modification using Adam-Williamson condition
c
        gr1=gr2
        rr1=rearth-dp0up(l)
        rr2=rearth-dp0lw(l)
        drho=(rho0up(l)-rho0lw(l))/(rr1-rr2)
        rho1=rho0lw(l)-drho*rr2
        mass=PI*(rr1-rr2)
     &      *((4.d0/3.d0)*rho1*(rr1**2+rr1*rr2+rr2**2)
     &        +drho*(rr1**3+rr1**2*rr2+rr1*rr2**2+rr2**3))
        gr2=(gr1*rr1**2-BIGG*mass)/rr2**2
        if(gr2.le.0.d0)then
          print *,' Error in qpcsublayer:'
     &          //' density profile inconsistent with surface gravity!'
          stop
        endif
c
	  h=dp0lw(l)-dp0up(l)
	  dkap=2.d0*dabs(kap0lw(l)-kap0up(l))/(kap0lw(l)+kap0up(l))
        if(mue0lw(l)+mue0up(l).gt.0.d0)then
	    dmue=2.d0*dabs(mue0lw(l)-mue0up(l))/(mue0lw(l)+mue0up(l))
        else
          dmue=0.d0
        endif
        drho=2.d0*dabs(rho0lw(l)-rho0up(l))/(rho0lw(l)+rho0up(l))
        i0=1+idint(dmax1(dkap/RESOLUT,dmue/RESOLUT,drho/RESOLUT))
        dkap=(kap0lw(l)-kap0up(l))/h
	  dmue=(mue0lw(l)-mue0up(l))/h
        drho=(rho0lw(l)-rho0up(l))/h
	  dh=h/dble(i0)
	  do i=1,i0
	    ly0=ly0+1
	    z=dble(i-1)*dh
	    rrup(ly0)=rearth-(zz+z)
          if(i.eq.1)then
            rhoup(ly0)=rho0up(l)
	      kapup(ly0)=kap0up(l)
	      mueup(ly0)=mue0up(l)
          else
            rhoup(ly0)=rholw(ly0-1)
	      kapup(ly0)=kaplw(ly0-1)
	      mueup(ly0)=muelw(ly0-1)
          endif
          z=z+dh
	    rrlw(ly0)=rearth-(zz+z)
          rholw(ly0)=rho0up(l)+drho*z
	    kaplw(ly0)=kap0up(l)+dkap*z
	    muelw(ly0)=mue0up(l)+dmue*z
	  enddo
        zz=zz+h
      enddo
c
c     lowest layer assumed to be a homogeneous sphere
c
      ly0=ly0+1
      rrup(ly0)=rearth-dp0up(l0)
      rhoup(ly0)=gr2/(4.d0*PI*BIGG*rrup(ly0)/3.d0)
	kapup(ly0)=0.5d0*(kap0up(l0)+kap0lw(l0))
	mueup(ly0)=0.5d0*(mue0up(l0)+mue0lw(l0))
c
      rrlw(ly0)=0.d0
      rholw(ly0)=rhoup(ly0)
      kaplw(ly0)=kapup(ly0)
      muelw(ly0)=mueup(ly0)
c
      rrcm=0.d0
      do ly=2,ly0
        if(mueup(ly).le.0.d0.and.muelw(ly-1).gt.0.d0)then
          rrcm=rrup(ly)
          goto 10
        endif
      enddo
10    continue
c
c     add source layers
c
      do ig=1,ngrn
        rrs=rearth-grndep(ig)
        if(rrs.lt.rrcm.or.rrs.gt.rearth)then
          stop ' Error in qpcsublayer: bad source depth!'
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 100
          endif
        enddo
100     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      kapup(ly+1)=kapup(ly)
	      mueup(ly+1)=mueup(ly)
	      rhoup(ly+1)=rhoup(ly)
c
            rrlw(ly+1)=rrlw(ly)
	      kaplw(ly+1)=kaplw(ly)
	      muelw(ly+1)=muelw(ly)
	      rholw(ly+1)=rholw(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
	    rholw(lys-1)=up*rhoup(lys-1)+lw*rholw(lys)
	    kaplw(lys-1)=up*kapup(lys-1)+lw*kaplw(lys)
	    muelw(lys-1)=up*mueup(lys-1)+lw*muelw(lys)
c
          rrup(lys)=rrs
	    kapup(lys)=kaplw(lys-1)
	    mueup(lys)=muelw(lys-1)
	    rhoup(lys)=rholw(lys-1)
c
          ly0=ly0+1       
        endif
c
        if(gdds(ig).le.0.d0)goto 104
c
        rrs=rearth-(grndep(ig)-0.5d0*gdds(ig))
        if(rrs.gt.rearth)then
          rrs=rearth
          goto 102
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 101
          endif
        enddo
101     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      kapup(ly+1)=kapup(ly)
	      mueup(ly+1)=mueup(ly)
	      rhoup(ly+1)=rhoup(ly)
c
            rrlw(ly+1)=rrlw(ly)
	      kaplw(ly+1)=kaplw(ly)
	      muelw(ly+1)=muelw(ly)
	      rholw(ly+1)=rholw(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
	    rholw(lys-1)=up*rhoup(lys-1)+lw*rholw(lys)
	    kaplw(lys-1)=up*kapup(lys-1)+lw*kaplw(lys)
	    muelw(lys-1)=up*mueup(lys-1)+lw*muelw(lys)
c
          rrup(lys)=rrs
	    kapup(lys)=kaplw(lys-1)
	    mueup(lys)=muelw(lys-1)
	    rhoup(lys)=rholw(lys-1)
c
          ly0=ly0+1       
        endif
c
102     continue
c
        rrs=rearth-(grndep(ig)+0.5d0*gdds(ig))
        if(rrs.lt.rrcm)then
          rrs=rrcm
          goto 104
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 103
          endif
        enddo
103     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      kapup(ly+1)=kapup(ly)
	      mueup(ly+1)=mueup(ly)
	      rhoup(ly+1)=rhoup(ly)
c
            rrlw(ly+1)=rrlw(ly)
	      kaplw(ly+1)=kaplw(ly)
	      muelw(ly+1)=muelw(ly)
	      rholw(ly+1)=rholw(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
	    rholw(lys-1)=up*rhoup(lys-1)+lw*rholw(lys)
	    kaplw(lys-1)=up*kapup(lys-1)+lw*kaplw(lys)
	    muelw(lys-1)=up*mueup(lys-1)+lw*muelw(lys)
c
          rrup(lys)=rrs
	    kapup(lys)=kaplw(lys-1)
	    mueup(lys)=muelw(lys-1)
	    rhoup(lys)=rholw(lys-1)
c
          ly0=ly0+1       
        endif
104     continue
      enddo
c
c     add receiver layers
c
      do iz=1,nz
        if(depr(iz).lt.0.d0)then
          rrr=rearth-depr(iz)
          do ly=1,ly0
            if(rrr.ge.rrlw(ly))then
              lyz=ly
              goto 200
            endif
          enddo
200       continue
          if(rrr.lt.rrup(lyz).and.rrr.gt.rrlw(lyz))then
            do ly=ly0,lyz,-1
              rrup(ly+1)=rrup(ly)
	        kapup(ly+1)=kapup(ly)
	        mueup(ly+1)=mueup(ly)
	        rhoup(ly+1)=rhoup(ly)
c
              rrlw(ly+1)=rrlw(ly)
	        kaplw(ly+1)=kaplw(ly)
	        muelw(ly+1)=muelw(ly)
	        rholw(ly+1)=rholw(ly)
            enddo
            lyz=lyz+1
            up=(rrr-rrlw(lyz))/(rrup(lyz-1)-rrlw(lyz))
            lw=1.d0-up
            rrlw(lyz-1)=rrr
	      rholw(lyz-1)=up*rhoup(lyz-1)+lw*rholw(lyz)
	      kaplw(lyz-1)=up*kapup(lyz-1)+lw*kaplw(lyz)
            muelw(lyz-1)=up*mueup(lyz-1)+lw*muelw(lyz)
c
            rrup(lyz)=rrr
	      kapup(lyz)=kaplw(lyz-1)
	      mueup(lyz)=muelw(lyz-1)
	      rhoup(lyz)=rholw(lyz-1)
            ly0=ly0+1      
          endif
        endif
      enddo
c
c     determine indices of main interfaces
c
      lycm=ly0+1
      do ly=2,ly0
        if(mueup(ly).le.0.d0.and.muelw(ly-1).gt.0.d0)then
          lycm=ly
          goto 300
        endif
      enddo
300   lycc=ly0+1
      do ly=max0(2,lycm+1),ly0
        if(mueup(ly).gt.0.d0.and.muelw(ly-1).le.0.d0)then
          lycc=ly
          goto 400
        endif
      enddo
400   continue
c
c     determine indices of receiver layer
c     izrly(ly): index iz of receiver located at the ly-th layer
c     lyr(iz): index of layer where the iz-th receiver is located
c
      do ly=1,ly0
        izrly(ly)=0
      enddo
c
      do iz=1,nz
        if(depr(iz).lt.0.d0)then
c
c        for receivers outside the earth
c
          lyr(iz)=0
        else
          rrr=rearth-depr(iz)
          j=0
          do ly=1,ly0
            if(rrr.ge.rrup(ly))then
              j=ly
              goto 500
            endif
          enddo
500       continue
          if(j.eq.0)then
            stop ' Error in qpcsublayer: wrong receiver layer number!'
          endif
          izrly(j)=iz
          lyr(iz)=j
        endif
      enddo
c
c     determine indices of source layers
c
      do ig=1,ngrn
        rrs=rearth-grndep(ig)
        lygrn(ig)=ly0+1
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn(ig)=ly
            if(lygrn(ig).ge.lycm)then
              stop ' Error in qpssublayer: source below the mantle!'
            endif
            goto 600
          endif
        enddo
600     continue
c
        rrs=dmin1(rearth,rearth-(grndep(ig)-0.5d0*gdds(ig)))
        lygrn1(ig)=ly0+1
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn1(ig)=ly
            goto 601
          endif
        enddo
601     continue
c
        rrs=dmax1(rrcm,rearth-(grndep(ig)+0.5d0*gdds(ig)))
        lygrn2(ig)=ly0+1
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn2(ig)=ly
            goto 602
          endif
        enddo
602     continue
      enddo
c
      grup(1)=gravity
      do ly=1,ly0-1
        if(ly.gt.1)then
          grup(ly)=grlw(ly-1)
        endif
        drho=(rhoup(ly)-rholw(ly))/(rrup(ly)-rrlw(ly))
        rho1=rholw(ly)-drho*rrlw(ly) 
        mass=PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*rho1
     &      *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &      +drho*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &      +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
        grlw(ly)=(grup(ly)*rrup(ly)**2-BIGG*mass)/rrlw(ly)**2
      enddo
      if(ly0.gt.1)then
        grup(ly0)=grlw(ly0-1)
      endif
      grlw(ly0)=0.d0
c
      do ly=1,lycm-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(5,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(6,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
      do ly=lycm,lycc-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=(1.d0,0.d0)
        cypnorm(3,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(4,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
      do ly=lycc,ly0
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=dcmplx(kapup(ly)+mueup(ly)*4.d0/3.d0,0.d0)
        cypnorm(5,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
        cypnorm(6,ly)=dcmplx(4.d0*PI*BIGG*Rhoup(ly),0.d0)
      enddo
c
	write(*,'(9a)')'  No','       R[km]',' rho[g/cm^3]',
     &    '   kappa[Pa]','     mue[Pa]','    g[m/s^2]'
c
      rrlw(ly0)=0.d0
c
      i=0
	do ly=1,ly0
        write(*,1001)ly,rrup(ly)/KM2M,rhoup(ly)/KM2M,
     &               kapup(ly),mueup(ly),grup(ly)
        j=0
        do ig=1,ngrn
          if(lygrn(ig).eq.ly)then
            j=j+1
            write(*,'(a3,$)')' S '
          endif
        enddo
        if(j.eq.0)write(*,'(a3,$)')'   '
        if(izrly(ly).gt.0)then
          write(*,'(a3)')' R '
        else
          write(*,'(a3)')'   '
        endif
	  write(*,1002)rrlw(ly)/1.d3,rholw(ly)/1.d3,
     &               kaplw(ly),muelw(ly),grlw(ly)
      enddo
c
1000  format(i4,3f10.3,f11.3)
1001	format(i4,2f12.4,2E12.4,f9.4,$)
1002  format(f16.4,f12.4,2E12.4,f9.4)
      return
	end
