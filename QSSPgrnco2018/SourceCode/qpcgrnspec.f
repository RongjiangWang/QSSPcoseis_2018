      subroutine qpcgrnspec(ig)
      use qpcalloc
      implicit none
      integer*4 ig
c
      integer*4 i,j,iz,il,istp,ly,lyrs,ldeg,ldegup,offset
      real*8 dll,expo,fac,ca,cb,cs1,cs2,cs3,cs4,ct1,ct2
      real*8 rrr,rrs,xs,a51,a55,a56,sum
c
      real*8 expos,diskcut
      data expos,diskcut/48.d0,1.0d-04/
c
      fac=rearth
      lyrs=lyr(1)
      do iz=1,nz
        if(depr(iz).ge.0.d0.and.fac.ge.dabs(depr(iz)-grndep(ig)))then
          fac=dabs(depr(iz)-grndep(ig))
          lyrs=lyr(iz)
        endif
      enddo
      fac=dsqrt(fac**2+(dism(1)+grrs(ig))**2)
      if(fac.le.0.d0)then
        ldegup=ldegcut
      else
        ldegup=min0(ldegcut,
     &              1000+idnint(10.d0*PI*rrup(lys)/fac))
      endif
c
      if(ldeggr.gt.0)then
        do ly=1,ly0
          do ldeg=0,ldegup
            ndruku(ldeg,ly)=10
          enddo
        enddo
      endif
c
      rrs=rrup(lys)
      gdds0=gdds(ig)
c
      if(lys2.gt.lys1)then
        mues=0.d0
        ksis=0.d0
        gddsnorm=0.d0
        do ly=lys1,lys-1
          ca=(rrup(lys)+0.5d0*gdds0-rrlw(ly))**3
          cb=(rrup(lys)+0.5d0*gdds0-rrup(ly))**3
          xs=(ca-cb)/3.d0
          gddsnorm=gddsnorm+xs
          mues=mues+0.5d0*(mueup(ly)+muelw(ly))*xs
          ksis=ksis+0.5d0*(kapup(ly)+kaplw(ly))*xs
        enddo
        do ly=lys,lys2-1
          ca=(rrup(ly)-rrup(lys)+0.5d0*gdds0)**3
          cb=(rrlw(ly)-rrup(lys)+0.5d0*gdds0)**3
          xs=(ca-cb)/3.d0
          gddsnorm=gddsnorm+xs
          mues=mues+0.5d0*(mueup(ly)+muelw(ly))*xs
          ksis=ksis+0.5d0*(kapup(ly)+kaplw(ly))*xs
        enddo
        mues=mues/gddsnorm
        ksis=ksis/gddsnorm+mues*4.d0/3.d0
      else
        mues=mueup(lys)
        ksis=kapup(lys)+mues*4.d0/3.d0
      endif
c
      if(grrs(ig).gt.0.d0)then
        xs=dcos(grrs(ig)/rrs)
        disk(0)=1.d0
        disk(1)=(3.d0+xs)/4.d0
        do ldeg=2,ldegup
          disk(ldeg)=(dble(2*ldeg-1)*xs*disk(ldeg-1)
     &               +dble(4-ldeg)*disk(ldeg-2))/dble(ldeg+3)
        enddo
        il=ldegup
        do ldeg=0,ldegup
          if(dabs(disk(ldeg)).gt.diskcut)il=ldeg
          disk(ldeg)=disk(ldeg)*dble(2*ldeg+1)/(4.d0*PI*rrs**2)
        enddo
        ldegup=il
      else
        do ldeg=0,ldegup
          disk(ldeg)=dble(2*ldeg+1)/(4.d0*PI*rrs**2)
        enddo
      endif
c
      write(*,*)' ==================================================='
      write(*,'(a,i3,a,f7.2,a)')'  Green functions for ',
     &        ig,'. source at depth ',grndep(ig)/KM2M,' km'
	write(*,'(a,i12)')'  cutoff harmonic degree = ',ldegup
      write(*,*)' ==================================================='
c
c     spherical coordinate system:
c     z = upward, s = co-latitude, e = longitude
c     put source at the north pole
c
      open(21,file=uspecfile(ig),form='unformatted',status='unknown')
      open(22,file=vspecfile(ig),form='unformatted',status='unknown')
      open(23,file=wspecfile(ig),form='unformatted',status='unknown')
      open(24,file=especfile(ig),form='unformatted',status='unknown')
      open(25,file=fspecfile(ig),form='unformatted',status='unknown')
      open(26,file=gspecfile(ig),form='unformatted',status='unknown')
      open(27,file=pspecfile(ig),form='unformatted',status='unknown')
      open(28,file=qspecfile(ig),form='unformatted',status='unknown')
      write(21)ldegup
      write(22)ldegup
      write(23)ldegup
      write(24)ldegup
      write(25)ldegup
      write(26)ldegup
      write(27)ldegup
      write(28)ldegup
c
      do istp=1,4
        do ldeg=0,1
          ul0(ldeg,istp)=0.d0
          vl0(ldeg,istp)=0.d0
          wl0(ldeg,istp)=0.d0
          el0(ldeg,istp)=0.d0
          fl0(ldeg,istp)=0.d0
          gl0(ldeg,istp)=0.d0
          ql0(ldeg,istp)=0.d0
          pl0(ldeg,istp)=0.d0
        enddo
      enddo
      if(.not.nogravity)then
        do ly=1,ly0
          do ldeg=0,ldegup
            ndruku(ldeg,ly)=0
          enddo
        enddo
      endif
c
      do ldeg=0,ldegup
        dll=dsqrt(dble(ldeg)*dble(ldeg+1))
c
c       determine degree dependent starting layer number
c       of sh solution
c
        lyupt(ldeg)=1
        expo=0.d0
        do ly=min0(lys1,lyr(1))-1,1,-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lyupt(ldeg)=ly
            goto 101
          endif
        enddo
101     continue
c
        lylwt(ldeg)=min0(lycm,ly0)
        expo=0.d0
        do ly=max0(lys2,lyr(nz))+1,min0(lycm,ly0)-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lylwt(ldeg)=ly
            goto 102
          endif
        enddo
102     continue
c
c       determine degree dependent starting layer number
c       of psv solution
c
        lyuppsv(ldeg)=1
        expo=0.d0
        do ly=min0(lys1,lyr(1))-1,1,-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lyuppsv(ldeg)=ly
            goto 201
          endif
        enddo
201     continue
c
        lylwpsv(ldeg)=ly0
        expo=0.d0
        do ly=max0(lys2,lyr(nz))+1,ly0-1
          expo=expo+dll*dlog(rrup(ly)/rrlw(ly))
          if(expo.gt.expos)then
            lylwpsv(ldeg)=ly
            goto 202
          endif
        enddo
202     continue
      enddo
c
c     determine layer dependent max. harmonic degree
c     of sh solution
c
      do ly=1,min0(lycm-1,ly0)
        ldegsh(ly)=1
        do ldeg=1,ldegup
          if(ly.ge.lyupt(ldeg).and.ly.le.lylwt(ldeg))then
            ldegsh(ly)=ldeg
          endif
        enddo
      enddo
c
c     determine layer dependent max. harmonic degree
c     of psv solution
c
      do ly=1,ly0
        ldegpsv(ly)=0
        do ldeg=0,ldegup
          if(ly.ge.min0(lyuppsv(ldeg),lyuppsv(ldeg)).and.
     &       ly.le.max0(lylwpsv(ldeg),lylwpsv(ldeg)))then
            ldegpsv(ly)=ldeg
          endif
        enddo
      enddo
c============================================================================
      open(30,file=specgrntmp,form='unformatted',status='unknown')
      do ldeg=0,ldegup
        if(nogravity)then
          call qpcpsvkern(ldeg)
        else
          fac=dble(ldeg)/dble(ldeggr)
          if(fac.le.1.d0)then
            call qpcpsvkerng(ldeg)
            do iz=1,nz
              do istp=1,4
                do i=1,6
                  ypsv(i,istp,iz)=ypsvg(i,istp,iz)
                enddo
              enddo
            enddo
          else if(fac.ge.1.d0+FLTAPER)then
            call qpcpsvkern(ldeg)
          else
            call qpcpsvkerng(ldeg)
            call qpcpsvkern(ldeg)
            ca=dsin(0.5d0*PI*(fac-1.d0)/FLTAPER)**2
            cb=1.d0-ca
            do iz=1,nz
              do istp=1,4
                do i=1,6
                  ypsv(i,istp,iz)=ca*ypsv(i,istp,iz)
     &                           +cb*ypsvg(i,istp,iz)
                enddo
              enddo
            enddo
          endif
        endif
c
        call qpcshkern(ldeg)
c
        do iz=1,nz
          write(30)((ypsv(i,istp,iz),i=1,6),istp=1,4),
     &             ((ysh(i,istp,iz),i=1,2),istp=1,2)
        enddo
      enddo
      close(30)
c============================================================================
c     iz: receiver depth index
c     il: 0/1/2 for saving coefficients of ldeg-2, ldeg-1 and ldeg
c     istp: 1/2/3/4 fundamental source inflation/strike-slip/dip-slip/clvd
c
      open(30,file=specgrntmp,form='unformatted',status='old')
      do iz=1,nz
        do ldeg=0,ldegup
          dll=dble(ldeg)*dble(ldeg+1)
c
c         each record includes 28 (i.e., 24 psv + 4 sh) times 8 byes data
c         and 2*4 = 8 byes alianment, together 29*8 = 232
c
          offset=232*(ldeg*nz+iz-1)
          call fseek(30,offset,0,i)
          read(30)((ypsv(i,istp,iz),i=1,6),istp=1,4),
     &            ((ysh(i,istp,iz),i=1,2),istp=1,2)
c
          rrr=rearth-depr(iz)
c
c         definitions:
c         (ul0,vl0,wl0)=spherical harmonic vectorial spectrum of displacement
c         (el0,fl0,gl0)=spherical harmonic vectorial spectrum of stress
c         yl0=spherical harmonic spectrum of vertical gravity (downwards positive!)
c         pl0=spherical harmonic spectrum of geopotential
c
          if(lyr(iz).gt.1)then
            a51=-4.d0*PI*BIGG*rholw(lyr(iz)-1)
          else
            a51=0.d0
          endif
          a55=dble(ldeg+1)/rrr
          a56=-1.d0
c
c         1. Explosion (M11=M22=M33=1)
c
          cs1=disk(ldeg)/ksis
          cs2=-4.d0*mues*disk(ldeg)/(ksis*rrs)
          cs4=2.d0*mues*disk(ldeg)/(ksis*rrs)
          ul0(ldeg,1)=cs1*ypsv(1,1,iz)+cs2*ypsv(1,2,iz)+cs4*ypsv(1,4,iz)
          vl0(ldeg,1)=cs1*ypsv(3,1,iz)+cs2*ypsv(3,2,iz)+cs4*ypsv(3,4,iz)
          wl0(ldeg,1)=0.d0
          el0(ldeg,1)=cs1*ypsv(2,1,iz)+cs2*ypsv(2,2,iz)+cs4*ypsv(2,4,iz)
          fl0(ldeg,1)=cs1*ypsv(4,1,iz)+cs2*ypsv(4,2,iz)+cs4*ypsv(4,4,iz)
          gl0(ldeg,1)=0.d0
          ql0(ldeg,1)=cs1*(a51*ypsv(1,1,iz)+a55*ypsv(5,1,iz)
     &                    +a56*ypsv(6,1,iz))
     &               +cs2*(a51*ypsv(1,2,iz)+a55*ypsv(5,2,iz)
     &                    +a56*ypsv(6,2,iz))
     &               +cs4*(a51*ypsv(1,4,iz)+a55*ypsv(5,4,iz)
     &                    +a56*ypsv(6,4,iz))
          pl0(ldeg,1)=cs1*ypsv(5,1,iz)+cs2*ypsv(5,2,iz)+cs4*ypsv(5,4,iz)
c
c         2. Strike-slip (M12=M21=1)
c
          if(ldeg.lt.2)then
            ul0(ldeg,2)=0.d0
            vl0(ldeg,2)=0.d0
            wl0(ldeg,2)=0.d0
            el0(ldeg,2)=0.d0
            fl0(ldeg,2)=0.d0
            gl0(ldeg,2)=0.d0
            ql0(ldeg,2)=0.d0
            pl0(ldeg,2)=0.d0
          else
            ct2=disk(ldeg)/(dll*rrs)
            cs4=-ct2
            ul0(ldeg,2)=cs4*ypsv(1,4,iz)
            vl0(ldeg,2)=cs4*ypsv(3,4,iz)
            wl0(ldeg,2)=ct2*ysh(1,2,iz)
            el0(ldeg,2)=cs4*ypsv(2,4,iz)
            fl0(ldeg,2)=cs4*ypsv(4,4,iz)
            gl0(ldeg,2)=ct2*ysh(2,2,iz)
            ql0(ldeg,2)=cs4*(a51*ypsv(1,4,iz)+a55*ypsv(5,4,iz)
     &                      +a56*ypsv(6,4,iz))
            pl0(ldeg,2)=cs4*ypsv(5,4,iz)
          endif
c
c         3. Dip-slip (M13=M31=1)
c
          if(ldeg.lt.1)then
            ul0(ldeg,3)=0.d0
            vl0(ldeg,3)=0.d0
            wl0(ldeg,3)=0.d0
            el0(ldeg,3)=0.d0
            fl0(ldeg,3)=0.d0
            gl0(ldeg,3)=0.d0
            ql0(ldeg,3)=0.d0
            pl0(ldeg,3)=0.d0
          else
            ct1=disk(ldeg)/(dll*mues)
            cs3=ct1
            ul0(ldeg,3)=cs3*ypsv(1,3,iz)
            vl0(ldeg,3)=cs3*ypsv(3,3,iz)
            wl0(ldeg,3)=ct1*ysh(1,1,iz)
            el0(ldeg,3)=cs3*ypsv(2,3,iz)
            fl0(ldeg,3)=cs3*ypsv(4,3,iz)
            gl0(ldeg,3)=ct1*ysh(2,1,iz)
            ql0(ldeg,3)=cs3*(a51*ypsv(1,3,iz)+a55*ypsv(5,3,iz)
     &                      +a56*ypsv(6,3,iz))
            pl0(ldeg,3)=cs3*ypsv(5,3,iz)
          endif
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c
          cs1=disk(ldeg)/ksis
          cs2=disk(ldeg)*(3.d0-4.d0*mues/ksis)/rrs
          cs4=-0.5d0*cs2
          ul0(ldeg,4)=cs1*ypsv(1,1,iz)+cs2*ypsv(1,2,iz)+cs4*ypsv(1,4,iz)
          vl0(ldeg,4)=cs1*ypsv(3,1,iz)+cs2*ypsv(3,2,iz)+cs4*ypsv(3,4,iz)
          wl0(ldeg,4)=0.d0
          el0(ldeg,4)=cs1*ypsv(2,1,iz)+cs2*ypsv(2,2,iz)+cs4*ypsv(2,4,iz)
          fl0(ldeg,4)=cs1*ypsv(4,1,iz)+cs2*ypsv(4,2,iz)+cs4*ypsv(4,4,iz)
          gl0(ldeg,4)=0.d0
          ql0(ldeg,4)=cs1*(a51*ypsv(1,1,iz)+a55*ypsv(5,1,iz)
     &                    +a56*ypsv(6,1,iz))
     &               +cs2*(a51*ypsv(1,2,iz)+a55*ypsv(5,2,iz)
     &                    +a56*ypsv(6,2,iz))
     &               +cs4*(a51*ypsv(1,4,iz)+a55*ypsv(5,4,iz)
     &                    +a56*ypsv(6,4,iz))
          pl0(ldeg,4)=cs1*ypsv(5,1,iz)+cs2*ypsv(5,2,iz)+cs4*ypsv(5,4,iz)
        enddo
c
c       neglect gravity potential of degree 0 and 1
c
        do istp=1,4
          do ldeg=0,1
            pl0(ldeg,istp)=0.d0
          enddo
        enddo
c
        write(21)((ul0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! Ypsv1
        write(22)((vl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! Ypsv3
        write(23)((wl0(ldeg,istp),ldeg=0,ldegup),istp=2,3)               ! Ysh1 (or Y7)
        write(24)((el0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! Ypsv2
        write(25)((fl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! Ypsv4
        write(26)((gl0(ldeg,istp),ldeg=0,ldegup),istp=2,3)               ! Ysh2 (or Y8)
        write(27)((pl0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! Ypsv5
        write(28)((ql0(ldeg,istp),ldeg=0,ldegup),istp=1,4)               ! dYpsv5/dr
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
c
      close(30)
      call system('del '//specgrntmp)
      return
      end
