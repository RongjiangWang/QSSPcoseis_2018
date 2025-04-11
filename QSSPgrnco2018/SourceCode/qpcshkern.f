      subroutine qpcshkern(ldeg)
      use qpcalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c
      integer*4 ldeg
c
      integer*4 i,iz,ly,istp,lyup,lylw
c
      do iz=1,nz
        do istp=1,2
          do i=1,2
            ysh(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      lyup=lyupt(ldeg)
      lylw=lylwt(ldeg)
c
      if(ldeg.le.1)then
        return
      else if(ldeg.gt.1)then
        do ly=lyup,lylw
          call qpctmat(ldeg,ly,lylw)
        enddo
        call qpctprop(ldeg,lyup,lylw)
      endif
c
      return
      end

