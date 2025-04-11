      subroutine qpcpsvkerng(ldeg)
      use qpcalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c
      integer*4 ldeg
c
      integer*4 i,iz,istp,lyup,lylw
c
      do iz=1,nz
        do istp=1,4
          do i=1,6
            ypsvg(i,istp,iz)=0.d0
          enddo
        enddo
      enddo
c
      lyup=lyuppsv(ldeg)
      lylw=lylwpsv(ldeg)
c
      if(ldeg.eq.0)then
        call qpcspropg0(lyup,lylw)
      else
        call qpcspropg(ldeg,lyup,lylw)
      endif
c
      return
      end