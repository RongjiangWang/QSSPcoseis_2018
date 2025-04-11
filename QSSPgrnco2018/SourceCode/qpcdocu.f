      subroutine qpcdocu(ierr)
      use qpcalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,k,ig,flen,iz
      real*8 vp,vs,rho
      character*24 fdate
c
      open(20,file=infofile,status='unknown')
      write(20,'(a)')'#   Space domain Green function database'
      write(20,'(a)')'#   using input file: '//inputfile
      write(20,'(a)')'#   calculated on '//fdate()
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   data_type size: real*8'
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   no of receiver distances'
      write(20,'(i8)')nr
      write(20,'(a)')'#   list of receiver distances [km]'
      j=0
      do i=1,nr/8
        write(20,'(8f10.3)')(dism(k)/KM2M,k=j+1,j+8)
        j=j+8
      enddo
      if(j.lt.nr)then
        write(20,'(8f10.3)')(dism(k)/KM2M,k=j+1,nr)
      endif
      write(20,'(a)')'#   no of receiver depths'
      write(20,'(i8)')nz
      write(20,'(a)')'#   receiver depth profile:'
      write(20,'(a)')'#   depth profile:'
      write(20,'(a)')'#   depth[km]   vp[km/s]   vs[km/s]   rho[g/cm^3]'
      do iz=1,nz
        if(depr(iz).ge.0.d0)then
          vp=dsqrt((kapup(lyr(iz))+mueup(lyr(iz))*4.d0/3.d0)
     &          /rhoup(lyr(iz)))
          vs=dsqrt(mueup(lyr(iz))/rhoup(lyr(iz)))
          rho=rhoup(lyr(iz))
        else
          vp=0.d0
          vs=0.d0
          rho=0.d0
        endif
        write(20,'(4f12.4)')depr(iz)/KM2M,vp/KM2M,vs/KM2M,rho/KM2M
      enddo
      write(20,'(a)')'#   no of source_depths / files,'
      write(20,'(i10)')ngrn
      write(20,'(a)')'#   depth[km] vp[km/s] vs[km/s] rho[g/cm^3]'
     &             //' Green_function_file'
      do ig=1,ngrn
        do flen=80,1,-1
          if(spatgrnfile(ig)(flen:flen).ne.' ')goto 100
        enddo
100     continue
        vp=dsqrt((kapup(lygrn(ig))+mueup(lygrn(ig))*4.d0/3.d0)
     &          /rhoup(lygrn(ig)))
        vs=dsqrt(mueup(lygrn(ig))/rhoup(lygrn(ig)))
        write(20,'(4f12.4,a)')grndep(ig)/KM2M,vp/KM2M,vs/KM2M,
     &                        rhoup(lygrn(ig))/KM2M,
     &                     '  '//''''//spatgrnfile(ig)(1:flen)//''''
      enddo
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   each file includes all fundamental'
      write(20,'(a)')'#   Green functions for all fundamental source'
      write(20,'(a)')'#   at a given depth'
      write(20,'(a)')'#   data need to be read in the following loops:'
      write(20,'(a)')'#   {i_distance = 1, ..., nr'
      write(20,'(a)')'#    {i_mechanism = expl,str-slip,dip-slip,clvd'
      write(20,'(a)')'#      {i_component = ur,ut,up*,err,ert,erp*,'
      write(20,'(a)')'#                     etr,ett,etp*,epr*,ept*,epp,'
      write(20,'(a)')'#                     gr,gt,gp*,gd'
      write(20,'(a)')'#       {list of nz data (depth profile)}}}}'
      write(20,'(a)')'#   note: each data profile'
      write(20,'(a)')'#         needs (4)+number_of_samples*8+(4)'
      write(20,'(a)')'#         bytes, where (4) bytes are for the'
      write(20,'(a)')'#         reserved places before and after each'
      write(20,'(a)')'#         record'
      write(20,'(a)')'#         *)no ue component for explosion'
      write(20,'(a)')'#           and clvd sources'
      write(20,'(a)')'#================================================'
      close(20)
      return
      end