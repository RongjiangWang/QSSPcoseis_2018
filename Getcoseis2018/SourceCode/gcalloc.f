      module gcalloc
c
      real*8 km2m,rearth,pi
      parameter(km2m=1.0d+03,rearth=6.371d+06,pi=3.14159265358979d+00)
c===================================================================
c     allocatable variables
c===================================================================
      real*8, allocatable:: rrg(:),zrg(:),rhorg(:),vprg(:),vsrg(:)
      real*8, allocatable:: zg(:),rhog(:),vpg(:),vsg(:)
      real*8, allocatable:: ur(:,:,:),ut(:,:,:),up(:,:,:)
      real*8, allocatable:: gr(:,:,:),gt(:,:,:),gp(:,:,:),po(:,:,:)
      real*8, allocatable:: err(:,:,:),ert(:,:,:),erp(:,:,:)
      real*8, allocatable:: etr(:,:,:),ett(:,:,:),etp(:,:,:)
      real*8, allocatable:: epr(:,:,:),ept(:,:,:),epp(:,:,:)
      real*8, allocatable:: dis(:,:),azi(:,:),bazi(:,:)
      logical*2, allocatable:: next(:)
c
      character*80, allocatable:: grnfile(:)
c
      end module