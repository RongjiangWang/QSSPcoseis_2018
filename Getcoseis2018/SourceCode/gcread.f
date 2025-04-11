      subroutine gcread(unit,nrg1,nrg2,nzrg)
      use gcalloc
      implicit none
c
      integer*4 unit,nrg1,nrg2,nzrg
c
      integer*4 izrg,irg
c
      do irg=nrg1,nrg2
        read(unit)(ur(izrg,irg,1),izrg=1,nzrg)
        read(unit)(ut(izrg,irg,1),izrg=1,nzrg)
        read(unit)(err(izrg,irg,1),izrg=1,nzrg)
        read(unit)(ert(izrg,irg,1),izrg=1,nzrg)
        read(unit)(etr(izrg,irg,1),izrg=1,nzrg)
        read(unit)(ett(izrg,irg,1),izrg=1,nzrg)
        read(unit)(epp(izrg,irg,1),izrg=1,nzrg)
        read(unit)(gr(izrg,irg,1),izrg=1,nzrg)
        read(unit)(gt(izrg,irg,1),izrg=1,nzrg)
        read(unit)(po(izrg,irg,1),izrg=1,nzrg)
c
        read(unit)(ur(izrg,irg,2),izrg=1,nzrg)
        read(unit)(ut(izrg,irg,2),izrg=1,nzrg)
        read(unit)(up(izrg,irg,2),izrg=1,nzrg)
        read(unit)(err(izrg,irg,2),izrg=1,nzrg)
        read(unit)(ert(izrg,irg,2),izrg=1,nzrg)
        read(unit)(erp(izrg,irg,2),izrg=1,nzrg)
        read(unit)(etr(izrg,irg,2),izrg=1,nzrg)
        read(unit)(ett(izrg,irg,2),izrg=1,nzrg)
        read(unit)(etp(izrg,irg,2),izrg=1,nzrg)
        read(unit)(epr(izrg,irg,2),izrg=1,nzrg)
        read(unit)(ept(izrg,irg,2),izrg=1,nzrg)
        read(unit)(epp(izrg,irg,2),izrg=1,nzrg)
        read(unit)(gr(izrg,irg,2),izrg=1,nzrg)
        read(unit)(gt(izrg,irg,2),izrg=1,nzrg)
        read(unit)(gp(izrg,irg,2),izrg=1,nzrg)
        read(unit)(po(izrg,irg,2),izrg=1,nzrg)
c
        read(unit)(ur(izrg,irg,3),izrg=1,nzrg)
        read(unit)(ut(izrg,irg,3),izrg=1,nzrg)
        read(unit)(up(izrg,irg,3),izrg=1,nzrg)
        read(unit)(err(izrg,irg,3),izrg=1,nzrg)
        read(unit)(ert(izrg,irg,3),izrg=1,nzrg)
        read(unit)(erp(izrg,irg,3),izrg=1,nzrg)
        read(unit)(etr(izrg,irg,3),izrg=1,nzrg)
        read(unit)(ett(izrg,irg,3),izrg=1,nzrg)
        read(unit)(etp(izrg,irg,3),izrg=1,nzrg)
        read(unit)(epr(izrg,irg,3),izrg=1,nzrg)
        read(unit)(ept(izrg,irg,3),izrg=1,nzrg)
        read(unit)(epp(izrg,irg,3),izrg=1,nzrg)
        read(unit)(gr(izrg,irg,3),izrg=1,nzrg)
        read(unit)(gt(izrg,irg,3),izrg=1,nzrg)
        read(unit)(gp(izrg,irg,3),izrg=1,nzrg)
        read(unit)(po(izrg,irg,3),izrg=1,nzrg)
c
        read(unit)(ur(izrg,irg,4),izrg=1,nzrg)
        read(unit)(ut(izrg,irg,4),izrg=1,nzrg)
        read(unit)(err(izrg,irg,4),izrg=1,nzrg)
        read(unit)(ert(izrg,irg,4),izrg=1,nzrg)
        read(unit)(etr(izrg,irg,4),izrg=1,nzrg)
        read(unit)(ett(izrg,irg,4),izrg=1,nzrg)
        read(unit)(epp(izrg,irg,4),izrg=1,nzrg)
        read(unit)(gr(izrg,irg,4),izrg=1,nzrg)
        read(unit)(gt(izrg,irg,4),izrg=1,nzrg)
        read(unit)(po(izrg,irg,4),izrg=1,nzrg)
      enddo
      return
      end
