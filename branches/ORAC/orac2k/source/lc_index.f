c--------------------------------------------------------------------------
c--------------------------------------------------------------------------
      subroutine lc_index(indmax,ncx,ncy,ncz,nind,indxi,indxj,indxk
     &     ,ctoff,co)

      implicit none
c----------------------------
c declarations des variables
c----------------------------
      real*8 dist_ijk

      real*8 ctoff,co(3,3)

      integer ncx,ncy,ncz,indmax
      integer nind
      integer  indxi(*),indxj(*),indxk(*)

      real*8 dx,dy,dz,sqcut,rmin

      integer i,j,k,n
      integer istart,jstart,kstart
      integer imax,jmax,kmax
      integer nxmax,nymax,nzmax
      integer warnx,warny,warnz

c----------------------------
c initialisations
c----------------------------
      sqcut=ctoff**2
      dx=2.d0/ncx
      dy=2.d0/ncy
      dz=2.d0/ncz
      imax=0
      jmax=0
      kmax=0

c----------------------------
c calcul de la demi-sphere
c des indices
c----------------------------
      nind=1
      indxi(1)=0
      indxj(1)=0
      indxk(1)=0

      istart=1-ncx
      do 10 i=istart,ncx-1
         jstart=1-ncy
         do 20 j=jstart,ncy-1
            kstart=1-ncz
               do 30 k=kstart,ncz-1
                  rmin=dist_ijk(i,j,k,dx,dy,dz,co)

                  if(rmin.lt.sqcut) then
                     nind=nind+1
                     if(nind.gt.INDMAX) then
                        print*,"INDMAX est trop petit"
                        stop
                     endif
                     
                     indxi(nind)=i
                     indxj(nind)=j
                     indxk(nind)=k

                     if(imax.lt.abs(i)) imax=abs(i)
                     if(jmax.lt.abs(j)) jmax=abs(j)
                     if(kmax.lt.abs(k)) kmax=abs(k)
                     if(i.eq.0 .and. j.eq.0 .and. k.eq.0)
     &               nind=nind-1
                  endif
30             continue
20          continue
10       continue
5     continue

c----------------------------
c pour eviter que, lors du
c comptage des paires d'atomes
c une cellule n'apparaisse
c deux fois, il faut ajouter
c le test suivant :
c si ncx est pair,
c il faut imax <= ncx/2
c si ncx est impair,
c il faut imax <= (ncx+1)/2
c (idem pour jmax et kmax)
c----------------------------

      nxmax=(ncx+1)/2
      nymax=(ncy+1)/2
      nzmax=(ncz+1)/2
      warnx=0
      warny=0
      warnz=0
      if(imax.ge.nxmax) warnx=1
      if(jmax.ge.nymax) warny=1
      if(kmax.ge.nzmax) warnz=1
      if(warnx.eq.1 .or. warny.eq.1 .or. warnz.eq.1) then
         print*,"des cellules risquent d'etre comptees deux fois"
         print*,"diminuez le cutoff ou :"
         if(warnx.eq.1) print*,"augmentez ncx"
         if(warny.eq.1) print*,"augmentez ncy"
         if(warnz.eq.1) print*,"augmentez ncz"
         stop
      endif
      return
      end

