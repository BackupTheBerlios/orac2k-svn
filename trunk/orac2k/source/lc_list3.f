      subroutine lc_list3(ncx,ncy,ncz,nind,indxi,indxj,indxk,rcut,co,xp
     &     ,yp,zp,natp,nstart,nend,node,nprocs,ncube,nnlpp,npp_loc
     &     ,npp_o,worka,kprint,n_slt)

************************************************************************
*  Written by Marc Souaille, CECAM 96. Compute linked cell list.       *
*  Cell indexing is done in lc_index.                                  *
************************************************************************

*======================= DECLARATIONS ==================================
      implicit none
      include 'parst.h' 
      include 'lc_list.h' 
  
*----------------------- ARGUMENTS -------------------------------------

      integer ncx,ncy,ncz,nstart,nend,npp_loc,npp_o,nprocs,node,ncube
      integer nind,worka(*)
      integer indxi(*),indxj(*),indxk(*)
      integer natp

      integer nnlpp(*)
      integer   kprint,n_slt

      real*8 xp(*),yp(*),zp(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      character*80 errmsg 
#if defined _CRAY_ | defined T3E
      integer headp(*)
      INTEGER ind_a,M_get_length,indxyz
      POINTER (ip_headp,headp)
#else
      integer headp(CELLMAX)
#endif

      integer chainp(tgroup),iret
      integer cellpi(tgroup),cellpj(tgroup),cellpk(tgroup)
      integer nx,ny,nz
      integer i,j,k,l,m,n,nv,nvtot,numcell
      integer iv,jv,kv,nmin
      integer npww,nppp,nppw,map,count

      real*8 dx,dy,dz,co(3,3),rcut
      real*8 x1,y1,z1,x2,y2,z2,xx,yy,zz 
      real*8 sqcut,d

c--   local stuff is dumped in a scratch common 

#if defined _CRAY_ | defined T3E
      COMMON /rag1/ ip_headp,chainp,cellpi,cellpj,cellpk
#else
      COMMON /rag1/ headp,chainp,cellpi,cellpj,cellpk
#endif

 
      INCLUDE 'pbc.h'
   
*=======================================================================
*                      Initialization
*=======================================================================

#if defined _CRAY_ | defined T3E

*=======================================================================
*---- Allocate memory for headp                                     ----
*=======================================================================

      indxyz=ncx*ncy*ncz
      ind_a=M_get_length(indxyz,4)
      CALL M_memory(ip_headp,ind_a)
#endif

      CALL zero0(worka,natp)
      iret=0
      errmsg=' '
      sqcut = rcut**2
      dx=2.d0/ncx
      dy=2.d0/ncy
      dz=2.d0/ncz

      npww=0
      nppp=0
      nppw=0

      do 10 n=1,ncx*ncy*ncz
         headp(n)=0
10    continue

*=======================================================================
*     Compute chain list for system
*=======================================================================

      do 30 n=1,natp
         x1=xp(n)/dx
         y1=yp(n)/dy
         z1=zp(n)/dz
         nx=int(x1)+(sign(1.d0,x1-int(x1))-1.)/2
         ny=int(y1)+(sign(1.d0,y1-int(y1))-1.)/2
         nz=int(z1)+(sign(1.d0,z1-int(z1))-1.)/2
         nx=mod(mod(nx,ncx)+ncx,ncx)
         ny=mod(mod(ny,ncy)+ncy,ncy)
         nz=mod(mod(nz,ncz)+ncz,ncz)
         cellpi(n)=nx
         cellpj(n)=ny
         cellpk(n)=nz
         numcell=nz+ncz*(ny+ncy*nx)+1
         chainp(n)=headp(numcell)
         headp(numcell)=n
30    continue

*=======================================================================
*     Compute neighbor list nnlpp 
*=======================================================================
      
      nppp=0
      nvtot=1
      do 60 n=nstart,nend
         x1=xp(n)
         y1=yp(n)
         z1=zp(n)
         nv=0         
         i=cellpi(n)
         j=cellpj(n)
         k=cellpk(n)
         do 70 m=1,nind
            iv=indxi(m)
            jv=indxj(m)
            kv=indxk(m)
            nx=mod(mod(i+iv,ncx)+ncx,ncx)
            ny=mod(mod(j+jv,ncy)+ncy,ncy)
            nz=mod(mod(k+kv,ncz)+ncz,ncz)
            numcell=nz+ncz*(ny+ncy*nx)+1
            l=headp(numcell)
            nmin=0
            if(m.eq.1) nmin=n
 
            do while(l.gt.nmin)
               if((l.gt. nend) .AND. (l.le.n_slt)) then
                  x2=x1-xp(l)
                  y2=y1-yp(l)
                  z2=z1-zp(l)
                  x2=x2-2.0*pbc(x2)
                  y2=y2-2.0*pbc(y2)
                  z2=z2-2.0*pbc(z2)
                  xx=co(1,1)*x2+co(1,2)*y2+co(1,3)*z2
                  yy=           co(2,2)*y2+co(2,3)*z2
                  zz=                      co(3,3)*z2
                  d=xx**2+yy**2+zz**2
                  if(d.lt.sqcut) then
                     nppp=nppp+1
                     nv=nv+1
                     nnlpp(nvtot+nv)=l
                  endif
               endif
               l=chainp(l)
            end do
 70      continue         

         if(nvtot .gt. npp_loc) then
            iret=1
            WRITE(errmsg,'('' Neighbor list dim. are'',
     &           '' insufficient. NCOUNT '',i8,'' PHYS_DIM '',
     &           i8,'' PE ='',i4)') nvtot,npp_loc,node
         endif
         nnlpp(nvtot)=nv
         nvtot=nvtot+nv+1
60    continue

#ifdef PARALLEL
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
#endif

#if defined _CRAY_ | defined T3E
      CALL M_free(ip_headp)
#endif
      npp_o=nvtot
      
      map=0
      DO i=nstart,nend
         m=nnlpp(map+1)
         map=map+m+1
         worka(i)=m
      END DO
#ifdef PARALLEL
      IF(nprocs .GT. 1) CALL P_merge_i(nvtot)
#endif
      WRITE(kprint,10000) nvtot
      
*================= END OF EXECUTABLE STATEMENTS ========================

10000 FORMAT(/'Neighbor Lists Dimensions *neighbor(',i8,')* '/)

      RETURN
      END
