        subroutine gpfft3d(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nsize,
     &             iproc,nproc,z,zt,isign,zw,ncache)
C	CALCULATES THE DISCRETE FOURIERTRANSFORM F(I1,I2,I3)=
C	S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
C       in parallel using MPI and BLAS library calls. 

C	OUTPUT:
C          Z: output array
C		real(F(i1,i2,i3))=Z(1,i1,i2,i3)
C		imag(F(i1,i2,i3))=Z(2,i1,i2,i3)
C	INPUT:
C          Z: input array
C		Z(1,i1,i2,i3)=real(R(i1,i2,i3))
C		Z(2,i1,i2,i3)=imag(R(i1,i2,i3))
C          nproc: number of processors used as returned by MPI_COMM_SIZE
C          iproc: [0:nproc-1] number of processor as returned by MPI_COMM_RANK
C	   n1,n2,n3: logical dimension of the transform. As transform lengths 
C		     most products of the prime factors 2,3,5 are allowed.
C                    The detailed table with allowed transform lengths can 
C                    be found in subroutine CTRIG
C                    If two ni's are equal it is recommended to place them 
C                    behind each other.
C          md1,md2,md3: Dimension of workarrays for MPI_ALLTOALL.
C                       They must be greater than the logical dimensions n1,n2,n3 
C                       and a multiple of the number of processors
C	   nd1,nd2,nd3: Dimension of intermediate single processor workareas. ndi must always 
C                       be greater or equal than ni. It is recomended to chose 
C                       ndi=ni if ni is odd and ndi=ni+1 if ni is even to obtain 
C                       optimal execution speed. 
C          k1=max(nd1,md1),k2=max(nd2,md2): 
C                       Physical dimensions of the input and output array Z,
C                       The input/output array Z is distributed along the z 
C                       direction among the different processors, i.e its 
C                       dimension on a single processor is 
C                       real*8 Z(2,k1,k2,md3/nproc), or complex*8 Z(k1,k2,md3/nproc)
C          k3=max(nd3,md3): needed for nsize
C          nsize=2*k1*k2*k3/nproc: total amount of memory that needs to be allocated 
C                       for both Z and ZT. Note that frequently nsize > 2*k1*k2*md3/nproc
C
C       WORKARRAYS:
C       ZT(nsize),ZW(ncache)
C
C       PERFORMANCE CONSIDERATIONS:
C       The maximum number of processors that can be used in the transformation 
C       steps along the y,z and x axis are n3,n1 and n2. Therefore in general 
C       it is reasonable to choose: nproc < min(n1,n2,n3).
C
C       It is very important to find the optimal 
C       value of NCACHE. NCACHE determines the size of the work array ZW, that
C       has to fit into cache. It has therefore to be chosen to equal roughly 
C	half the size of the physical cache in units of real*8 numbers.
C       The optimal value of ncache can easily be determined by numerical 
C       experimentation. A too large value of ncache leads to a dramatic 
C       and sudden decrease of performance, a too small value to a to a 
C       slow and less dramatic decrease of performance. If NCACHE is set 
C       to a value so small, that not even a single one dimensional transform 
C       can be done in the workarray zw, the program stops with an error message.
C
C       RSTRICTIONS on USAGE
C 	Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
C 	modified by Stefan Goedecker, Stuttgart, Germany, November 5, 1999
C       Distribution  and commercial use is prohibited without the explicit 
C       permission of the author.

        implicit real*8 (a-h,o-z)
#ifdef PARALLEL
        include 'mpif.h'
#endif
	dimension z(nsize),zt(nsize),zw(ncache)

c check input
	if (nd1.lt.n1) stop 'ERROR:nd1'
	if (nd2.lt.n2) stop 'ERROR:nd2'
	if (nd3.lt.n3) stop 'ERROR:nd3'
	if (md1.lt.n1) stop 'ERROR:md1'
	if (md2.lt.n2) stop 'ERROR:md2'
	if (md3.lt.n3) stop 'ERROR:md3'
	if (mod(md1,nproc).ne.0) stop 'ERROR:md1'
	if (mod(md2,nproc).ne.0) stop 'ERROR:md2'
	if (mod(md3,nproc).ne.0) stop 'ERROR:md3'
	if (2*max(md1,nd1)*max(md2,nd2)*max(md3,nd3)/nproc.gt.nsize)
     &      stop 'ERROR:nsize'

c input: i1,i2,i3
        call pfft(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc,
     &             z,zt,isign,zw,ncache)
c ioput: I2,i3,i1
        call pfft(n2,n3,n1,nd2,nd3,nd1,md2,md3,md1,iproc,nproc,
     &             z,zt,isign,zw,ncache)
c ioput: I3,i1,I2
        call pfft(n3,n1,n2,nd3,nd1,nd2,md3,md1,md2,iproc,nproc,
     &             z,zt,isign,zw,ncache)
c output: I1,I2,I3

	return
	end


	subroutine pfft(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,iproc,nproc,
     &             z,zt,isign,zw,ncache)
	implicit real*8 (a-h,o-z)
#ifdef PARALLEL
        include 'mpif.h'
#endif
	integer after,before,now
	dimension trig(2,1024),zw(2,ncache/4,2),
     &	z(2,nd1,nd2,md3/nproc),zt(2,nd2,nd1,md3/nproc),
     &  after(20),now(20),before(20)
	if (max(n1,n2).gt.1024) stop '1024'


	do 2500,j3=1,md3/nproc
	if (iproc*(md3/nproc)+j3.le.n3) then

        lot=max(1,ncache/(4*n2))
c	print*,'lot',lot
	if (2*n2*lot*2.gt.ncache) stop 'enlarge ncache :2'
	call ctrig(n2,trig,after,before,now,isign,ic)

c input: i1,i2,j3,(jp3)
      if (ic.eq.1) then
	i=ic
	nfft=n1
        call fftrot(nd1,nfft,nd2,nd1,nd2,
     &       z(1,1,1,j3),zt(1,1,1,j3),
     &	     trig,after(i),now(i),before(i),isign)

      else

        do 2000,j=1,n1,lot
        ma=j
        mb=min(j+(lot-1),n1)
	nfft=mb-ma+1

	i=1
	inzee=2
	call fftstp(nd1,nfft,nd2,lot,n2,
     &       z(1,j,1,j3),zw(1,1,3-inzee),
     &	     trig,after(i),now(i),before(i),isign)
	inzee=1

	do 2093,i=2,ic-1
	call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee),
     1	trig,after(i),now(i),before(i),isign)
2093	inzee=3-inzee

	i=ic
	call fftrot(lot,nfft,n2,nd1,nd2,
     &       zw(1,1,inzee),zt(1,1,j,j3),
     &	     trig,after(i),now(i),before(i),isign)
2000	continue
      endif
      endif
2500	continue
c output: I2,i1,j3,(jp3)

c input: I2,i1,j3,(jp3)
        call switch(n1,n2,nd1,nd2,md1,md2,md3,nproc,zt,z)
c output: I2,j1,j3,jp1,(jp3)


c input: I2,j1,j3,jp1,(jp3)
	if (nproc.gt.1) then
        call MPI_ALLTOALL(z,2*md2*(md1/nproc)*(md3/nproc), 
     &                    MPI_double_precision,
     &                   zt,2*md2*(md1/nproc)*(md3/nproc),
     &                    MPI_double_precision,MPI_COMM_WORLD,ierr)
	else
	call MYCOPY(2*md2*md1*md3,z,zt)
	endif
c output: I2,j1,j3,jp3,(jp1)

c input: I2,j1,j3,jp3,(jp1)
	call twitch(n2,n3,nd2,nd3,md1,md2,md3,nproc,zt,z)
c output: I2,i3,j1,(jp1)

	return
	end

	subroutine mycopy(n,a,b)
	implicit real*8 (a-h,o-z)
	dimension a(n),b(n)
	do i=1,n-3,4
	b(i+0)=a(i+0)
	b(i+1)=a(i+1)
	b(i+2)=a(i+2)
	b(i+3)=a(i+3)
	enddo
	do i=i,n
	b(i)=a(i)
	enddo
	return
	end


	subroutine switch(n1,n2,nd1,nd2,md1,md2,md3,nproc,zt,z)
	implicit real*8 (a-h,o-z)
	dimension zt(2,nd2,nd1,md3/nproc),
     &           z(2,md2,md1/nproc,md3/nproc,nproc)

	if (mod(n2,2).eq.0.d0) then

	i1=0
	do 200,jp1=1,nproc
	do 200,j1=1,md1/nproc
	i1=i1+1
	if (i1.le.n1) then
        do 210,j3=1,md3/nproc
	do 211,i2=1,n2-1,2
	z(1,i2,j1,j3,jp1)=zt(1,i2,i1,j3)
	z(2,i2,j1,j3,jp1)=zt(2,i2,i1,j3)
	z(1,i2+1,j1,j3,jp1)=zt(1,i2+1,i1,j3)
	z(2,i2+1,j1,j3,jp1)=zt(2,i2+1,i1,j3)
211	continue
210	continue
	endif
200	continue

	else

	i1=0
	do 100,jp1=1,nproc
	do 100,j1=1,md1/nproc
	i1=i1+1
	if (i1.le.n1) then
        do 110,j3=1,md3/nproc
	do 111,i2=1,n2-1,2
	z(1,i2,j1,j3,jp1)=zt(1,i2,i1,j3)
	z(2,i2,j1,j3,jp1)=zt(2,i2,i1,j3)
	z(1,i2+1,j1,j3,jp1)=zt(1,i2+1,i1,j3)
	z(2,i2+1,j1,j3,jp1)=zt(2,i2+1,i1,j3)
111	continue
	z(1,i2,j1,j3,jp1)=zt(1,i2,i1,j3)
	z(2,i2,j1,j3,jp1)=zt(2,i2,i1,j3)
110	continue
	endif
100	continue

	endif

	return
	end


	subroutine twitch(n2,n3,nd2,nd3,md1,md2,md3,nproc,zt,z)
	implicit real*8 (a-h,o-z)
	dimension zt(2,md2,md1/nproc,md3/nproc,nproc),
     &           z(2,nd2,nd3,md1/nproc)

	if (mod(n2,2).eq.0.d0) then

	i3=0
	do 200,jp3=1,nproc
	do 200,j3=1,md3/nproc
	i3=i3+1
	if (i3.le.n3) then
	do 210,j1=1,md1/nproc
	do 211,i2=1,n2-1,2
	z(1,i2,i3,j1)=zt(1,i2,j1,j3,jp3)
	z(2,i2,i3,j1)=zt(2,i2,j1,j3,jp3)
	z(1,i2+1,i3,j1)=zt(1,i2+1,j1,j3,jp3)
	z(2,i2+1,i3,j1)=zt(2,i2+1,j1,j3,jp3)
211	continue
210	continue
	endif
200	continue

	else

	i3=0
	do 100,jp3=1,nproc
	do 100,j3=1,md3/nproc
	i3=i3+1
	if (i3.le.n3) then
	do 110,j1=1,md1/nproc
	do 111,i2=1,n2-1,2
	z(1,i2,i3,j1)=zt(1,i2,j1,j3,jp3)
	z(2,i2,i3,j1)=zt(2,i2,j1,j3,jp3)
	z(1,i2+1,i3,j1)=zt(1,i2+1,j1,j3,jp3)
	z(2,i2+1,i3,j1)=zt(2,i2+1,j1,j3,jp3)
111	continue
	z(1,i2,i3,j1)=zt(1,i2,j1,j3,jp3)
	z(2,i2,i3,j1)=zt(2,i2,j1,j3,jp3)
110	continue
	endif
100	continue

	endif

	return
	end
        subroutine pout(n1,n2,n3,md1,md2,md3,iproc,nproc,z,zin)
	implicit real*8 (a-h,o-z)
	dimension z(2,md1,md2,md3/nproc),zin(2,n1,n2,md3/nproc)
	write(6,*) '----------------------------------'
	scl=1.d0/(n1*n2*n3)
	do 9764,i3=1,md3/nproc
	if (i3+iproc*md3/nproc.le.n3) then
	do 9763,i2=1,n2
	do 9763,i1=1,n1
22	format(3(i4),4(e12.5))
        ttr=abs(z(1,i1,i2,i3)*scl-zin(1,i1,i2,i3))
	tti=abs(z(2,i1,i2,i3)*scl-zin(2,i1,i2,i3))
	if (ttr+tti.gt.1.d-10) then
	write(6,22) i1,i2,i3+iproc*md3/nproc,
     &              scl*z(1,i1,i2,i3),scl*z(2,i1,i2,i3),
     &              zin(1,i1,i2,i3),zin(2,i1,i2,i3)
	endif
9763	continue
	endif
9764	continue
	write(6,*) '++++++++++++++++++++++++++++++++++'
	return
	end



C FFT PART -----------------------------------------------------------------

	subroutine ctrig(n,trig,after,before,now,isign,ic)
C       Copyright by Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
C       modified by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
C       modified by Stefan Goedecker, Stuttgart, Germany, October 6, 1995
C       Commercial use is prohibited without the explicit permission of the author.

	implicit real*8 (a-h,o-z)
	integer after,before
	dimension idata(7,82),now(7),after(7),before(7),
     1  trig(2,1024)

	data idata /
     &      3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,
     &      5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1,
     &      8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,
     &     12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1,
     &     16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,
     &     20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,
     &     25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,
     &     30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1,
     &     36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,
     &     45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,
     &     54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,
     &     64,   4, 4, 4, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1,
     &     75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,
     &     81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,
     &     96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,
     &    108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1,
     &    125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,
     &    135,   5, 3, 3, 3, 1, 1,     144,   4, 4, 3, 3, 1, 1,
     &    150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,
     &    162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,
     &    192,   4, 4, 4, 3, 1, 1,     200,   8, 5, 5, 1, 1, 1,
     &    216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1,
     &    240,   5, 4, 4, 3, 1, 1,     243,   3, 3, 3, 3, 3, 1,
     &    256,   4, 4, 4, 4, 1, 1,     270,   6, 5, 3, 3, 1, 1,
     &    288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,
     &    320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1,
     &    360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,
     &    384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,
     &    405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,
     &    450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1,
     &    486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,
     &    512,   8, 4, 4, 4, 1, 1,     540,   5, 4, 3, 3, 3, 1,
     &    576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,
     &    625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1,
     &    648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,
     &    720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,
     &    750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,
     &    800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1,
     &    864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,
     &    960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,
     &   1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1  /
	do 111,i=1,82
	if (n.eq.idata(1,i)) then
	ic=0
	do 11,j=1,6
	itt=idata(1+j,i)
	if (itt.gt.1) then
	ic=ic+1
	now(j)=idata(1+j,i)
	else
	goto 1000
	endif
11	continue
	goto 1000
	endif
111	continue
	print*,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
37	format(15(i5))
	write(6,37) (idata(1,i),i=1,82)
	stop
1000	continue

	after(1)=1
	before(ic)=1
	do 22,i=2,ic
	after(i)=after(i-1)*now(i-1)
22	before(ic-i+1)=before(ic-i+2)*now(ic-i+2)

12	format(6(i3))
C	write(6,12) (after(i),i=1,ic)
C	write(6,12) (now(i),i=1,ic)
C	write(6,12) (before(i),i=1,ic)

	twopi=8.d0*datan(1.d0)
	angle=isign*twopi/n
	trig(1,1)=1.d0
	trig(2,1)=0.d0
	do 40,i=1,n-1
	trig(1,i+1)=cos(i*angle)
	trig(2,i+1)=sin(i*angle)
40	continue
check
C	if (mod(n,2).eq.0) then
C	trig(2,n/2+1)=-100000.d0
C	endif

	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccc


 	subroutine fftstp(mm,nfft,m,nn,n,zin,zout,
     1	                 trig,after,now,before,isign)
C    Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
C    modified by Stefan Goedecker, Stuttgart, Germany, October 15, 1995
C    Commercial use is prohibited without the explicit permission of the author.
c xeon FFTSTP

	implicit real*8 (a-h,o-z)
	integer after,before,atn,atb
	dimension trig(2,1024),zin(2,mm,m),zout(2,nn,n)
	atn=after*now
	atb=after*before

c	 sqrt(.5d0)
	rt2i=0.7071067811865475d0
	if (now.eq.2) then
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 2001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nout1=nout1+atn
	nout2=nout1+after
	do 2001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	zout(1,j,nout1)= r2 + r1
	zout(2,j,nout1)= s2 + s1
	zout(1,j,nout2)= r1 - r2
	zout(2,j,nout2)= s1 - s2
2001	continue
	do 2000,ia=2,after
	ias=ia-1
	if (2*ias.eq.after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(2,j,nin2)
			s2=zin(1,j,nin2)
			zout(1,j,nout1)= r1 - r2
			zout(2,j,nout1)= s2 + s1
			zout(1,j,nout2)= r2 + r1
			zout(2,j,nout2)= s1 - s2
2010			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(2,j,nin2)
			s2=zin(1,j,nin2)
			zout(1,j,nout1)= r2 + r1
			zout(2,j,nout1)= s1 - s2
			zout(1,j,nout2)= r1 - r2
			zout(2,j,nout2)= s2 + s1
2020			continue
		endif
	else if (4*ias.eq.after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2030,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2030,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r - s)*rt2i
			s2=(r + s)*rt2i
			zout(1,j,nout1)= r2 + r1
			zout(2,j,nout1)= s2 + s1
			zout(1,j,nout2)= r1 - r2
			zout(2,j,nout2)= s1 - s2
2030			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2040,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2040,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			zout(1,j,nout1)= r2 + r1
			zout(2,j,nout1)= s2 + s1
			zout(1,j,nout2)= r1 - r2
			zout(2,j,nout2)= s1 - s2
2040			continue
		endif
	else if (4*ias.eq.3*after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2050,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2050,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(r - s)*rt2i
			zout(1,j,nout1)= r1 - r2
			zout(2,j,nout1)= s2 + s1
			zout(1,j,nout2)= r2 + r1
			zout(2,j,nout2)= s1 - s2
2050			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2060,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2060,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(s - r)*rt2i
			s2=(r + s)*rt2i
			zout(1,j,nout1)= r2 + r1
			zout(2,j,nout1)= s1 - s2
			zout(1,j,nout2)= r1 - r2
			zout(2,j,nout2)= s2 + s1
2060			continue
		endif
	else
		itrig=ias*before+1
		cr2=trig(1,itrig)
		ci2=trig(2,itrig)
		nin1=ia-after
		nout1=ia-atn
		do 2090,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nout1=nout1+atn
		nout2=nout1+after
		do 2090,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=r*cr2 - s*ci2
		s2=r*ci2 + s*cr2
		zout(1,j,nout1)= r2 + r1
		zout(2,j,nout1)= s2 + s1
		zout(1,j,nout2)= r1 - r2
		zout(2,j,nout2)= s1 - s2
2090		continue
	endif
2000	continue
	else if (now.eq.4) then
	if (isign.eq.1) then 
		ia=1
		nin1=ia-after
		nout1=ia-atn
		do 4001,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		do 4001,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(1,j,nin2)
		s2=zin(2,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r4=zin(1,j,nin4)
		s4=zin(2,j,nin4)
		r=r1 + r3
		s=r2 + r4
		zout(1,j,nout1) = r + s
		zout(1,j,nout3) = r - s
		r=r1 - r3
		s=s2 - s4
		zout(1,j,nout2) = r - s 
		zout(1,j,nout4) = r + s
		r=s1 + s3
		s=s2 + s4
		zout(2,j,nout1) = r + s 
		zout(2,j,nout3) = r - s
		r=s1 - s3
		s=r2 - r4
		zout(2,j,nout2) = r + s 
		zout(2,j,nout4) = r - s
4001		continue
		do 4000,ia=2,after
		ias=ia-1
		if (2*ias.eq.after) then
			nin1=ia-after
			nout1=ia-atn
			do 4010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r-s)*rt2i
			s2=(r+s)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(r + s)*rt2i
			s4=(r - s)*rt2i
			r=r1 - r3
			s=r2 - r4
			zout(1,j,nout1) = r + s
			zout(1,j,nout3) = r - s
			r=r1 + r3
			s=s2 - s4
			zout(1,j,nout2) = r - s 
			zout(1,j,nout4) = r + s
			r=s1 + s3
			s=s2 + s4
			zout(2,j,nout1) = r + s 
			zout(2,j,nout3) = r - s
			r=s1 - s3
			s=r2 + r4
			zout(2,j,nout2) = r + s 
			zout(2,j,nout4) = r - s
4010			continue
		else
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 4020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=r1 + r3
			s=r2 + r4
			zout(1,j,nout1) = r + s
			zout(1,j,nout3) = r - s
			r=r1 - r3
			s=s2 - s4
			zout(1,j,nout2) = r - s 
			zout(1,j,nout4) = r + s
			r=s1 + s3
			s=s2 + s4
			zout(2,j,nout1) = r + s 
			zout(2,j,nout3) = r - s
			r=s1 - s3
			s=r2 - r4
			zout(2,j,nout2) = r + s 
			zout(2,j,nout4) = r - s
4020			continue
		endif
4000		continue
	else
		ia=1
		nin1=ia-after
		nout1=ia-atn
		do 4101,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		do 4101,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(1,j,nin2)
		s2=zin(2,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r4=zin(1,j,nin4)
		s4=zin(2,j,nin4)
		r=r1 + r3
		s=r2 + r4
		zout(1,j,nout1) = r + s
		zout(1,j,nout3) = r - s
		r=r1 - r3
		s=s2 - s4
		zout(1,j,nout2) = r + s
		zout(1,j,nout4) = r - s
		r=s1 + s3
		s=s2 + s4
		zout(2,j,nout1) = r + s
		zout(2,j,nout3) = r - s
		r=s1 - s3
		s=r2 - r4
		zout(2,j,nout2) = r - s
		zout(2,j,nout4) = r + s
4101		continue
		do 4100,ia=2,after
		ias=ia-1
		if (2*ias.eq.after) then
			nin1=ia-after
			nout1=ia-atn
			do 4110,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4110,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(s - r)*rt2i
			s4=(r + s)*rt2i
			r=r1 + r3
			s=r2 + r4
			zout(1,j,nout1) = r + s
			zout(1,j,nout3) = r - s
			r=r1 - r3
			s=s2 + s4
			zout(1,j,nout2) = r + s
			zout(1,j,nout4) = r - s
			r=s1 - s3
			s=s2 - s4
			zout(2,j,nout1) = r + s
			zout(2,j,nout3) = r - s
			r=s1 + s3
			s=r2 - r4
			zout(2,j,nout2) = r - s
			zout(2,j,nout4) = r + s
4110			continue
		else
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 4120,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4120,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=r1 + r3
			s=r2 + r4
			zout(1,j,nout1) = r + s
			zout(1,j,nout3) = r - s
			r=r1 - r3
			s=s2 - s4
			zout(1,j,nout2) = r + s
			zout(1,j,nout4) = r - s
			r=s1 + s3
			s=s2 + s4
			zout(2,j,nout1) = r + s
			zout(2,j,nout3) = r - s
			r=s1 - s3
			s=r2 - r4
			zout(2,j,nout2) = r - s
			zout(2,j,nout4) = r + s
4120			continue
		endif
4100		continue
	endif
	else if (now.eq.8) then
	if (isign.eq.-1) then 
		ia=1
			nin1=ia-after
			nout1=ia-atn
			do 8120,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8120,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(1,j,nin2)
			s2=zin(2,j,nin2)
			r3=zin(1,j,nin3)
			s3=zin(2,j,nin3)
			r4=zin(1,j,nin4)
			s4=zin(2,j,nin4)
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r6=zin(1,j,nin6)
			s6=zin(2,j,nin6)
			r7=zin(1,j,nin7)
			s7=zin(2,j,nin7)
			r8=zin(1,j,nin8)
			s8=zin(2,j,nin8)
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,j,nout1) = ap + bp
			zout(2,j,nout1) = cp + dp
			zout(1,j,nout5) = ap - bp
			zout(2,j,nout5) = cp - dp
			zout(1,j,nout3) = am + dm
			zout(2,j,nout3) = cm - bm
			zout(1,j,nout7) = am - dm
			zout(2,j,nout7) = cm + bm
			r=r1 - r5
			s=s3 - s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r3 - r7
			bp=r + s
			bm=r - s
			r=s4 - s8
			s=r2 - r6
			cp=r + s
			cm=r - s
			r=s2 - s6
			s=r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( dm - cp)*rt2i
			cp= ( cm + dp)*rt2i
			dp = ( cm - dp)*rt2i
			zout(1,j,nout2) = ap + r
			zout(2,j,nout2) = bm + s
			zout(1,j,nout6) = ap - r
			zout(2,j,nout6) = bm - s
			zout(1,j,nout4) = am + cp
			zout(2,j,nout4) = bp + dp
			zout(1,j,nout8) = am - cp
			zout(2,j,nout8) = bp - dp
8120			continue
		do 8000,ia=2,after
		ias=ia-1
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			itrig=itrig+itt
			cr5=trig(1,itrig)
			ci5=trig(2,itrig)
			itrig=itrig+itt
			cr6=trig(1,itrig)
			ci6=trig(2,itrig)
			itrig=itrig+itt
			cr7=trig(1,itrig)
			ci7=trig(2,itrig)
			itrig=itrig+itt
			cr8=trig(1,itrig)
			ci8=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 8020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=zin(1,j,nin5)
			s=zin(2,j,nin5)
			r5=r*cr5 - s*ci5
			s5=r*ci5 + s*cr5
			r=zin(1,j,nin6)
			s=zin(2,j,nin6)
			r6=r*cr6 - s*ci6
			s6=r*ci6 + s*cr6
			r=zin(1,j,nin7)
			s=zin(2,j,nin7)
			r7=r*cr7 - s*ci7
			s7=r*ci7 + s*cr7
			r=zin(1,j,nin8)
			s=zin(2,j,nin8)
			r8=r*cr8 - s*ci8
			s8=r*ci8 + s*cr8
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,j,nout1) = ap + bp
			zout(2,j,nout1) = cp + dp
			zout(1,j,nout5) = ap - bp
			zout(2,j,nout5) = cp - dp
			zout(1,j,nout3) = am + dm
			zout(2,j,nout3) = cm - bm
			zout(1,j,nout7) = am - dm
			zout(2,j,nout7) = cm + bm
			r=r1 - r5
			s=s3 - s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r3 - r7
			bp=r + s
			bm=r - s
			r=s4 - s8
			s=r2 - r6
			cp=r + s
			cm=r - s
			r=s2 - s6
			s=r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( dm - cp)*rt2i
			cp= ( cm + dp)*rt2i
			dp = ( cm - dp)*rt2i
			zout(1,j,nout2) = ap + r
			zout(2,j,nout2) = bm + s
			zout(1,j,nout6) = ap - r
			zout(2,j,nout6) = bm - s
			zout(1,j,nout4) = am + cp
			zout(2,j,nout4) = bp + dp
			zout(1,j,nout8) = am - cp
			zout(2,j,nout8) = bp - dp

8020			continue
8000		continue

	else
		ia=1
			nin1=ia-after
			nout1=ia-atn
			do 8121,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8121,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(1,j,nin2)
			s2=zin(2,j,nin2)
			r3=zin(1,j,nin3)
			s3=zin(2,j,nin3)
			r4=zin(1,j,nin4)
			s4=zin(2,j,nin4)
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r6=zin(1,j,nin6)
			s6=zin(2,j,nin6)
			r7=zin(1,j,nin7)
			s7=zin(2,j,nin7)
			r8=zin(1,j,nin8)
			s8=zin(2,j,nin8)
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,j,nout1) = ap + bp
			zout(2,j,nout1) = cp + dp
			zout(1,j,nout5) = ap - bp
			zout(2,j,nout5) = cp - dp
			zout(1,j,nout3) = am - dm
			zout(2,j,nout3) = cm + bm
			zout(1,j,nout7) = am + dm
			zout(2,j,nout7) = cm - bm
			r= r1 - r5
			s=-s3 + s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r7 - r3
			bp=r + s
			bm=r - s
			r=-s4 + s8
			s= r2 - r6
			cp=r + s
			cm=r - s
			r=-s2 + s6
			s= r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( cp - dm)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( dp - cm)*rt2i
			zout(1,j,nout2) = ap + r
			zout(2,j,nout2) = bm + s
			zout(1,j,nout6) = ap - r
			zout(2,j,nout6) = bm - s
			zout(1,j,nout4) = am + cp
			zout(2,j,nout4) = bp + dp
			zout(1,j,nout8) = am - cp
			zout(2,j,nout8) = bp - dp
8121			continue

		do 8001,ia=2,after
		ias=ia-1
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			itrig=itrig+itt
			cr5=trig(1,itrig)
			ci5=trig(2,itrig)
			itrig=itrig+itt
			cr6=trig(1,itrig)
			ci6=trig(2,itrig)
			itrig=itrig+itt
			cr7=trig(1,itrig)
			ci7=trig(2,itrig)
			itrig=itrig+itt
			cr8=trig(1,itrig)
			ci8=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 8021,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8021,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=zin(1,j,nin5)
			s=zin(2,j,nin5)
			r5=r*cr5 - s*ci5
			s5=r*ci5 + s*cr5
			r=zin(1,j,nin6)
			s=zin(2,j,nin6)
			r6=r*cr6 - s*ci6
			s6=r*ci6 + s*cr6
			r=zin(1,j,nin7)
			s=zin(2,j,nin7)
			r7=r*cr7 - s*ci7
			s7=r*ci7 + s*cr7
			r=zin(1,j,nin8)
			s=zin(2,j,nin8)
			r8=r*cr8 - s*ci8
			s8=r*ci8 + s*cr8
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,j,nout1) = ap + bp
			zout(2,j,nout1) = cp + dp
			zout(1,j,nout5) = ap - bp
			zout(2,j,nout5) = cp - dp
			zout(1,j,nout3) = am - dm
			zout(2,j,nout3) = cm + bm
			zout(1,j,nout7) = am + dm
			zout(2,j,nout7) = cm - bm
			r= r1 - r5
			s=-s3 + s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r7 - r3
			bp=r + s
			bm=r - s
			r=-s4 + s8
			s= r2 - r6
			cp=r + s
			cm=r - s
			r=-s2 + s6
			s= r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( cp - dm)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( dp - cm)*rt2i
			zout(1,j,nout2) = ap + r
			zout(2,j,nout2) = bm + s
			zout(1,j,nout6) = ap - r
			zout(2,j,nout6) = bm - s
			zout(1,j,nout4) = am + cp
			zout(2,j,nout4) = bp + dp
			zout(1,j,nout8) = am - cp
			zout(2,j,nout8) = bp - dp
8021			continue
8001		continue

	endif
	else if (now.eq.3) then 
c	 .5d0*sqrt(3.d0)
	bb=isign*0.8660254037844387d0
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 3001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	do 3001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	r3=zin(1,j,nin3)
	s3=zin(2,j,nin3)
	r=r2 + r3
	s=s2 + s3
	zout(1,j,nout1) = r + r1
	zout(2,j,nout1) = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r2=bb*(r2-r3)
	s2=bb*(s2-s3)
	zout(1,j,nout2) = r1 - s2 
	zout(2,j,nout2) = s1 + r2
	zout(1,j,nout3) = r1 + s2 
	zout(2,j,nout3) = s1 - r2
3001	continue
	do 3000,ia=2,after
	ias=ia-1
	if (4*ias.eq.3*after) then
	if (isign.eq.1) then
		nin1=ia-after
		nout1=ia-atn
		do 3010,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3010,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(2,j,nin2)
		s2=zin(1,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r=r3 + r2
		s=s2 - s3
		zout(1,j,nout1) = r1 - r
		zout(2,j,nout1) = s + s1
		r1=r1 + .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2-r3)	
		s2=bb*(s2+s3)
		zout(1,j,nout2) = r1 - s2 
		zout(2,j,nout2) = s1 - r2
		zout(1,j,nout3) = r1 + s2 
		zout(2,j,nout3) = s1 + r2
3010		continue
	else
		nin1=ia-after
		nout1=ia-atn
		do 3020,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3020,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(2,j,nin2)
		s2=zin(1,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r=r2 - r3
		s=s2 + s3
		zout(1,j,nout1) = r + r1
		zout(2,j,nout1) = s1 - s
		r1=r1 - .5d0*r
		s1=s1 + .5d0*s	
		r2=bb*(r2+r3)	
		s2=bb*(s2-s3)
		zout(1,j,nout2) = r1 + s2 
		zout(2,j,nout2) = s1 + r2
		zout(1,j,nout3) = r1 - s2 
		zout(2,j,nout3) = s1 - r2
3020		continue
	endif
	else if (8*ias.eq.3*after) then
	if (isign.eq.1) then
		nin1=ia-after
		nout1=ia-atn
		do 3030,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3030,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=(r - s)*rt2i
		s2=(r + s)*rt2i
		r3=zin(2,j,nin3)
		s3=zin(1,j,nin3) 
		r=r2 - r3
		s=s2 + s3
		zout(1,j,nout1) = r + r1
		zout(2,j,nout1) = s + s1
		r1=r1 - .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2+r3)	
		s2=bb*(s2-s3)
		zout(1,j,nout2) = r1 - s2 
		zout(2,j,nout2) = s1 + r2
		zout(1,j,nout3) = r1 + s2 
		zout(2,j,nout3) = s1 - r2
3030		continue
	else
		nin1=ia-after
		nout1=ia-atn
		do 3040,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3040,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=(r + s)*rt2i
		s2=(s - r)*rt2i
		r3=zin(2,j,nin3)
		s3=zin(1,j,nin3)
		r=r2 + r3
		s=s2 - s3
		zout(1,j,nout1) = r + r1
		zout(2,j,nout1) = s + s1
		r1=r1 - .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2-r3)	
		s2=bb*(s2+s3)
		zout(1,j,nout2) = r1 - s2 
		zout(2,j,nout2) = s1 + r2
		zout(1,j,nout3) = r1 + s2 
		zout(2,j,nout3) = s1 - r2
3040		continue
	endif
	else
	itt=ias*before
	itrig=itt+1
	cr2=trig(1,itrig)
	ci2=trig(2,itrig)
	itrig=itrig+itt
	cr3=trig(1,itrig)
	ci3=trig(2,itrig)
	nin1=ia-after
	nout1=ia-atn
	do 3090,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	do 3090,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r=zin(1,j,nin2)
	s=zin(2,j,nin2)
	r2=r*cr2 - s*ci2
	s2=r*ci2 + s*cr2
	r=zin(1,j,nin3)
	s=zin(2,j,nin3)
	r3=r*cr3 - s*ci3
	s3=r*ci3 + s*cr3
	r=r2 + r3
	s=s2 + s3
	zout(1,j,nout1) = r + r1
	zout(2,j,nout1) = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r2=bb*(r2-r3)
	s2=bb*(s2-s3)
	zout(1,j,nout2) = r1 - s2 
	zout(2,j,nout2) = s1 + r2
	zout(1,j,nout3) = r1 + s2 
	zout(2,j,nout3) = s1 - r2
3090	continue
	endif
3000	continue
	else if (now.eq.5) then
c	 cos(2.d0*pi/5.d0)
	cos2=0.3090169943749474d0
c	 cos(4.d0*pi/5.d0)
	cos4=-0.8090169943749474d0
c	sin(2.d0*pi/5.d0)
	sin2=isign*0.9510565162951536d0
c	 sin(4.d0*pi/5.d0)
	sin4=isign*0.5877852522924731d0
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 5001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nin4=nin3+atb
	nin5=nin4+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	nout4=nout3+after
	nout5=nout4+after
	do 5001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	r3=zin(1,j,nin3)
	s3=zin(2,j,nin3)
	r4=zin(1,j,nin4)
	s4=zin(2,j,nin4)
	r5=zin(1,j,nin5)
	s5=zin(2,j,nin5)
	r25 = r2 + r5
	r34 = r3 + r4
	s25 = s2 - s5
	s34 = s3 - s4
	zout(1,j,nout1) = r1 + r25 + r34
	r = r1 + cos2*r25 + cos4*r34
	s = sin2*s25 + sin4*s34
	zout(1,j,nout2) = r - s
	zout(1,j,nout5) = r + s
	r = r1 + cos4*r25 + cos2*r34
	s = sin4*s25 - sin2*s34
	zout(1,j,nout3) = r - s
	zout(1,j,nout4) = r + s
	r25 = r2 - r5
	r34 = r3 - r4
	s25 = s2 + s5
	s34 = s3 + s4
	zout(2,j,nout1) = s1 + s25 + s34
	r = s1 + cos2*s25 + cos4*s34
	s = sin2*r25 + sin4*r34
	zout(2,j,nout2) = r + s
	zout(2,j,nout5) = r - s
	r = s1 + cos4*s25 + cos2*s34
	s = sin4*r25 - sin2*r34
	zout(2,j,nout3) = r + s
	zout(2,j,nout4) = r - s
5001	continue
	do 5000,ia=2,after
	ias=ia-1
	if (8*ias.eq.5*after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 5010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb	
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			do 5010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r - s)*rt2i
			s2=(r + s)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3) 
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(r + s)*rt2i
			s4=(r - s)*rt2i
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r25 = r2 - r5
			r34 = r3 + r4
			s25 = s2 + s5
			s34 = s3 - s4
			zout(1,j,nout1) = r1 + r25 - r34
			r = r1 + cos2*r25 - cos4*r34 
			s = sin2*s25 + sin4*s34
			zout(1,j,nout2) = r - s
			zout(1,j,nout5) = r + s
			r = r1 + cos4*r25 - cos2*r34 
			s = sin4*s25 - sin2*s34
			zout(1,j,nout3) = r - s
			zout(1,j,nout4) = r + s
			r25 = r2 + r5
			r34 = r4 - r3
			s25 = s2 - s5
			s34 = s3 + s4
			zout(2,j,nout1) = s1 + s25 + s34
			r = s1 + cos2*s25 + cos4*s34
			s = sin2*r25 + sin4*r34
			zout(2,j,nout2) = r + s
			zout(2,j,nout5) = r - s
			r = s1 + cos4*s25 + cos2*s34
			s = sin4*r25 - sin2*r34
			zout(2,j,nout3) = r + s
			zout(2,j,nout4) = r - s
5010			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 5020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb	
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			do 5020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(s - r)*rt2i
			s4=(r + s)*rt2i
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r25 = r2 - r5
			r34 = r3 + r4
			s25 = s2 + s5
			s34 = s4 - s3
			zout(1,j,nout1) = r1 + r25 + r34
			r = r1 + cos2*r25 + cos4*r34
			s = sin2*s25 + sin4*s34
			zout(1,j,nout2) = r - s
			zout(1,j,nout5) = r + s
			r = r1 + cos4*r25 + cos2*r34
			s = sin4*s25 - sin2*s34
			zout(1,j,nout3) = r - s
			zout(1,j,nout4) = r + s
			r25 = r2 + r5
			r34 = r3 - r4
			s25 = s2 - s5
			s34 = s3 + s4
			zout(2,j,nout1) = s1 + s25 - s34
			r = s1 + cos2*s25 - cos4*s34
			s = sin2*r25 + sin4*r34
			zout(2,j,nout2) = r + s
			zout(2,j,nout5) = r - s
			r = s1 + cos4*s25 - cos2*s34
			s = sin4*r25 - sin2*r34
			zout(2,j,nout3) = r + s
			zout(2,j,nout4) = r - s
5020			continue
		endif
	else
		ias=ia-1
		itt=ias*before
		itrig=itt+1
		cr2=trig(1,itrig)
		ci2=trig(2,itrig)
		itrig=itrig+itt
		cr3=trig(1,itrig)
		ci3=trig(2,itrig)
		itrig=itrig+itt
		cr4=trig(1,itrig)
		ci4=trig(2,itrig)
		itrig=itrig+itt
		cr5=trig(1,itrig)
		ci5=trig(2,itrig)
		nin1=ia-after
		nout1=ia-atn
		do 5100,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nin5=nin4+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		nout5=nout4+after
		do 5100,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=r*cr2 - s*ci2
		s2=r*ci2 + s*cr2
		r=zin(1,j,nin3)
		s=zin(2,j,nin3)
		r3=r*cr3 - s*ci3
		s3=r*ci3 + s*cr3
		r=zin(1,j,nin4)
		s=zin(2,j,nin4)
		r4=r*cr4 - s*ci4
		s4=r*ci4 + s*cr4
		r=zin(1,j,nin5)
		s=zin(2,j,nin5)
		r5=r*cr5 - s*ci5
		s5=r*ci5 + s*cr5
		r25 = r2 + r5
		r34 = r3 + r4
		s25 = s2 - s5
		s34 = s3 - s4
		zout(1,j,nout1) = r1 + r25 + r34
		r = r1 + cos2*r25 + cos4*r34
		s = sin2*s25 + sin4*s34
		zout(1,j,nout2) = r - s
		zout(1,j,nout5) = r + s
		r = r1 + cos4*r25 + cos2*r34
		s = sin4*s25 - sin2*s34
		zout(1,j,nout3) = r - s
		zout(1,j,nout4) = r + s
		r25 = r2 - r5
		r34 = r3 - r4
		s25 = s2 + s5
		s34 = s3 + s4
		zout(2,j,nout1) = s1 + s25 + s34
		r = s1 + cos2*s25 + cos4*s34
		s = sin2*r25 + sin4*r34
		zout(2,j,nout2) = r + s
		zout(2,j,nout5) = r - s
		r = s1 + cos4*s25 + cos2*s34
		s = sin4*r25 - sin2*r34
		zout(2,j,nout3) = r + s
		zout(2,j,nout4) = r - s
5100		continue
	endif
5000	continue
       else if (now.eq.6) then
c	 .5d0*sqrt(3.d0)
	bb=isign*0.8660254037844387d0

	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 6120,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nin4=nin3+atb
	nin5=nin4+atb
	nin6=nin5+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	nout4=nout3+after
	nout5=nout4+after
	nout6=nout5+after
	do 6120,j=1,nfft
	r2=zin(1,j,nin3)
	s2=zin(2,j,nin3)
	r3=zin(1,j,nin5)
	s3=zin(2,j,nin5)
	r=r2 + r3
	s=s2 + s3
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	ur1 = r + r1
	ui1 = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r=r2-r3
	s=s2-s3
	ur2 = r1 - s*bb
	ui2 = s1 + r*bb
	ur3 = r1 + s*bb
	ui3 = s1 - r*bb

	r2=zin(1,j,nin6)
	s2=zin(2,j,nin6)
	r3=zin(1,j,nin2)
	s3=zin(2,j,nin2)
	r=r2 + r3
	s=s2 + s3
	r1=zin(1,j,nin4)
	s1=zin(2,j,nin4)
	vr1 = r + r1
	vi1 = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r=r2-r3
	s=s2-s3
	vr2 = r1 - s*bb
	vi2 = s1 + r*bb
	vr3 = r1 + s*bb
	vi3 = s1 - r*bb

	zout(1,j,nout1)=ur1+vr1
	zout(2,j,nout1)=ui1+vi1
	zout(1,j,nout5)=ur2+vr2
	zout(2,j,nout5)=ui2+vi2
	zout(1,j,nout3)=ur3+vr3
	zout(2,j,nout3)=ui3+vi3
	zout(1,j,nout4)=ur1-vr1
	zout(2,j,nout4)=ui1-vi1
	zout(1,j,nout2)=ur2-vr2
	zout(2,j,nout2)=ui2-vi2
	zout(1,j,nout6)=ur3-vr3
	zout(2,j,nout6)=ui3-vi3
6120	continue

	else 
	stop 'error fftstp'
	endif

	return
	end



 	subroutine fftrot(mm,nfft,m,nn,n,zin,zout,
     1	                 trig,after,now,before,isign)
C    Copyright by Stefan Goedecker, Cornell, Ithaca, USA, March 25, 1994
C    modified by Stefan Goedecker, Stuttgart, Germany, October 15, 1995
C    Commercial use is prohibited without the explicit permission of the author.
c xeon FFTSTP

	implicit real*8 (a-h,o-z)
	integer after,before,atn,atb
	dimension trig(2,1024),zin(2,mm,m),zout(2,n,nn)
	atn=after*now
	atb=after*before

c	 sqrt(.5d0)
	rt2i=0.7071067811865475d0
	if (now.eq.2) then
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 2001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nout1=nout1+atn
	nout2=nout1+after
	do 2001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	zout(1,nout1,j)= r2 + r1
	zout(2,nout1,j)= s2 + s1
	zout(1,nout2,j)= r1 - r2
	zout(2,nout2,j)= s1 - s2
2001	continue
	do 2000,ia=2,after
	ias=ia-1
	if (2*ias.eq.after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(2,j,nin2)
			s2=zin(1,j,nin2)
			zout(1,nout1,j)= r1 - r2
			zout(2,nout1,j)= s2 + s1
			zout(1,nout2,j)= r2 + r1
			zout(2,nout2,j)= s1 - s2
2010			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(2,j,nin2)
			s2=zin(1,j,nin2)
			zout(1,nout1,j)= r2 + r1
			zout(2,nout1,j)= s1 - s2
			zout(1,nout2,j)= r1 - r2
			zout(2,nout2,j)= s2 + s1
2020			continue
		endif
	else if (4*ias.eq.after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2030,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2030,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r - s)*rt2i
			s2=(r + s)*rt2i
			zout(1,nout1,j)= r2 + r1
			zout(2,nout1,j)= s2 + s1
			zout(1,nout2,j)= r1 - r2
			zout(2,nout2,j)= s1 - s2
2030			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2040,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2040,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			zout(1,nout1,j)= r2 + r1
			zout(2,nout1,j)= s2 + s1
			zout(1,nout2,j)= r1 - r2
			zout(2,nout2,j)= s1 - s2
2040			continue
		endif
	else if (4*ias.eq.3*after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 2050,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2050,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(r - s)*rt2i
			zout(1,nout1,j)= r1 - r2
			zout(2,nout1,j)= s2 + s1
			zout(1,nout2,j)= r2 + r1
			zout(2,nout2,j)= s1 - s2
2050			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 2060,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nout1=nout1+atn
			nout2=nout1+after
			do 2060,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(s - r)*rt2i
			s2=(r + s)*rt2i
			zout(1,nout1,j)= r2 + r1
			zout(2,nout1,j)= s1 - s2
			zout(1,nout2,j)= r1 - r2
			zout(2,nout2,j)= s2 + s1
2060			continue
		endif
	else
		itrig=ias*before+1
		cr2=trig(1,itrig)
		ci2=trig(2,itrig)
		nin1=ia-after
		nout1=ia-atn
		do 2090,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nout1=nout1+atn
		nout2=nout1+after
		do 2090,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=r*cr2 - s*ci2
		s2=r*ci2 + s*cr2
		zout(1,nout1,j)= r2 + r1
		zout(2,nout1,j)= s2 + s1
		zout(1,nout2,j)= r1 - r2
		zout(2,nout2,j)= s1 - s2
2090		continue
	endif
2000	continue
	else if (now.eq.4) then
	if (isign.eq.1) then 
		ia=1
		nin1=ia-after
		nout1=ia-atn
		do 4001,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		do 4001,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(1,j,nin2)
		s2=zin(2,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r4=zin(1,j,nin4)
		s4=zin(2,j,nin4)
		r=r1 + r3
		s=r2 + r4
		zout(1,nout1,j) = r + s
		zout(1,nout3,j) = r - s
		r=r1 - r3
		s=s2 - s4
		zout(1,nout2,j) = r - s 
		zout(1,nout4,j) = r + s
		r=s1 + s3
		s=s2 + s4
		zout(2,nout1,j) = r + s 
		zout(2,nout3,j) = r - s
		r=s1 - s3
		s=r2 - r4
		zout(2,nout2,j) = r + s 
		zout(2,nout4,j) = r - s
4001		continue
		do 4000,ia=2,after
		ias=ia-1
		if (2*ias.eq.after) then
			nin1=ia-after
			nout1=ia-atn
			do 4010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r-s)*rt2i
			s2=(r+s)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(r + s)*rt2i
			s4=(r - s)*rt2i
			r=r1 - r3
			s=r2 - r4
			zout(1,nout1,j) = r + s
			zout(1,nout3,j) = r - s
			r=r1 + r3
			s=s2 - s4
			zout(1,nout2,j) = r - s 
			zout(1,nout4,j) = r + s
			r=s1 + s3
			s=s2 + s4
			zout(2,nout1,j) = r + s 
			zout(2,nout3,j) = r - s
			r=s1 - s3
			s=r2 + r4
			zout(2,nout2,j) = r + s 
			zout(2,nout4,j) = r - s
4010			continue
		else
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 4020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=r1 + r3
			s=r2 + r4
			zout(1,nout1,j) = r + s
			zout(1,nout3,j) = r - s
			r=r1 - r3
			s=s2 - s4
			zout(1,nout2,j) = r - s 
			zout(1,nout4,j) = r + s
			r=s1 + s3
			s=s2 + s4
			zout(2,nout1,j) = r + s 
			zout(2,nout3,j) = r - s
			r=s1 - s3
			s=r2 - r4
			zout(2,nout2,j) = r + s 
			zout(2,nout4,j) = r - s
4020			continue
		endif
4000		continue
	else
		ia=1
		nin1=ia-after
		nout1=ia-atn
		do 4101,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		do 4101,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(1,j,nin2)
		s2=zin(2,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r4=zin(1,j,nin4)
		s4=zin(2,j,nin4)
		r=r1 + r3
		s=r2 + r4
		zout(1,nout1,j) = r + s
		zout(1,nout3,j) = r - s
		r=r1 - r3
		s=s2 - s4
		zout(1,nout2,j) = r + s
		zout(1,nout4,j) = r - s
		r=s1 + s3
		s=s2 + s4
		zout(2,nout1,j) = r + s
		zout(2,nout3,j) = r - s
		r=s1 - s3
		s=r2 - r4
		zout(2,nout2,j) = r - s
		zout(2,nout4,j) = r + s
4101		continue
		do 4100,ia=2,after
		ias=ia-1
		if (2*ias.eq.after) then
			nin1=ia-after
			nout1=ia-atn
			do 4110,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4110,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(s - r)*rt2i
			s4=(r + s)*rt2i
			r=r1 + r3
			s=r2 + r4
			zout(1,nout1,j) = r + s
			zout(1,nout3,j) = r - s
			r=r1 - r3
			s=s2 + s4
			zout(1,nout2,j) = r + s
			zout(1,nout4,j) = r - s
			r=s1 - s3
			s=s2 - s4
			zout(2,nout1,j) = r + s
			zout(2,nout3,j) = r - s
			r=s1 + s3
			s=r2 - r4
			zout(2,nout2,j) = r - s
			zout(2,nout4,j) = r + s
4110			continue
		else
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 4120,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			do 4120,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=r1 + r3
			s=r2 + r4
			zout(1,nout1,j) = r + s
			zout(1,nout3,j) = r - s
			r=r1 - r3
			s=s2 - s4
			zout(1,nout2,j) = r + s
			zout(1,nout4,j) = r - s
			r=s1 + s3
			s=s2 + s4
			zout(2,nout1,j) = r + s
			zout(2,nout3,j) = r - s
			r=s1 - s3
			s=r2 - r4
			zout(2,nout2,j) = r - s
			zout(2,nout4,j) = r + s
4120			continue
		endif
4100		continue
	endif
	else if (now.eq.8) then
	if (isign.eq.-1) then 
		ia=1
			nin1=ia-after
			nout1=ia-atn
			do 8120,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8120,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(1,j,nin2)
			s2=zin(2,j,nin2)
			r3=zin(1,j,nin3)
			s3=zin(2,j,nin3)
			r4=zin(1,j,nin4)
			s4=zin(2,j,nin4)
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r6=zin(1,j,nin6)
			s6=zin(2,j,nin6)
			r7=zin(1,j,nin7)
			s7=zin(2,j,nin7)
			r8=zin(1,j,nin8)
			s8=zin(2,j,nin8)
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,nout1,j) = ap + bp
			zout(2,nout1,j) = cp + dp
			zout(1,nout5,j) = ap - bp
			zout(2,nout5,j) = cp - dp
			zout(1,nout3,j) = am + dm
			zout(2,nout3,j) = cm - bm
			zout(1,nout7,j) = am - dm
			zout(2,nout7,j) = cm + bm
			r=r1 - r5
			s=s3 - s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r3 - r7
			bp=r + s
			bm=r - s
			r=s4 - s8
			s=r2 - r6
			cp=r + s
			cm=r - s
			r=s2 - s6
			s=r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( dm - cp)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( cm - dp)*rt2i
			zout(1,nout2,j) = ap + r
			zout(2,nout2,j) = bm + s
			zout(1,nout6,j) = ap - r
			zout(2,nout6,j) = bm - s
			zout(1,nout4,j) = am + cp
			zout(2,nout4,j) = bp + dp
			zout(1,nout8,j) = am - cp
			zout(2,nout8,j) = bp - dp
8120			continue
		do 8000,ia=2,after
		ias=ia-1
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			itrig=itrig+itt
			cr5=trig(1,itrig)
			ci5=trig(2,itrig)
			itrig=itrig+itt
			cr6=trig(1,itrig)
			ci6=trig(2,itrig)
			itrig=itrig+itt
			cr7=trig(1,itrig)
			ci7=trig(2,itrig)
			itrig=itrig+itt
			cr8=trig(1,itrig)
			ci8=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 8020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=zin(1,j,nin5)
			s=zin(2,j,nin5)
			r5=r*cr5 - s*ci5
			s5=r*ci5 + s*cr5
			r=zin(1,j,nin6)
			s=zin(2,j,nin6)
			r6=r*cr6 - s*ci6
			s6=r*ci6 + s*cr6
			r=zin(1,j,nin7)
			s=zin(2,j,nin7)
			r7=r*cr7 - s*ci7
			s7=r*ci7 + s*cr7
			r=zin(1,j,nin8)
			s=zin(2,j,nin8)
			r8=r*cr8 - s*ci8
			s8=r*ci8 + s*cr8
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,nout1,j) = ap + bp
			zout(2,nout1,j) = cp + dp
			zout(1,nout5,j) = ap - bp
			zout(2,nout5,j) = cp - dp
			zout(1,nout3,j) = am + dm
			zout(2,nout3,j) = cm - bm
			zout(1,nout7,j) = am - dm
			zout(2,nout7,j) = cm + bm
			r=r1 - r5
			s=s3 - s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r3 - r7
			bp=r + s
			bm=r - s
			r=s4 - s8
			s=r2 - r6
			cp=r + s
			cm=r - s
			r=s2 - s6
			s=r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( dm - cp)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( cm - dp)*rt2i
			zout(1,nout2,j) = ap + r
			zout(2,nout2,j) = bm + s
			zout(1,nout6,j) = ap - r
			zout(2,nout6,j) = bm - s
			zout(1,nout4,j) = am + cp
			zout(2,nout4,j) = bp + dp
			zout(1,nout8,j) = am - cp
			zout(2,nout8,j) = bp - dp

8020			continue
8000		continue

	else
		ia=1
			nin1=ia-after
			nout1=ia-atn
			do 8121,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8121,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r2=zin(1,j,nin2)
			s2=zin(2,j,nin2)
			r3=zin(1,j,nin3)
			s3=zin(2,j,nin3)
			r4=zin(1,j,nin4)
			s4=zin(2,j,nin4)
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r6=zin(1,j,nin6)
			s6=zin(2,j,nin6)
			r7=zin(1,j,nin7)
			s7=zin(2,j,nin7)
			r8=zin(1,j,nin8)
			s8=zin(2,j,nin8)
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,nout1,j) = ap + bp
			zout(2,nout1,j) = cp + dp
			zout(1,nout5,j) = ap - bp
			zout(2,nout5,j) = cp - dp
			zout(1,nout3,j) = am - dm
			zout(2,nout3,j) = cm + bm
			zout(1,nout7,j) = am + dm
			zout(2,nout7,j) = cm - bm
			r= r1 - r5
			s=-s3 + s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r7 - r3
			bp=r + s
			bm=r - s
			r=-s4 + s8
			s= r2 - r6
			cp=r + s
			cm=r - s
			r=-s2 + s6
			s= r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( cp - dm)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( dp - cm)*rt2i
			zout(1,nout2,j) = ap + r
			zout(2,nout2,j) = bm + s
			zout(1,nout6,j) = ap - r
			zout(2,nout6,j) = bm - s
			zout(1,nout4,j) = am + cp
			zout(2,nout4,j) = bp + dp
			zout(1,nout8,j) = am - cp
			zout(2,nout8,j) = bp - dp
8121			continue

		do 8001,ia=2,after
		ias=ia-1
			itt=ias*before
			itrig=itt+1
			cr2=trig(1,itrig)
			ci2=trig(2,itrig)
			itrig=itrig+itt
			cr3=trig(1,itrig)
			ci3=trig(2,itrig)
			itrig=itrig+itt
			cr4=trig(1,itrig)
			ci4=trig(2,itrig)
			itrig=itrig+itt
			cr5=trig(1,itrig)
			ci5=trig(2,itrig)
			itrig=itrig+itt
			cr6=trig(1,itrig)
			ci6=trig(2,itrig)
			itrig=itrig+itt
			cr7=trig(1,itrig)
			ci7=trig(2,itrig)
			itrig=itrig+itt
			cr8=trig(1,itrig)
			ci8=trig(2,itrig)
			nin1=ia-after
			nout1=ia-atn
			do 8021,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb
			nin6=nin5+atb
			nin7=nin6+atb
			nin8=nin7+atb
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			nout6=nout5+after
			nout7=nout6+after
			nout8=nout7+after
			do 8021,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=r*cr2 - s*ci2
			s2=r*ci2 + s*cr2
			r=zin(1,j,nin3)
			s=zin(2,j,nin3)
			r3=r*cr3 - s*ci3
			s3=r*ci3 + s*cr3
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=r*cr4 - s*ci4
			s4=r*ci4 + s*cr4
			r=zin(1,j,nin5)
			s=zin(2,j,nin5)
			r5=r*cr5 - s*ci5
			s5=r*ci5 + s*cr5
			r=zin(1,j,nin6)
			s=zin(2,j,nin6)
			r6=r*cr6 - s*ci6
			s6=r*ci6 + s*cr6
			r=zin(1,j,nin7)
			s=zin(2,j,nin7)
			r7=r*cr7 - s*ci7
			s7=r*ci7 + s*cr7
			r=zin(1,j,nin8)
			s=zin(2,j,nin8)
			r8=r*cr8 - s*ci8
			s8=r*ci8 + s*cr8
			r=r1 + r5
			s=r3 + r7
			ap=r + s
			am=r - s
			r=r2 + r6
			s=r4 + r8
			bp=r + s
			bm=r - s
			r=s1 + s5
			s=s3 + s7
			cp=r + s
			cm=r - s
			r=s2 + s6
			s=s4 + s8
			dp=r + s
			dm=r - s
			zout(1,nout1,j) = ap + bp
			zout(2,nout1,j) = cp + dp
			zout(1,nout5,j) = ap - bp
			zout(2,nout5,j) = cp - dp
			zout(1,nout3,j) = am - dm
			zout(2,nout3,j) = cm + bm
			zout(1,nout7,j) = am + dm
			zout(2,nout7,j) = cm - bm
			r= r1 - r5
			s=-s3 + s7
			ap=r + s
			am=r - s
			r=s1 - s5
			s=r7 - r3
			bp=r + s
			bm=r - s
			r=-s4 + s8
			s= r2 - r6
			cp=r + s
			cm=r - s
			r=-s2 + s6
			s= r4 - r8
			dp=r + s
			dm=r - s
			r = ( cp + dm)*rt2i
			s = ( cp - dm)*rt2i
			cp= ( cm + dp)*rt2i
			dp= ( dp - cm)*rt2i
			zout(1,nout2,j) = ap + r
			zout(2,nout2,j) = bm + s
			zout(1,nout6,j) = ap - r
			zout(2,nout6,j) = bm - s
			zout(1,nout4,j) = am + cp
			zout(2,nout4,j) = bp + dp
			zout(1,nout8,j) = am - cp
			zout(2,nout8,j) = bp - dp
8021			continue
8001		continue

	endif
	else if (now.eq.3) then 
c	 .5d0*sqrt(3.d0)
	bb=isign*0.8660254037844387d0
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 3001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	do 3001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	r3=zin(1,j,nin3)
	s3=zin(2,j,nin3)
	r=r2 + r3
	s=s2 + s3
	zout(1,nout1,j) = r + r1
	zout(2,nout1,j) = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r2=bb*(r2-r3)
	s2=bb*(s2-s3)
	zout(1,nout2,j) = r1 - s2 
	zout(2,nout2,j) = s1 + r2
	zout(1,nout3,j) = r1 + s2 
	zout(2,nout3,j) = s1 - r2
3001	continue
	do 3000,ia=2,after
	ias=ia-1
	if (4*ias.eq.3*after) then
	if (isign.eq.1) then
		nin1=ia-after
		nout1=ia-atn
		do 3010,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3010,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(2,j,nin2)
		s2=zin(1,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r=r2 + r3
		s=s2 - s3
		zout(1,nout1,j) = r1 - r
		zout(2,nout1,j) = s + s1
		r1=r1 + .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2-r3)	
		s2=bb*(s2+s3)
		zout(1,nout2,j) = r1 - s2 
		zout(2,nout2,j) = s1 - r2
		zout(1,nout3,j) = r1 + s2 
		zout(2,nout3,j) = s1 + r2
3010		continue
	else
		nin1=ia-after
		nout1=ia-atn
		do 3020,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3020,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r2=zin(2,j,nin2)
		s2=zin(1,j,nin2)
		r3=zin(1,j,nin3)
		s3=zin(2,j,nin3)
		r=r2 - r3
		s=s2 + s3
		zout(1,nout1,j) = r + r1
		zout(2,nout1,j) = s1 - s
		r1=r1 - .5d0*r
		s1=s1 + .5d0*s	
		r2=bb*(r2+r3)	
		s2=bb*(s2-s3)
		zout(1,nout2,j) = r1 + s2 
		zout(2,nout2,j) = s1 + r2
		zout(1,nout3,j) = r1 - s2 
		zout(2,nout3,j) = s1 - r2
3020		continue
	endif
	else if (8*ias.eq.3*after) then
	if (isign.eq.1) then
		nin1=ia-after
		nout1=ia-atn
		do 3030,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3030,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=(r - s)*rt2i
		s2=(r + s)*rt2i
		r3=zin(2,j,nin3)
		s3=zin(1,j,nin3) 
		r=r2 - r3
		s=s2 + s3
		zout(1,nout1,j) = r + r1
		zout(2,nout1,j) = s + s1
		r1=r1 - .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2+r3)	
		s2=bb*(s2-s3)
		zout(1,nout2,j) = r1 - s2 
		zout(2,nout2,j) = s1 + r2
		zout(1,nout3,j) = r1 + s2 
		zout(2,nout3,j) = s1 - r2
3030		continue
	else
		nin1=ia-after
		nout1=ia-atn
		do 3040,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		do 3040,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=(r + s)*rt2i
		s2=(s - r)*rt2i
		r3=zin(2,j,nin3)
		s3=zin(1,j,nin3)
		r=r2 + r3
		s=s2 - s3
		zout(1,nout1,j) = r + r1
		zout(2,nout1,j) = s + s1
		r1=r1 - .5d0*r
		s1=s1 - .5d0*s	
		r2=bb*(r2-r3)	
		s2=bb*(s2+s3)
		zout(1,nout2,j) = r1 - s2 
		zout(2,nout2,j) = s1 + r2
		zout(1,nout3,j) = r1 + s2 
		zout(2,nout3,j) = s1 - r2
3040		continue
	endif
	else
	itt=ias*before
	itrig=itt+1
	cr2=trig(1,itrig)
	ci2=trig(2,itrig)
	itrig=itrig+itt
	cr3=trig(1,itrig)
	ci3=trig(2,itrig)
	nin1=ia-after
	nout1=ia-atn
	do 3090,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	do 3090,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r=zin(1,j,nin2)
	s=zin(2,j,nin2)
	r2=r*cr2 - s*ci2
	s2=r*ci2 + s*cr2
	r=zin(1,j,nin3)
	s=zin(2,j,nin3)
	r3=r*cr3 - s*ci3
	s3=r*ci3 + s*cr3
	r=r2 + r3
	s=s2 + s3
	zout(1,nout1,j) = r + r1
	zout(2,nout1,j) = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r2=bb*(r2-r3)
	s2=bb*(s2-s3)
	zout(1,nout2,j) = r1 - s2 
	zout(2,nout2,j) = s1 + r2
	zout(1,nout3,j) = r1 + s2 
	zout(2,nout3,j) = s1 - r2
3090	continue
	endif
3000	continue
	else if (now.eq.5) then
c	 cos(2.d0*pi/5.d0)
	cos2=0.3090169943749474d0
c	 cos(4.d0*pi/5.d0)
	cos4=-0.8090169943749474d0
c	sin(2.d0*pi/5.d0)
	sin2=isign*0.9510565162951536d0
c	 sin(4.d0*pi/5.d0)
	sin4=isign*0.5877852522924731d0
	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 5001,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nin4=nin3+atb
	nin5=nin4+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	nout4=nout3+after
	nout5=nout4+after
	do 5001,j=1,nfft
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	r2=zin(1,j,nin2)
	s2=zin(2,j,nin2)
	r3=zin(1,j,nin3)
	s3=zin(2,j,nin3)
	r4=zin(1,j,nin4)
	s4=zin(2,j,nin4)
	r5=zin(1,j,nin5)
	s5=zin(2,j,nin5)
	r25 = r2 + r5
	r34 = r3 + r4
	s25 = s2 - s5
	s34 = s3 - s4
	zout(1,nout1,j) = r1 + r25 + r34
	r = r1 + cos2*r25 + cos4*r34
	s = sin2*s25 + sin4*s34
	zout(1,nout2,j) = r - s
	zout(1,nout5,j) = r + s
	r = r1 + cos4*r25 + cos2*r34
	s = sin4*s25 - sin2*s34
	zout(1,nout3,j) = r - s
	zout(1,nout4,j) = r + s
	r25 = r2 - r5
	r34 = r3 - r4
	s25 = s2 + s5
	s34 = s3 + s4
	zout(2,nout1,j) = s1 + s25 + s34
	r = s1 + cos2*s25 + cos4*s34
	s = sin2*r25 + sin4*r34
	zout(2,nout2,j) = r + s
	zout(2,nout5,j) = r - s
	r = s1 + cos4*s25 + cos2*s34
	s = sin4*r25 - sin2*r34
	zout(2,nout3,j) = r + s
	zout(2,nout4,j) = r - s
5001	continue
	do 5000,ia=2,after
	ias=ia-1
	if (8*ias.eq.5*after) then
		if (isign.eq.1) then
			nin1=ia-after
			nout1=ia-atn
			do 5010,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb	
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			do 5010,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r - s)*rt2i
			s2=(r + s)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3) 
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(r + s)*rt2i
			s4=(r - s)*rt2i
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r25 = r2 - r5
			r34 = r3 + r4
			s25 = s2 + s5
			s34 = s3 - s4
			zout(1,nout1,j) = r1 + r25 - r34
			r = r1 + cos2*r25 - cos4*r34
			s = sin2*s25 + sin4*s34
			zout(1,nout2,j) = r - s
			zout(1,nout5,j) = r + s
			r = r1 + cos4*r25 - cos2*r34
			s = sin4*s25 - sin2*s34
			zout(1,nout3,j) = r - s
			zout(1,nout4,j) = r + s
			r25 = r2 + r5
			r34 = r4 - r3
			s25 = s2 - s5
			s34 = s3 + s4
			zout(2,nout1,j) = s1 + s25 + s34
			r = s1 + cos2*s25 + cos4*s34
			s = sin2*r25 + sin4*r34
			zout(2,nout2,j) = r + s
			zout(2,nout5,j) = r - s
			r = s1 + cos4*s25 + cos2*s34
			s = sin4*r25 - sin2*r34
			zout(2,nout3,j) = r + s
			zout(2,nout4,j) = r - s
5010			continue
		else
			nin1=ia-after
			nout1=ia-atn
			do 5020,ib=1,before
			nin1=nin1+after
			nin2=nin1+atb
			nin3=nin2+atb
			nin4=nin3+atb
			nin5=nin4+atb	
			nout1=nout1+atn
			nout2=nout1+after
			nout3=nout2+after
			nout4=nout3+after
			nout5=nout4+after
			do 5020,j=1,nfft
			r1=zin(1,j,nin1)
			s1=zin(2,j,nin1)
			r=zin(1,j,nin2)
			s=zin(2,j,nin2)
			r2=(r + s)*rt2i
			s2=(s - r)*rt2i
			r3=zin(2,j,nin3)
			s3=zin(1,j,nin3)
			r=zin(1,j,nin4)
			s=zin(2,j,nin4)
			r4=(s - r)*rt2i
			s4=(r + s)*rt2i
			r5=zin(1,j,nin5)
			s5=zin(2,j,nin5)
			r25 = r2 - r5
			r34 = r3 + r4
			s25 = s2 + s5
			s34 = s4 - s3
			zout(1,nout1,j) = r1 + r25 + r34
			r = r1 + cos2*r25 + cos4*r34
			s = sin2*s25 + sin4*s34
			zout(1,nout2,j) = r - s
			zout(1,nout5,j) = r + s
			r = r1 + cos4*r25 + cos2*r34
			s = sin4*s25 - sin2*s34
			zout(1,nout3,j) = r - s
			zout(1,nout4,j) = r + s
			r25 = r2 + r5
			r34 = r3 - r4
			s25 = s2 - s5
			s34 = s3 + s4
			zout(2,nout1,j) = s1 + s25 - s34
			r = s1 + cos2*s25 - cos4*s34
			s = sin2*r25 + sin4*r34
			zout(2,nout2,j) = r + s
			zout(2,nout5,j) = r - s
			r = s1 + cos4*s25 - cos2*s34
			s = sin4*r25 - sin2*r34
			zout(2,nout3,j) = r + s
			zout(2,nout4,j) = r - s
5020			continue
		endif
	else
		ias=ia-1
		itt=ias*before
		itrig=itt+1
		cr2=trig(1,itrig)
		ci2=trig(2,itrig)
		itrig=itrig+itt
		cr3=trig(1,itrig)
		ci3=trig(2,itrig)
		itrig=itrig+itt
		cr4=trig(1,itrig)
		ci4=trig(2,itrig)
		itrig=itrig+itt
		cr5=trig(1,itrig)
		ci5=trig(2,itrig)
		nin1=ia-after
		nout1=ia-atn
		do 5100,ib=1,before
		nin1=nin1+after
		nin2=nin1+atb
		nin3=nin2+atb
		nin4=nin3+atb
		nin5=nin4+atb
		nout1=nout1+atn
		nout2=nout1+after
		nout3=nout2+after
		nout4=nout3+after
		nout5=nout4+after
		do 5100,j=1,nfft
		r1=zin(1,j,nin1)
		s1=zin(2,j,nin1)
		r=zin(1,j,nin2)
		s=zin(2,j,nin2)
		r2=r*cr2 - s*ci2
		s2=r*ci2 + s*cr2
		r=zin(1,j,nin3)
		s=zin(2,j,nin3)
		r3=r*cr3 - s*ci3
		s3=r*ci3 + s*cr3
		r=zin(1,j,nin4)
		s=zin(2,j,nin4)
		r4=r*cr4 - s*ci4
		s4=r*ci4 + s*cr4
		r=zin(1,j,nin5)
		s=zin(2,j,nin5)
		r5=r*cr5 - s*ci5
		s5=r*ci5 + s*cr5
		r25 = r2 + r5
		r34 = r3 + r4
		s25 = s2 - s5
		s34 = s3 - s4
		zout(1,nout1,j) = r1 + r25 + r34
		r = r1 + cos2*r25 + cos4*r34
		s = sin2*s25 + sin4*s34
		zout(1,nout2,j) = r - s
		zout(1,nout5,j) = r + s
		r = r1 + cos4*r25 + cos2*r34
		s = sin4*s25 - sin2*s34
		zout(1,nout3,j) = r - s
		zout(1,nout4,j) = r + s
		r25 = r2 - r5
		r34 = r3 - r4
		s25 = s2 + s5
		s34 = s3 + s4
		zout(2,nout1,j) = s1 + s25 + s34
		r = s1 + cos2*s25 + cos4*s34
		s = sin2*r25 + sin4*r34
		zout(2,nout2,j) = r + s
		zout(2,nout5,j) = r - s
		r = s1 + cos4*s25 + cos2*s34
		s = sin4*r25 - sin2*r34
		zout(2,nout3,j) = r + s
		zout(2,nout4,j) = r - s
5100		continue
	endif
5000	continue
       else if (now.eq.6) then
c	 .5d0*sqrt(3.d0)
	bb=isign*0.8660254037844387d0

	ia=1
	nin1=ia-after
	nout1=ia-atn
	do 6120,ib=1,before
	nin1=nin1+after
	nin2=nin1+atb
	nin3=nin2+atb
	nin4=nin3+atb
	nin5=nin4+atb
	nin6=nin5+atb
	nout1=nout1+atn
	nout2=nout1+after
	nout3=nout2+after
	nout4=nout3+after
	nout5=nout4+after
	nout6=nout5+after
	do 6120,j=1,nfft
	r2=zin(1,j,nin3)
	s2=zin(2,j,nin3)
	r3=zin(1,j,nin5)
	s3=zin(2,j,nin5)
	r=r2 + r3
	s=s2 + s3
	r1=zin(1,j,nin1)
	s1=zin(2,j,nin1)
	ur1 = r + r1
	ui1 = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r=r2-r3
	s=s2-s3
	ur2 = r1 - s*bb
	ui2 = s1 + r*bb
	ur3 = r1 + s*bb
	ui3 = s1 - r*bb

	r2=zin(1,j,nin6)
	s2=zin(2,j,nin6)
	r3=zin(1,j,nin2)
	s3=zin(2,j,nin2)
	r=r2 + r3
	s=s2 + s3
	r1=zin(1,j,nin4)
	s1=zin(2,j,nin4)
	vr1 = r + r1
	vi1 = s + s1
	r1=r1 - .5d0*r
	s1=s1 - .5d0*s
	r=r2-r3
	s=s2-s3
	vr2 = r1 - s*bb
	vi2 = s1 + r*bb
	vr3 = r1 + s*bb
	vi3 = s1 - r*bb

	zout(1,nout1,j)=ur1+vr1
	zout(2,nout1,j)=ui1+vi1
	zout(1,nout5,j)=ur2+vr2
	zout(2,nout5,j)=ui2+vi2
	zout(1,nout3,j)=ur3+vr3
	zout(2,nout3,j)=ui3+vi3
	zout(1,nout4,j)=ur1-vr1
	zout(2,nout4,j)=ui1-vi1
	zout(1,nout2,j)=ur2-vr2
	zout(2,nout2,j)=ui2-vi2
	zout(1,nout6,j)=ur3-vr3
	zout(2,nout6,j)=ui3-vi3
6120	continue

       else
	stop 'error fftrot'
       endif

	return
	end
