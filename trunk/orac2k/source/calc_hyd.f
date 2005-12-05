
      SUBROUTINE calc_hyd(kout,ninst,nstep,fstep,nres,nores,nbtype,pnbd1
     &                   ,coeff_hyd,min_hyd,max_hyd,nnlpp,nstart,nend
     &                   ,grppt,co,xpa,ypa,zpa,nato_slt
     &                   ,nmol_slv,res_time,kprot_rest)

*=======================================================================
      implicit none
 

      INTEGER kout,ninst,nstep,min_hyd,max_hyd,nores,nstart,nend
      INTEGER nmol_slv,nato_slt,nres(*),nbtype(*),grppt(2,*)
      INTEGER kprot_rest
      INTEGER nnlpp(*)
      REAL*8 fstep,xpa(*),ypa(*),zpa(*),co(3,3),coeff_hyd,pnbd1(*)
      LOGICAL res_time

c---- LOCAL VARIABLES

      INTEGER max_res,max_occ
      PARAMETER(max_res=1100,max_occ=100)
      REAL*8 nhyd(max_res),norma,dist,dist_LJ,lj_wt
      REAL*8 xv,yv,zv,xs,ys,zs,xg,yg,zg,xc,yc,zc,xmap0,ymap0,zmap0
      INTEGER i,i1,j,j1,k,m,first,sumt,idx_res,prot_res,idx_nnl,old_res
      INTEGER idx_occ,occ(max_occ)
      LOGICAL OCCUPIED
      DATA first/1/

      INCLUDE 'pbc.h'
*=======================================================================





      IF ( first .EQ. 1) THEN
         first = first + 1
         lj_wt = pnbd1(nbtype(nato_slt + 1))

         sumt = 0

         do i=1,max_res
            nhyd(i) = 0.0d0
         end do
c-----   effective number of residues in the solute
         prot_res = nores - nmol_slv

         if(max_hyd .eq. 0) max_hyd = prot_res

         if ( max_res .lt. max_hyd ) THEN
              WRITE(*,*) 
              WRITE(*,*) 'WARNING'
              WRITE(*,*) 'calc_hyd: max_res insufficient'
              WRITE(*,*) 
              STOP
         end if

      END IF

*=======================================================================

      idx_nnl = 0
      idx_occ=0
      do k=1,max_occ
         occ(k) = 0
      end do
c      old_res = nres(grppt(1,nstart))
      old_res = min_hyd

      do i=nstart,nend

         m = nnlpp(idx_nnl+1)

         idx_res = nres(grppt(1,i))

         if(idx_res .lt. min_hyd) go to 10
         if(idx_res .gt. max_hyd) go to 20

         if(idx_res .ne. old_res) then
            if(res_time) then
               write(kprot_rest,200) old_res,idx_occ
               write(kprot_rest,201) (occ(k),k=1,idx_occ)
c               write(kprot_rest) old_res,idx_occ
c               write(kprot_rest) (occ(k),k=1,idx_occ)
            end if
            idx_occ=0
            do k=1,max_occ
               occ(k) = 0
            end do
            old_res = idx_res
         end if
           
         if(m .eq. 0) go to 10

         do j=1,m
            j1=nnlpp(idx_nnl+1+j)

            do k=1,idx_occ
               if(occ(k).eq.(j1-nend)) GO TO 30
            end do

            xv = xpa(grppt(1,j1))
            yv = ypa(grppt(1,j1))
            zv = zpa(grppt(1,j1))

            do i1=grppt(1,i),grppt(2,i)
               xs = xpa(i1)
               ys = ypa(i1)
               zs = zpa(i1)
               dist_LJ = lj_wt + pnbd1(nbtype(i1))
               dist_LJ = dist_lj*coeff_hyd

               xg = xv - xs
               yg = yv - ys
               zg = zv - zs

               xmap0 = 2.0d0*PBC(xg)
               ymap0 = 2.0d0*PBC(yg)
               zmap0 = 2.0d0*PBC(zg)

               xg = xg - xmap0
               yg = yg - ymap0
               zg = zg - zmap0

               xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
               yc=           co(2,2)*yg+co(2,3)*zg
               zc=                      co(3,3)*zg

               dist = xc**2+yc**2+zc**2
               dist = DSQRT(dist)

               IF( dist. LE. dist_LJ ) THEN
                  nhyd(idx_res) = nhyd(idx_res) + 1.0d0
                  idx_occ=idx_occ+1
                  occ(idx_occ)=j1-nend
                  GO TO 30
               END IF

            end do

 30         continue
         end do

 10      continue
         idx_nnl = idx_nnl + 1 + m
      end do

 20   sumt = sumt + 1

c---- write last residue
      if(res_time) then
         write(kprot_rest,200) old_res,idx_occ
         write(kprot_rest,201) (occ(k),k=1,idx_occ)
c         write(kprot_rest) old_res,idx_occ
c         write(kprot_rest) (occ(k),k=1,idx_occ)
      end if

*=======================================================================
            
      IF( MOD(nstep,ninst) .EQ. 0) THEN
          norma = 1.0d0/DFLOAT(sumt)
          WRITE(kout,'('' Tstep = '',f12.2)') fstep
          WRITE(kout,'('' Residue   Hydration number'')')
          DO j=min_hyd,max_hyd
             WRITE(kout,100) j, nhyd(j)*norma
          END DO
          WRITE(kout,*)
      END IF

*=======================================================================

      RETURN
100   FORMAT(2x,i4,2x,f12.6)
200   FORMAT(i4,i4)
201   FORMAT(100(i5))
      END
