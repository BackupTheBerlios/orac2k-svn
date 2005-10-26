      SUBROUTINE analyse_voronoi(nstart_ah,nend_ah,nlocal_ah,nstart_uh
     &     ,nend_uh,nlocal_uh,node,nprocs,ncube,fstep,volume_co
     &     ,voronoi_volume,voronoi_access,voronoi_contact
     &     ,voronoi_neighbor,ncontact_slt,contact_slt,voronoi_res
     &     ,kvoronoi,ss_index,ss_point,grppt,mend,protl,nprot,nato,nbun
     &     ,beta,atomp,res1,res2,mres,prsymb,iret,errmsg)

************************************************************************
*   Time-stamp: <99/02/12 12:51:59 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jun 22 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart_ah,nend_ah,nlocal_ah,nstart_uh,nend_uh,nlocal_uh
     &     ,node,nprocs,ncube
      INTEGER iret,kvoronoi,nato,nprot,nbun,ncontact_slt,contact_slt(2,
     &     *),voronoi_res(*)
      INTEGER res1(*),res2(*),ss_index(*),protl(*),atomp(*),mres(2,*)
     &     ,grppt(2,*),mend(*),ss_point(*)
      REAL*8  fstep,volume_co
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)
      CHARACTER*80 errmsg
      LOGICAL voronoi_volume,voronoi_access,voronoi_contact
     &     ,voronoi_neighbor

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'voronoi.h'
      REAL*8  vol_res(m1),vol_slv,vol_slt(npm),area_slt(nores)
     &     ,area_slv(nores),area_frac(nores),area_tot(nores)
     &     ,area_contact(pcontact_slt),area_prot(npm),vol_type(nores)
      INTEGER nvol_type,index(nores),mmap(3),vol_typep(nores)
      COMMON /rag1/ vol_res,vol_slv,vol_slt,area_slt,area_slv
     &     ,area_contact,area_prot,vol_type,index,mmap,vol_typep

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER igo,ign,map,i1,j,count,m,n,i,typei,typeij,p1,p2,j1,ia,ib
     &     ,ic,o
      REAL*8  vol,volume,area2,tot,sum
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      igo=0
      map=0

#ifdef PARALLEL
      IF(nprocs .NE. 1) CALL P_expand_r8(volume_vor,nstart_ah,nend_ah
     &     ,nlocal_ah,node,nprocs)
#endif

*=======================================================================
*--- Compute Voronoi volume per residue --------------------------------
*=======================================================================


      IF(node .EQ. 0) THEN
         IF(voronoi_volume) THEN
            volume=0.0D0
            DO i1=1,nato
               ign=res1(i1)
               volume=volume+volume_vor(i1)
               IF(i1 .EQ. 1) THEN
                  vol=volume_vor(i1)
               ELSE IF(i1 .EQ. nato) THEN
                  map=map+1
                  vol=vol+volume_vor(i1)
                  vol_res(map)=vol
               ELSE 
                  IF(igo .EQ. ign) THEN
                     vol=vol+volume_vor(i1)
                  ELSE IF(igo .NE. ign) THEN
                     map=map+1
                     vol_res(map)=vol
                     vol=volume_vor(i1)
                  END IF
               END IF
               igo=ign
            END DO
            IF(voronoi_res(1) .EQ. 0) THEN
               IF(map .NE. 0) THEN
                  DO i1=1,map,4
                     n=3
                     IF(map-i1 .LT. 3) n=map-i1
                     WRITE(kvoronoi,'(a1,f11.2,1x,a6,4(f10.3,1x,i4)
     &                    )') 'T',fstep,'VolRes',(vol_res(j),j,j=i1,i1+n
     &                    )
                  END DO
               END IF
            ELSE
               IF(map .NE. 0) THEN
                  DO i1=1,voronoi_res(1),4
                     n=3
                     IF(voronoi_res(1)-i1 .LT. 3) n=voronoi_res(1)-i1
                     WRITE(kvoronoi,'(a1,f11.2,1x,a6,4(f10.3,1x,i4)
     &                    )') 'T',fstep,'VolRes',(vol_res(voronoi_res(j
     &                    +1)),voronoi_res(j+1),j=i1,i1+n)
                  END DO
               END IF
            END IF
         END IF


*=======================================================================
*--- Compute Voronoi volume of solvent and solute ----------------------
*=======================================================================

         count=0
         map=0
         vol_slv=0.0D0
         DO n=1,nprot
            m=protl(count+1)
            vol_slt(n)=0.0D0
            DO i=1,m
               i1=protl(count+1+i)
               typei=ss_index(i1)
               IF(typei .EQ. 1) THEN
                  IF(i .EQ. 1) map=map+1
                  vol_slt(n)=vol_slt(n)+volume_vor(i1)
               ELSE
                  vol_slv=vol_slv+volume_vor(i1)
               END IF
            END DO
            count=count+m+1
         END DO
         volume=0.0D0
         DO i1=1,map
            volume=volume+vol_slt(i1)
         END DO
         WRITE(kvoronoi,'(''T'',f11.2,1x,''VolSlt'',f12.3,
     &        '' VolSlv'',f12.3)') fstep,volume,vol_slv
         IF(node .EQ. 0) WRITE(*,1000) volume+vol_slv,volume_co
         IF(map .NE. 0) THEN
            DO i1=1,map,4
               n=3
               IF(map-i1 .LT. 3) n=map-i1
               WRITE(kvoronoi,'(a1,f11.2,1x,a9,4(f10.3,1x,i4))') 'T'
     &              ,fstep,'VolMolSlt',(vol_slt(j),j,j=i1,i1+n)
            END DO
         END IF

*=======================================================================
*--- Compute Voronoi Volume for types of resisues ----------------------
*=======================================================================

         nvol_type=0
         DO i=1,nbun
            IF(nvol_type .LT. mend(i)) nvol_type=mend(i)
         END DO

         DO i=1,nvol_type
            vol_type(i)=0.0D0
            vol_typep(i)=0
         END DO
         DO i=1,nbun
            vol_typep(mend(i))=vol_typep(mend(i))+1
         END DO
         DO i=1,nato
            ia=res2(i)
            vol_type(ia)=vol_type(ia)+volume_vor(i)
         END DO

         map=0
         DO ia=1,nvol_type
            IF(vol_typep(ia) .NE. 0) THEN
               map=map+1
               vol_type(map)=vol_type(ia)/DBLE(vol_typep(ia))
               vol_typep(map)=ia
            END IF
         END DO
         IF(map .NE. 0) THEN
            DO i1=1,map,3
               n=2
               IF(map-i1 .LT. 2) n=map-i1
               WRITE(kvoronoi,'(a1,f11.2,1x,a7,3(f10.4,1x,a8))') 'T'
     &              ,fstep,'VolType',(vol_type(j),prsymb(vol_typep(j)
     &              ),j=i1,i1+n)
            END DO
         END IF

      END IF
#ifdef PARALLEL
      CALL P_barrier
#endif

*=======================================================================
*--- Compute Voronoi interface between solute and solvent---------------
*=======================================================================

      IF(voronoi_access) THEN
         o=0
         DO ia=nstart_uh,nend_uh
            area_slt(ia)=0.0D0
            area_slv(ia)=0.0D0
            DO ib=mres(1,ia),mres(2,ia)
               DO i1=grppt(1,ib),grppt(2,ib)
                  o=o+1
                  typei=ss_index(i1)
                  IF(typei .EQ. 1) THEN
                     m=nnlpp_vor(1,o)
                     DO j=1,m
                        j1=nnlpp_vor(j+1,o)
                        typeij=ss_index(j1)+typei-1
                        IF(typeij .EQ. 1) THEN
                           
*---- Compute inter solute area 
                           
                           IF(res1(j1) .NE. ia) THEN
                              area_slt(ia)=area_slt(ia)+area_vor(j+1,o)
                           END IF
                           
                        ELSE IF(typeij .EQ. 2) THEN
                           
*---- Compute solute-solvent interface
                           
                           area_slv(ia)=area_slv(ia)+area_vor(j+1,o)
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         END DO

#ifdef PARALLEL
         CALL P_expand_r8(area_slt,nstart_uh,nend_uh,nlocal_uh,node
     &     ,nprocs)
         CALL P_expand_r8(area_slv,nstart_uh,nend_uh,nlocal_uh,node
     &     ,nprocs)
#endif
         map=0
         DO ia=1,nbun
            ib=mres(1,ia)
            i1=grppt(1,ib)
            typei=ss_index(i1)
            IF(typei .EQ. 1) THEN
               map=map+1
               index(map)=ia
               area_tot(map)=area_slt(ia)+area_slv(ia)
               IF(area_tot(ia) .EQ. 0.0D0) THEN
                  area_frac(map)=-1.0D0
               ELSE
                  area_frac(map)=area_slv(ia)/area_tot(ia)
               END IF
            END IF
         END DO
         IF(node .EQ. 0) THEN
            IF(voronoi_res(1) .EQ. 0) THEN
               IF(map .NE. 0) THEN
                  DO i1=1,map,3
                     n=2
                     IF(map-i1 .LT. 2) n=map-i1
                     WRITE(kvoronoi,'(a1,f11.2,1x,a4,3(f10.3,f10.3,i6)
     &                    )') 'T',fstep,'Area',(area_frac(ia)
     &                    ,area_tot(ia),index(ia),ia=i1,i1+n)
                  END DO
               END IF
            ELSE
               IF(map .NE. 0) THEN
                  DO i1=1,map
                     ok=.FALSE.
                     j1=1
                     DO WHILE(.NOT. ok .AND. j1 .LT. voronoi_res(i))
                        j1=j1+1
                        IF(index(i1) .EQ. voronoi_res(1+j1)) ok=.TRUE.
                     END DO
                     IF(ok) WRITE(kvoronoi
     &                    ,'(a1,f11.2,1x,a4,f10.3,f10.3,i6)')'T',fstep
     &                    ,'Area',area_frac(i1),area_tot(i1),index(i1)
                  END DO
               END IF
            END IF
         END IF
#ifdef PARALLEL
         CALL P_barrier
#endif
      END IF

*=======================================================================
*--- Compute Voronoi contact surface among molecules of the ------------
*--- solute ------------------------------------------------------------
*=======================================================================

      IF(nprocs .EQ. 1) THEN
         IF(voronoi_contact) THEN
            DO ic=1,ncontact_slt
               area_contact(ic)=0.0D0
            END DO
            DO ic=1,nprot
               area_prot(ic)=0.0D0
            END DO
            ia=ss_point(1)
            DO ib=1,ia
               i1=ss_point(1+ib)
               m=nnlpp_vor(1,i1)
               p1=atomp(i1)
               DO j=1,m
                  j1=nnlpp_vor(j+1,i1)
                  p2=atomp(j1)
                  DO ic=1,ncontact_slt
                     IF(p1 .EQ. contact_slt(1,ic) .AND. p2 .EQ.
     &                    contact_slt(2,ic)) THEN
                        
*---- Compute contact area between given molecules
                        
                        area_contact(ic)=area_contact(ic)+area_vor(j+1
     &                       ,i1)
                     END IF
                  END DO
                  IF(p1 .NE. p2) THEN
                     
*---- Compute total area for the molecule p1        
                     area_prot(p1)=area_prot(p1)+area_vor(j+1,i1)
                  END IF
               END DO
            END DO
            IF(map .NE. 0) THEN
               DO ic=1,ncontact_slt
                  WRITE(kvoronoi,'(a1,f11.2,a14,i6,a5,i6,3f10.3)') 'T'
     &                 ,fstep,' Contact Mol.',contact_slt(1,ic),' Mol '
     &                 ,contact_slt(2,ic),area_contact(ic)
     &                 ,area_prot(contact_slt(1,ic))
     &                 ,area_prot(contact_slt(2,ic))
               END DO
            END IF
         END IF
      END IF

*=======================================================================
*--- Compute number of neighbors among solute and solvent --------------
*=======================================================================

      IF(voronoi_neighbor) THEN
         mmap(1)=0
         mmap(2)=0
         mmap(3)=0
         o=0
         DO i1=nstart_ah,nend_ah
            o=o+1
            typei=ss_index(i1)
            m=nnlpp_vor(1,o)
            DO j=1,m
               j1=nnlpp_vor(j+1,o)
               typeij=typei+ss_index(j1)-1
               IF(typeij .EQ. 1) THEN
                  mmap(1)=mmap(1)+1
               ELSE IF(typeij .EQ. 3) THEN
                  mmap(3)=mmap(3)+1
               ELSE IF(typeij .EQ. 2) THEN
                  mmap(2)=mmap(2)+1
               END IF
            END DO
         END DO

#ifdef PARALLEL
         IF(nprocs .NE. 1) THEN
            CALL P_merge_r8(mmap(1))
            CALL P_merge_r8(mmap(2))
            CALL P_merge_r8(mmap(3))
         END IF
#endif
         IF(node .EQ. 0) WRITE(kvoronoi,'(a1,f11.2,a16,i7,a9,i7,a9,i7)')
     &        'T',fstep,' Neigh. SLT-SLT ',mmap(1),' SLV-SLV
     &        ',mmap(3),' SLT-SLV ',mmap(2)

#ifdef PARALLEL
         CALL P_barrier
#endif
      END IF
1000  FORMAT('|***** TotVol. Check: Voronoi ',f12.3,
     &        ' Traj. = ',f12.3,'*****|') 
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
