      SUBROUTINE rdcmac(pdb_unit,kprint,res,beta,x0,y0,z0,nato
     &     ,multiple,mol,nmol,iret,errmsg)
************************************************************************
*      Read pdb file                                                   *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER pdb_unit,kprint,nato,res(nato),iret,mol,nmol
      REAL*8 x0(*),y0(*),z0(*)
      CHARACTER*7 beta(*)
      CHARACTER*80 errmsg
      LOGICAL multiple

      include 'parst.h'
      INTEGER  mmm
      parameter     (mmm=nores+slvunit) 

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8   xtp,ytp,ztp,xa,ya,za
      INTEGER  ires,i,ires_old,res_sequence,pdb_res,n,nwx,nwy,nwz,iline
     &     ,nword,idot,dot(3),pdb_res_old,ibeg(mmm),iend(mmm),jatm,j,is
     &     ,nwx_first,nwy_first,nwz_first,ibgdt(3),ieddt(3),nerror,i1,i2
     &     ,nh,pdb_res_count
      LOGICAL*1     lok,found,ltmp(sitslu+tsites)
      LOGICAL ok
      CHARACTER*10  integers
      CHARACTER*5   field(3)
      CHARACTER*4   beta4,char 
      CHARACTER*80 line
      CHARACTER*1  beta1

*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      COMMON /rag2/ ibeg,iend,ltmp

c=======================================================================
c---  finds the start-of-residue end-of-residue pointers in 
c---  topologic array. Also find total expected number of residue 
c=======================================================================

      ires = 0
      res_sequence=0
      nh = 0
      jatm=0
      do i=1,nato
c---     this coordinates are overwritten if hydrogen are in the pdb
         beta1=beta(i)
	 if(beta1.eq.'h') then 
	     nh = nh + 1
	     x0(i)=1.d10 
	     y0(i)=1.d10 
	     z0(i)=1.d10
	 end if 
         ltmp(i)=.false.
         ires_old = ires
         ires=res(i)
         if(ires.ne.ires_old) then 
            res_sequence = res_sequence + 1
            ibeg(res_sequence) = i 
            if (i.gt.1) THEN 
               iend(res_sequence-1) = i-1
            END IF
         END IF
      END DO 

      iend(res_sequence) = nato

      pdb_res= 1
      pdb_res_count=1
      ok=.FALSE.

c======================================================================
c     Start parsing the PDB file 
c======================================================================

      nwx=0
      nwy=0
      nwz=0
      ibgdt(1)=31
      ieddt(1)=38
      ibgdt(2)=39
      ieddt(2)=46
      ibgdt(3)=47
      ieddt(3)=54
      field(1) = '31-38'
      field(2) = '39-46'
      field(3) = '47-45'
      integers = '0123456789' 
      nerror=0
       
      n=0
      iline=0

      IF(mol .EQ. 1) WRITE(kprint,'(5x,a)')
     &     'Reading the PDB file               ---->'

 1    iline = iline + 1 
c-----   read line and bring to lowercase
         READ(pdb_unit,'(a80)',END=2,ERR=1001) line
         CALL up_low(line,80)
         IF(line(1:4) .eq. 'atom' .or. line(1:6) .eq. 'hetatm'
     &        ) THEN 

c======================================================================
c           Check line to see if .pdb 
c======================================================================

c===========check position of (.) character for coordinates============

            idot=0
            dot(1) = 0
            dot(2) = 0
            dot(3) = 0

c---        if position are wrong, abort

            do is=31,54
               if(line(is:is).eq.'.') then 
                  idot = idot+1 
                  if(idot.gt.3) THEN 
                     WRITE(kprint,100) iline
                     nerror=nerror+1       
                  end if 
                  dot(idot) = is
               endif   
            end do

            if(idot.lt.3) THEN 
               WRITE(kprint,210) iline
               nerror=nerror+1       
            end if 

c---        if position are misaligned give warning 

            do i=1,3
               lok = dot(i).ge.ibgdt(i).and.dot(i).le.ieddt(i)
               IF (.not.lok) THEN 
                  WRITE(kprint,130) iline,field(i)
                  nerror=nerror+1       
               END IF
            END DO
            IF(line(35:35).ne.'.') THEN
               nwx = nwx + 1
               if (nwx.eq.1) nwx_first=iline
            ENDIF 
            IF (line(43:43).ne.'.') THEN 
               nwy = nwy + 1
               if (nwy.eq.1) nwy_first=iline
            ENDIF 
            IF (line(51:51).ne.'.') THEN 
               nwz = nwz + 1
               if (nwz.eq.1) nwz_first=iline
            ENDIF 

c===========check if field 26-26 is a number================================= 

            lok = .false.
            do i=1,10
               if(line(26:26).eq.integers(i:i)) THEN
                  lok = .true. 
               end if
            end do
            IF (.not.lok) THEN 
               WRITE(kprint,120) iline
               STOP
            END IF

            if(nerror.gt.100) THEN 
               write(kprint,310) 
               STOP
            END IF
               
c----       tests passed: Good chance that line is in genuin pdb format  

c======================================================================
c           Now do internal read of LINE in PBD format
c======================================================================

            IF(ok) pdb_res_old = pdb_res
            READ(line,4,err=3) char,pdb_res,xtp,ytp,ztp 
            IF(.NOT. ok) pdb_res_old=pdb_res
            ok=.TRUE.

            IF (pdb_res .NE. pdb_res_old)THEN
               pdb_res_count=pdb_res_count+1
            END IF
            
c----       check if .pdb residue number is something meaningful

            if(pdb_res_count .GT. res_sequence) then
               IF(.NOT. multiple) THEN
                  write(kprint,10) iline,pdb_res_count
                  STOP
               END IF
               BACKSPACE pdb_unit
               GOTO 2
            END IF

c--         justify left character string 

            call jusleft(char)

c---        Now we check the Residue sequence of the .pdb file. 
c----       check if CHAR matches some BETA in residue ires

            found = .false.
            do j=ibeg(pdb_res_count),iend(pdb_res_count) 
               beta4=beta(j)
               if (char.eq.beta4) then 
                  found = .true. 

c----             j-th atom of residue pdb_res matches!!!!

                  if(.not.ltmp(j)) then 
                     ltmp(j)=.true.
                   else

c----                 unique labels must be specified in the .pdb 

                      WRITE(kprint,75) iline,char
                      STOP
                   END IF 
                  go to 101
               end if
            end do

c---        AAAAArgh label not found: write error message

            WRITE(kprint,20) iline,char
            STOP
 101        continue
            if (found) then 
               n = n + 1 
               x0(j)=xtp
               y0(j)=ytp
               z0(j)=ztp
            end if   

c---        if the atom was an hydrogen decrement number of h atoms in
c---        tipological sequence and increment nh of pbd sequence            

            beta1=char
            if(beta1.eq.'h') then 
               nh = nh-1
            end if
            go to 1
         END IF
         if(pdb_res_count .lt.res_sequence) go to 1
      CONTINUE

 2    continue

c--   check for bad alignment 

      IF(nerror.gt.0) THEN
         WRITE(kprint,66)  nerror
         STOP
      END IF

c---  check number of atoms 
      
      n = n + nh 

      IF (n.ne.nato.and.n.gt.0) THEN  
         WRITE(kprint,40) pdb_res_count, res_sequence
         STOP
      END IF

      IF (n.eq.0) THEN
         WRITE(kprint,50) 
         STOP 
      END IF

c======================================================================
c     Atoms number is OK; trasform coordinates and exit; 
c======================================================================

      IF(mol .EQ. nmol) WRITE(kprint,290) n*nmol,nh

c======================================================================
c     Signals misaligned floating point in .pdb coordinates
c======================================================================

c---  jump here because of internal reading error; set iret to 1
c---  prints lines where possible misalignment of reals could 
c---  have been found and exit 

      IF(nwx.gt.0) THEN
         WRITE (kprint,65) 
         WRITE(kprint,70) nwx,nwx_first
         WRITE (kprint,65) 
      END IF
      IF(nwy.gt.0) THEN
         WRITE (kprint,65) 
         WRITE(kprint,76) nwy,nwy_first
         
         WRITE (kprint,65) 
      END IF
      IF(nwz.gt.0) THEN
         WRITE (kprint,65) 
         WRITE(kprint,80) nwz,nwz_first
         
         WRITE (kprint,65) 
      END IF

      RETURN

 3    WRITE(kprint,55) iline
      STOP

 1001 WRITE(kprint,60) iline
*======================= END statements  ================================
      STOP
4     FORMAT(12x,a4,6x,i4,4x,3f8.3)
10    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,
     &     '; RESIDUE NUMBER IS ',i5,' ?????'/)   
20    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,
     &     '; LABEL ',A4,' UNKNOWN'/)
40    FORMAT(//'****ERROR IN PDB FILE:',i4,' residues found while ',
     &     'expecting ',i4,/,4x,'residues from sequence given in JOIN',
     &     '(&PARAMETERS)')
50    FORMAT(//'****ERROR IN PDB FILE: no ATOM or HETATM',
     &     ' string found')
55    FORMAT(/'****ERROR in PDB: internal reading error at line:',i5)
60    FORMAT('****ERROR in reading PDB at line:',i5) 
65    FORMAT('***************************************************')
66    FORMAT( / '-----  Stop while reading the PDB file:',    
     &     / '----- ',i3,' fatal errors were found')
70    FORMAT('***WARNING:IN PDB FILE X coordinates are misaligned:'
     &     , /,'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned X, starting at line: ',i5)     
      
75    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,
     &     '; ATOM ',a4, ' already assigned')
76    FORMAT('***WARNING:IN PDB FILE Y coordinates are misaligned:'
     &     ,/,'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned Y, starting at line: ',i5)     
80    FORMAT('***WARNING:IN PDB FILE Z coordinates are misaligned:'
     &     ,/'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned Z, starting at line: ',i5)     
100   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' unexpected FORMAT: more than 3 reals')
120   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' - at column 26 a number is expected')
130   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' Character (.) not found in field ', A5)
210   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' unexpected FORMAT: less  than 3 reals')
280   FORMAT
     &     (/'                    Reading PDB file.......') 
310   FORMAT(// ' --- more than 100 errors. Check alignment'/)
      
290   FORMAT
     &     (5x,'PDB file parsed:',i6,' atoms; - ',i6,' hydrogens to be'
     &     ,' assigned')

      END
************************************************************************
      SUBROUTINE jusleft(STRING)
      implicit none 
C     justify to the left 4 character string  
*======================= DECLARATIONS ==================================
      CHARACTER*4 STRING
      INTEGER  i,istart 
*======================= Exec statements  ==============================
      i=0
      istart=1
 1    i=i+1
      if (string(i:i).ne.' ') then
         istart = i
         go to 2
      end if
      go to 1
2     continue
      string = string(istart:4)		   	 	 

*======================= END statements  ================================

      RETURN
      END 
