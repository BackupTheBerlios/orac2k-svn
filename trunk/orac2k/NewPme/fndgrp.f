      SUBROUTINE fndgrp(nprot,ngrp,protl,grppt,atomg,protg,groupp,atomp
     &     ,ndim,iret,errmsg)

************************************************************************
*                                                                      *
*          => nprot,ngrp,prot,grppt                                    *
*         <=  protg,groupp,atomp                                       *
*                                                                      *
*                                                                      *
*          protg(1,i)   primo e...                                     *
*          protg(2,i)   ...ultimo gruppo della proteina i              *
*                                                                      *
*          groupp(i)  proteina a cui appartiene il gruppo i            *
*                                                                      *
*          atomp(i)   proteina a cui appartiene l'atomo i              *
*                                                                      *
************************************************************************

           IMPLICIT none
           INCLUDE 'parst.h'
           INTEGER i,j,k,m,n,nprot,ngrp,iret,ndim,j1,j2,count,mg
     &          ,count2,i1
           INTEGER protg(*),protl(*),grppt(2,*),groupp(*),atomp(*)
     &          ,atomg(*)
           CHARACTER*80 errmsg

           iret=0.0D0
           count=0
           count2=0
           DO i=1,nprot
              m=protl(count+1)
              mg=0
              j2=0
              DO i1=1,m
                 j=protl(count+1+i1)
                 j1=atomg(j)
                 IF(i1 .NE. 1) THEN
                    IF(j1 .NE. j2) THEN
                       mg=mg+1
                       protg(count2+1+mg)=j1
                    END IF
                 ELSE IF(i1 .EQ. 1) THEN
                    mg=mg+1
                    protg(count2+1+mg)=j1
                 END IF
                 IF(i1 .EQ. m) THEN
                    protg(1+count2)=mg
                    count2=count2+mg+1
                    If(count2 > m15) Then
                       iret=1
                       errmsg='Protg underdimensioned, increase m15 in '
     &                      //'parst.h.  Abort '
                       Return
                    End If
                 END IF
                 j2=j1
              END DO
              count=count+1+m
           END DO
           count=0
           DO i=1,nprot
              m=protg(count+1)
              count=count+m+1
           END DO

           count=0
           DO i=1,nprot
              m=protl(count+1)
              count=count+m+1
           END DO

           count=0
           DO i=1,nprot
              m=protg(count+1)
              DO i1=1,m
                 j=protg(count+1+i1)
                 groupp(j) = i 
              END DO
              count=count+m+1
           END DO

           count=0
           DO i=1,nprot
              m=protl(count+1)
              DO i1=1,m
                 j=protl(count+1+i1)
                 atomp(j)=i
              END DO
              count=count+m+1
           END DO

           RETURN
           END
