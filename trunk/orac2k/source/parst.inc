
#include "config.h"


*------------------ PARAMETER  STATEMENT -------------------------------

      INTEGER sitslv,sitslu,types,typep,tsites,tsitep,
     &     elsiz,secsiz,mfold,nores,tgroup,hoov,npm,mspline,mpp,nbyte
     &     ,rbyte,types_gofr,pcontact_slt,ntopol,n_mat,cheb_order,maxsk
      PARAMETER(sitslv=_SIT_SOLV_,
     &          sitslu=_SIT_SOLU_,
     &          types=_TYP_SOLV_,
     &          types_gofr=_TYP_SOLV_GOFR_,
     &          typep=_TYP_SOLU_,
     &          nores=_NRES_,
     &          tgroup=_TGROUP_,
     &          hoov=_HOOVER_,
     &          npm=_N_PROT_MAX_,
     &          mspline=_DIR_SPLINE_,
     &          nbyte=_NBYTE_,
     &          rbyte=_RBYTE_,
     &          pcontact_slt=_CONTACT_VORONOI_,
     &          ntopol=_ANL_TOPO_,
     &          n_mat=_DYNAMIC_DIM_,
     &          cheb_order=_CHEB_ORDER_,
     &          tsites=sitslv,tsitep=sitslu,
     &          elsiz=1,secsiz=1)
      INTEGER n1,n2,n3,n4,n5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13
     &     ,m14,maxint,g1,g2,f1,f2,f3,ndan1,mex,slv2,slv3,slv4,slv5,slv6
     &     ,slv_cns1,slv_cns2,slvunit,slvsit,slvatm,slv_group,nmap_slv
     &     ,m15

      PARAMETER(n1=sitslv,n2=types*(types+1)/2,n3=10,n4=tsites,
     &          n5=_TOT_RES_,
     &     m1=tsitep,m2=m1*3,m3=3*m1,m4=m1/3,m5=typep*(typep+1)/2,
     &          m6=typep,m7=2,m8=tsitep*10,
     &          maxint=_MAX_BIN_,
     &          maxsk=_MAX_BIN_*5,
     &          g1=types_gofr*(types_gofr+1)/2+3,
     &          g2=100,m9=m1+600,
     &          m10=_CONN_,
     &          m11=tgroup,m12=tgroup*2,
     &          m13=2,m14=2,mfold=10,f1=_F1_,f2=_F2_,
     &          f3=_F3_,ndan1=_SUBUN_, mpp=_NAT_PP_*tgroup,
     &          mex=20*sitslu,
     &          m15=2*tgroup+tgroup/2)

      PARAMETER(slvsit=_SIT_SOLV_,
     &          slvunit=_MOL_NRES_,
     &          slvatm=slvsit*slvunit,
     &          slv2=slvatm+20,
     &          slv3=4*slvatm,
     &          slv4=4*slvatm,
     &          slv5=slvatm,
     &          slv6=_CONN_,
     &          slv_cns1=10,
     &          slv_cns2=2,
     &          slv_group=_SOLV_GROUP_,
     &          nmap_slv=tsites*10)
