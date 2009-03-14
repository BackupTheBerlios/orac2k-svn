
      INTEGER mapnl(m8),mapnl_slv(nmap_slv)
      INTEGER na,nb,nc,nd,ne,nf,ng,la
      PARAMETER(na=tsites, nb=tsitep, nc=elsiz, nd=tgroup, ne=hoov,
     x          nf=elsiz, ng=npm, la=1)
      REAL*8  xp0(nb),yp0(nb),zp0(nb),xpcm(ng),ypcm(ng),zpcm(ng),xpg(nd)
     &     ,ypg(nd),zpg(nd),gh(hoov)
      LOGICAL mapp(m1),mapp_slv(slvatm)
      SAVE xp0,yp0,zp0
