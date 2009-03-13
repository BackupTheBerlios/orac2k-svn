#ifdef __Fast
                           lij=type(nbti,nbtype(j))
                           xg=xpi-xp0(j)-xmap1(jj)
                           yg=ypi-yp0(j)-ymap1(jj)
                           zg=zpi-zp0(j)-zmap1(jj)
                           xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                           yc=           co(2,2)*yg+co(2,3)*zg
                           zc=                      co(3,3)*zg
                           rsq=xc*xc+yc*yc+zc*zc
                           rsqi=1.0D0/rsq


                           qforce=0.0D0; emvir=0.0D0; aux1=0.0D0
                           furpar= chrgei*charge(j)
                           ok_LJ=rsq < eccc(lij)
                           ok_Co=(rsq < Ewald_cut .AND. ABS(furpar) >
     &                          Tol_q)
              
                           IF((.NOT. ok_LJ) .AND. (.NOT. ok_Co)) CYCLE

                           IF(ok_LJ) THEN
                              r6=rsqi*rsqi*rsqi
                              r12=r6*r6
                              ssvir=12.0d0*ecc12(lij)*r12-6.0d0
     &                             *ecc6(lij)*r6
                              qforce=ssvir*rsqi*swrs(jj)
                              conf=ecc12(lij)*r12-ecc6(lij)*r6
                              cmap2(jj)=cmap2(jj)+conf
                              uconf(typeij)=uconf(typeij)+swrs(jj)*conf
                           End If
                           
                           IF(ok_Co) THEN
                              rsp=DSQRT(rsq)
                              rspi=one/rsp
                              rspqi=rsqi*rspi
                              alphar=alphal*rsp
                           
                              xx=alphar
                              ind = erftbdns*xx + 1
                              dxx = xx - (ind-1)*erfc_bin
                              aux1=dxx*(erfc_arr(3,ind)+dxx*erfc_arr(4
     &                             ,ind)*threei)
                              derfcst = -(erfc_arr(2,ind)+dxx
     &                             *(erfc_arr(3,ind) + dxx*erfc_arr(4
     &                             ,ind)*0.5d0))
                              erfcst = erfc_arr(1,ind)+dxx*(erfc_arr(2
     &                             ,ind)+aux1*twoi)

                              ucoul(typeij)=ucoul(typeij)+swrs(jj)
     &                             *furpar*erfcst*rspi
                              aux1  = furpar*(erfcst+alphar*derfcst)
     &                             *rspqi
                              qforce= qforce+aux1*swrs(jj)
                              cmap2(jj)=cmap2(jj)+furpar*erfcst*rspi
                           End If

                           emvir = qforce
                           
                           fppx(i1)=fppx(i1)+qforce*xc
                           fppy(i1)=fppy(i1)+qforce*yc
                           fppz(i1)=fppz(i1)+qforce*zc
                           fppx(j)=fppx(j)-qforce*xc
                           fppy(j)=fppy(j)-qforce*yc
                           fppz(j)=fppz(j)-qforce*zc
                           phi(i1)=phi(i1)+swrs(jj)*charge(j)*erfcst*rspi
                           phi(j)=phi(j)+swrs(jj)*chrgei*erfcst*rspi
#ifdef PRESSURE
                           st1 = st1+emvir*xc*xg
                           st2 = st2+emvir*xc*yg
                           st3 = st3+emvir*xc*zg
                           st4 = st4+emvir*yc*xg
                           st5 = st5+emvir*yc*yg
                           st6 = st6+emvir*yc*zg
                           st7 = st7+emvir*zc*xg
                           st8 = st8+emvir*zc*yg
                           st9 = st9+emvir*zc*zg
#endif
#else
                           lij=type(nbti,nbtype(j))
                           xg=xpi-xp0(j)-xmap1(jj)
                           yg=ypi-yp0(j)-ymap1(jj)
                           zg=zpi-zp0(j)-zmap1(jj)
                           xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                           yc=           co(2,2)*yg+co(2,3)*zg
                           zc=                      co(3,3)*zg
                           rsq=xc*xc+yc*yc+zc*zc
                           rsqi=1.0D0/rsq
                           r6=rsqi*rsqi*rsqi
                           r12=r6*r6
                           ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)
     &                          *r6
                           qforce=ssvir*rsqi*swrs(jj)
                           conf=ecc12(lij)*r12-ecc6(lij)*r6
                           cmap2(jj)=cmap2(jj)+conf
                           uconf(typeij)=uconf(typeij)+swrs(jj)*conf
                           rsp=DSQRT(rsq)
                           rspi=one/rsp
                           rspqi=rsqi*rspi
                           alphar=alphal*rsp
                           
                           xx=alphar
                           ind = erftbdns*xx + 1
                           dxx = xx - (ind-1)*erfc_bin
                           aux1=dxx*(erfc_arr(3,ind)+dxx*erfc_arr(4,ind)
     &                          *threei)
                           derfcst = -(erfc_arr(2,ind)+dxx*(erfc_arr(3
     &                          ,ind) + dxx*erfc_arr(4,ind)*0.5d0))
                           erfcst = erfc_arr(1,ind)+dxx*(erfc_arr(2,ind)
     &                          +aux1*twoi)
                           furpar= chrgei*charge(j)
                           ucoul(typeij)=ucoul(typeij)+swrs(jj)*furpar
     &                          *erfcst*rspi
                           aux1  = furpar*(erfcst+alphar*derfcst)*rspqi
                           qforce= qforce+aux1*swrs(jj)
                           cmap2(jj)=cmap2(jj)+furpar*erfcst*rspi
                           emvir = qforce
                           
                           fppx(i1)=fppx(i1)+qforce*xc
                           fppy(i1)=fppy(i1)+qforce*yc
                           fppz(i1)=fppz(i1)+qforce*zc
                           fppx(j)=fppx(j)-qforce*xc
                           fppy(j)=fppy(j)-qforce*yc
                           fppz(j)=fppz(j)-qforce*zc
                           phi(i1)=phi(i1)+swrs(jj)*charge(j)*erfcst*rspi
                           phi(j)=phi(j)+swrs(jj)*chrgei*erfcst*rspi
#ifdef PRESSURE
                           st1 = st1+emvir*xc*xg
                           st2 = st2+emvir*xc*yg
                           st3 = st3+emvir*xc*zg
                           st4 = st4+emvir*yc*xg
                           st5 = st5+emvir*yc*yg
                           st6 = st6+emvir*yc*zg
                           st7 = st7+emvir*zc*xg
                           st8 = st8+emvir*zc*yg
                           st9 = st9+emvir*zc*zg
#endif
#endif
