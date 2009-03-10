#define PRESSURE
                        lij=type(nbti,nbtype(j))
                        xg=xd1-xp0(j)
                        yg=yd1-yp0(j)
                        zg=zd1-zp0(j)
                        xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                        yc=           co(2,2)*yg+co(2,3)*zg
                        zc=                      co(3,3)*zg
                        rsq=xc*xc+yc*yc+zc*zc
                        rsqi=1.0d0/rsq
                        r6=rsqi*rsqi*rsqi
                        r12=r6*r6
                        ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)
     &                       *r6
                        qforce=ssvir*rsqi
                        emvir = ssvir*rsqi
                        
                        conf=ecc12(lij)*r12-ecc6(lij)*r6
                        uconfa=uconfa+conf
                        
                        fppx(i1)=fppx(i1)+qforce*xc
                        fppy(i1)=fppy(i1)+qforce*yc
                        fppz(i1)=fppz(i1)+qforce*zc
                        fppx(j)=fppx(j)-qforce*xc
                        fppy(j)=fppy(j)-qforce*yc
                        fppz(j)=fppz(j)-qforce*zc
#ifdef PRESSURE
                        qfx=emvir*xc
                        qfy=emvir*yc
                        qfz=emvir*zc
                        st1 = st1+qfx*xg
                        st2 = st2+qfx*yg
                        st3 = st3+qfx*zg
                        st4 = st4+qfy*xg
                        st5 = st5+qfy*yg
                        st6 = st6+qfy*zg
                        st7 = st7+qfz*xg
                        st8 = st8+qfz*yg
                        st9 = st9+qfz*zg
#endif
