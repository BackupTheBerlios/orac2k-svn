
        if(polar) then  
           crij(1) = xc
           crij(2) = yc
           crij(3) = zc
           arg = rsp/plrzbij(lij)
           fact = 1.0d0 + arg + 0.50d0*arg**2
           cexp = exp(-arg)

           do k=1,3
              TDi(k) = 0.0d0
              TDj(k) = 0.0d0
              do l=1,3
                 Tij(k,l) = fudge*3.0d0*crij(k)*crij(l)*rspqi*rsqi
                 Tij(k,l) = Tij(k,l)*(1.0d0-(fact+arg**3/6.0d0)*cexp)
              end do
              Tij(k,k) = Tij(k,k) - fudge*(1.0d0-fact*cexp)*rspqi

              TDj(k) = TDj(k) + Tij(k,1)*dipx(j)
              TDj(k) = TDj(k) + Tij(k,2)*dipy(j)
              TDj(k) = TDj(k) + Tij(k,3)*dipz(j)
              TDi(k) = TDi(k) + Tij(k,1)*dipx(i1)
              TDi(k) = TDi(k) + Tij(k,2)*dipy(i1)
              TDi(k) = TDi(k) + Tij(k,3)*dipz(i1)

           end do

           Edx(i1) = Edx(i1) + TDj(1)
           Edy(i1) = Edy(i1) + TDj(2)
           Edz(i1) = Edz(i1) + TDj(3)
           Edx(j) = Edx(j) + TDi(1)
           Edy(j) = Edy(j) + TDi(2)
           Edz(j) = Edz(j) + TDi(3)

           auxd = dipx(i1)*TDj(1)+dipy(i1)*TDj(2)+dipz(i1)*TDj(3)
           Udd  = Udd  - auxd

           auxi = dipx(j)*crij(1)+dipy(j)*crij(2)+dipz(j)*crij(3)
           auxj = dipx(i1)*crij(1)+dipy(i1)*crij(2)+dipz(i1)*crij(3)
           enedip(i1) = enedip(i1) + auxi*rspqi
           enedip(j)  = enedip(j)  - auxj*rspqi

        end if
