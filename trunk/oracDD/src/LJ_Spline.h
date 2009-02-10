/* LJ_Spline.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _LJ_SPLINE_H
#define _LJ_SPLINE_H

#define  __dlj(r,vout)			\
  vout=-(12.0D0*(1.0D0/(r))**12-6.0D0*(1.0D0/(r))**6)/(r);

#define  __lj(r,vout)			\
  vout=(1.0D0/(r))**12-(1.0D0/(r))**6;

#define _Third         0.333333333333333333_8
#define _Half          0.50_8

#define _Lj(i,h)     (s(1,i)+h*(s(2,i)+h*(s(3,i)+h*s(4,i)*_Third)*_Half))
#define _DLj(i,h)    (s(2,i)+h*(s(3,i)+h*s(4,i)*_Half))

#endif
