/* BoxGeometry.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _BOXGEOMETRY_H
#define _BOXGEOMETRY_H

#define __MaxSide	      6
#define _x_        1
#define _y_        2
#define _z_        3

#define __Segment(n)					\
  0.5D0*(1.0D0+DBLE(n));

#define __Fraction				\
  0.10_8

#define __Half					\
  0.50_8
#define __Pbc(p,pout)				\
  pout%x=-2.0D0*ANINT(0.5D0*p%x);		\
  pout%y=-2.0D0*ANINT(0.5D0*p%y);		\
  pout%z=-2.0D0*ANINT(0.5D0*p%z);	       
#endif
