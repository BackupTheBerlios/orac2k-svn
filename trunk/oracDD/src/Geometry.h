/* Geometry.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#define __subtract_v(v1,v2,vout)		\
  vout%x=v2%x-v1%x;				\
  vout%y=v2%y-v1%y;				\
  vout%z=v2%z-v1%z;				

#define __add_v(v1,v2,vout)			\
  vout%x=v2%x+v1%x;				\
  vout%y=v2%y+v1%y;				\
  vout%z=v2%z+v1%z;				

#define __crossprod(v1,v2,vout)			\
  vout%x=v1%y*v2%z-v1%z*v2%y;			\
  vout%y=v1%z*v2%x-v1%x*v2%z;			\
  vout%z=v1%x*v2%y-v1%y*v2%x;

#define __dotprod(v1,v2,vout)			\
  vout=v1%x*v2%x+v1%y*v2%y+v1%z*v2%z;

#define __norm_v(v1)      			\
  SQRT(v1%x**2+v1%y**2+v1%z**2)

#define __normalize_v(v1,a)      			\
  a=__norm_v(v1);					\
  v1%x=v1%x/a;						\
  v1%y=v1%y/a;						\
  v1%z=v1%z/a;
    
#define __prod_v(m,v1,vout)				\
  vout%x=m*v1%x;					\
  vout%y=m*v1%y;					\
  vout%z=m*v1%z;
 
#endif
