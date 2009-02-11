/* Integrate.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _INTEGRATE_H
#define _INTEGRATE_H
#define _Half    0.5_8

#define __correct(dt,fp)				\
  Do nn=1,Size(IndBox_a_p);				\
  m=IndBox_a_p(nn);							\
  tfact=_Half*dt/Atoms(m)%mass;					\
  Atoms(m)%vx=Atoms(m)%vx+tfact*fp(m) % x;					\
  Atoms(m)%vy=Atoms(m)%vy+tfact*fp(m) % y;					\
  Atoms(m)%vz=Atoms(m)%vz+tfact*fp(m) % z;					\
  End Do;

#define __verlet(dt,dp)				\
  Do nn=1,Size(IndBox_a_p);			\
  m=IndBox_a_p(nn);							\
  tfact=_Half*dt/Atoms(m)%mass;					\
  Atoms(m)%vx=Atoms(m)%vx+tfact*fp_n0(m)%x;		                        \
  Atoms(m)%x=Atoms(m)%x+dt*Atoms(m)%vx;			                \
  Atoms(m)%vy=Atoms(m)%vy+tfact*fp_n0(m)%y;		                        \
  Atoms(m)%y=Atoms(m)%y+dt*Atoms(m)%vy;		                        \
  Atoms(m)%vz=Atoms(m)%vz+tfact*fp_n0(m)%z;		                        \
  Atoms(m)%z=Atoms(m)%z+dt*Atoms(m)%vz;		                        \
  End Do;

#endif
