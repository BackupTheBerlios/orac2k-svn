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

#define __correct_vp(dt,fp)				\
  Do n=1,natom_P;							\
     tfact=_Half*dt*xmass(n);						\
     v0(n)%x=v0(n)%x+tfact*fp(n) % x;					\
     v0(n)%y=v0(n)%y+tfact*fp(n) % y;					\
     v0(n)%z=v0(n)%z+tfact*fp(n) % z;					\
  End Do;

#define __verlet_vp(dt,fp)						\
  Do n=1,natom_P;								\
     tfact=_Half*dt*xmass(n);						\
     v0(n)%x=v0(n)%x+tfact*fp(n)%x;		                        \
     p0(n)%x=p0(n)%x+dt*v0(n)%x;			                \
     v0(n)%y=v0(n)%y+tfact*fp(n)%y;		                        \
     p0(n)%y=p0(n)%y+dt*v0(n)%y;		                        \
     v0(n)%z=v0(n)%z+tfact*fp(n)%z;		                        \
     p0(n)%z=p0(n)%z+dt*v0(n)%z;		                        \
  End Do;


#endif
