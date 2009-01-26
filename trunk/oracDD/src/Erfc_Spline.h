/* Erfc_Spline.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _ERFC_SPLINE_H
#define _ERFC_SPLINE_H

#define MyErfc    Erfc

#define _Third         0.333333333333333333_8
#define _Half          0.50_8
#define _Erfc(i,h)     c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)*_Third)*_Half)
#define _DErfc(i,h)    c(2,i)+h*(c(3,i)+h*c(4,i)*_Half)

#endif
