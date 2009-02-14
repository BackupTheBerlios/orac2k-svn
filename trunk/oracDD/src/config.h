/* config.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */

#ifndef _CONFIG_H
#define _CONFIG_H

/*
  Define type of pointt-point communications: 
        __Blocking       using mpi_sendrecv  
	__NonBlocking    using mpi_isend and mpi_irecv
	
#define __Blocking
*/

#define __NonBlocking

#endif
