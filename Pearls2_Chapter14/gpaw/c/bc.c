/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Copyright (C) 2005       CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

// Copyright (C) 2003  CAMP
// Please see the accompanying LICENSE file for further information.

#include <string.h>
#include <assert.h>
#include "bc.h"
#include "extensions.h"
#include <stdio.h>
#include <stdlib.h>


boundary_conditions* bc_init(const long size1[3],
           const long padding[3][2],
           const long npadding[3][2],
           const long neighbors[3][2],
           MPI_Comm comm, bool real, bool cfd)
{
  boundary_conditions* bc = GPAW_MALLOC(boundary_conditions, 1);

  for (int i = 0; i < 3; i++)
    {
      bc->size1[i] = size1[i];
      bc->size2[i] = size1[i] + padding[i][0] + padding[i][1];
      bc->padding[i] = padding[i][0];
    }

  bc->comm = comm;
  bc->ndouble = (real ? 1 : 2);
  bc->cfd = cfd;

  int rank = 0;
  if (comm != MPI_COMM_NULL)
    MPI_Comm_rank(comm, &rank);

  int start[3];
  int size[3];
  for (int i = 0; i < 3; i++)
    {
      start[i] = padding[i][0];
      size[i] = size1[i];
    }

  for (int i = 0; i < 3; i++)
    {
      int n = bc->ndouble;
      for (int j = 0; j < 3; j++)
      if (j != i)
        n *= size[j];

      for (int d = 0; d < 2; d++)
        {
          int ds = npadding[i][d];
          int dr = padding[i][d];
          for (int j = 0; j < 3; j++)
            {
              bc->sendstart[i][d][j] = start[j];
              bc->sendsize[i][d][j] = size[j];
              bc->recvstart[i][d][j] = start[j];
              bc->recvsize[i][d][j] = size[j];
            }
          if (d == 0)
            {
              bc->sendstart[i][d][i] = dr;
              bc->recvstart[i][d][i] = 0;
            }
          else
            {
              bc->sendstart[i][d][i] = padding[i][0] + size1[i] - ds;
              bc->recvstart[i][d][i] = padding[i][0] + size1[i];
            }
          bc->sendsize[i][d][i] = ds;
          bc->recvsize[i][d][i] = dr;

          bc->sendproc[i][d] = DO_NOTHING;
          bc->recvproc[i][d] = DO_NOTHING;
          bc->nsend[i][d] = 0;
          bc->nrecv[i][d] = 0;

          int p = neighbors[i][d];
          if (p == rank)
            {
              if (ds > 0)
                bc->sendproc[i][d] = COPY_DATA;
              if (dr > 0)
                bc->recvproc[i][d] = COPY_DATA;
            }
          else if (p >= 0)
            {
              // Communication required:
              if (ds > 0)
                {
                  bc->sendproc[i][d] = p;
                  bc->nsend[i][d] = n * ds;
                }
                    if (dr > 0)
                {
                  bc->recvproc[i][d] = p;
                  bc->nrecv[i][d] = n * dr;
                }
            }
        }

      if (cfd == 0)
        {
          start[i] = 0;
          size[i] = bc->size2[i];
        }
      // If the two neighboring processors along the
      // i'th axis are the same, then we join the two communications
      // into one:
      bc->rjoin[i] = ((bc->recvproc[i][0] == bc->recvproc[i][1]) &&
          bc->recvproc[i][0] >= 0);
      bc->sjoin[i] = ((bc->sendproc[i][0] == bc->sendproc[i][1]) &&
          bc->sendproc[i][0] >= 0);
    }

  bc->maxsend = 0;
  bc->maxrecv = 0;
  for (int i = 0; i < 3; i++)
    {
      int n = bc->nsend[i][0] + bc->nsend[i][1];
      if (n > bc->maxsend)
        bc->maxsend = n;
      n = bc->nrecv[i][0] + bc->nrecv[i][1];
      if (n > bc->maxrecv)
        bc->maxrecv = n;
    }

  return bc;
}


void bc_unpack1(const boundary_conditions* bc,
                const double* aa1, double* aa2, int i,
                MPI_Request recvreq[2],
                MPI_Request sendreq[2],
                double* rbuff, double* sbuff,
                const double_complex phases[2], int thd, int nin)
{

  int ng = bc->ndouble * bc->size1[0] * bc->size1[1] * bc->size1[2];
  int ng2 = bc->ndouble * bc->size2[0] * bc->size2[1] * bc->size2[2];
  bool real = (bc->ndouble == 1);
  for (int m = 0; m < nin; m++)
    // Copy data:
    if (i == 0)
      {
        // Zero all of a2 array.  We should only zero the bounaries
        // that are not periodic, but it's simpler to zero everything!
        // XXX
        memset(aa2 + m * ng2, 0, ng2 * sizeof(double));

        // Copy data from a1 to central part of a2:
        if (real)
          bmgs_paste(aa1 + m * ng, bc->size1, aa2 + m * ng2,
		     bc->size2, bc->sendstart[0][0]);
        else
          bmgs_pastez((const double_complex*)(aa1 + m * ng), bc->size1,
		      (double_complex*)(aa2 + m * ng2),
		      bc->size2, bc->sendstart[0][0]);
      }

#ifdef PARALLEL
  // Start receiving.
  for (int d = 0; d < 2; d++)
    {
      int p = bc->recvproc[i][d];
      if (p >= 0)
        {
          if (bc->rjoin[i])
            {
              if (d == 0)
                MPI_Irecv(rbuff, (bc->nrecv[i][0] + bc->nrecv[i][1]) * nin,
			  MPI_DOUBLE, p,
                          10 * thd + 1000 * i + 100000,
                          bc->comm, &recvreq[0]);
            }
          else
          {
            MPI_Irecv(rbuff, bc->nrecv[i][d] * nin, MPI_DOUBLE, p,
		      d + 10 * thd + 1000 * i,
                      bc->comm, &recvreq[d]);
	    rbuff += bc->nrecv[i][d] * nin;
          }
        }
    }
  // Prepare send-buffers and start sending:
  double* sbuf = sbuff;
  double* sbuf0 = sbuff;
  for (int d = 0; d < 2; d++)
    {
      sendreq[d] = 0;
      int p = bc->sendproc[i][d];
      if (p >= 0)
        {
          const int* start = bc->sendstart[i][d];
          const int* size = bc->sendsize[i][d];

	  for (int m = 0; m < nin; m++)
	    if (real)
	      bmgs_cut(aa2 + m * ng2, bc->size2, start,
		       sbuf + m * bc->nsend[i][d],
		       size);
	    else
	      bmgs_cutmz((const double_complex*)(aa2 + m * ng2),
			 bc->size2, start,
			 (double_complex*)(sbuf + m * bc->nsend[i][d]),
			 size, phases[d]);

          if (bc->sjoin[i])
            {
              if (d == 1)
                {
                  MPI_Isend(sbuf0, (bc->nsend[i][0] + bc->nsend[i][1]) * nin,
			    MPI_DOUBLE, p,
                            10 * thd + 1000 * i + 100000,
                            bc->comm, &sendreq[0]);
                }
            }
          else
            {
              MPI_Isend(sbuf, bc->nsend[i][d] * nin, MPI_DOUBLE, p,
                        1 - d + 10 * thd + 1000 * i, bc->comm, &sendreq[d]);
            }
          sbuf += bc->nsend[i][d] * nin;
        }
    }
#endif // Parallel
  for (int m = 0; m < nin; m++)
    {
      // Copy data for periodic boundary conditions:
      for (int d = 0; d < 2; d++)
        if (bc->sendproc[i][d] == COPY_DATA)
          {
            if (real)
              bmgs_translate(aa2 + m * ng2, bc->size2, bc->sendsize[i][d],
                 bc->sendstart[i][d], bc->recvstart[i][1 - d]);
            else
              bmgs_translatemz((double_complex*)(aa2 + m * ng2), bc->size2,
                   bc->sendsize[i][d],
                   bc->sendstart[i][d], bc->recvstart[i][1 - d],
                       phases[d]);
          }
    }
}


void bc_unpack2(const boundary_conditions* bc,
    double* a2, int i,
    MPI_Request recvreq[2],
    MPI_Request sendreq[2],
    double* rbuf, int nin)
{
#ifdef PARALLEL
  int ng2 = bc->ndouble * bc->size2[0] * bc->size2[1] * bc->size2[2];

  // Store data from receive-buffer:
  bool real = (bc->ndouble == 1);

  double* rbuf0 = rbuf;
  for (int d = 0; d < 2; d++)
    if (bc->recvproc[i][d] >= 0)
      {
        if (bc->rjoin[i])
          {
            if (d == 0)
              {
                MPI_Wait(&recvreq[0], MPI_STATUS_IGNORE);
                rbuf += bc->nrecv[i][1] * nin;
              }
            else
              rbuf = rbuf0;
	  }
	else
	  MPI_Wait(&recvreq[d], MPI_STATUS_IGNORE);
	
	for (int m = 0; m < nin; m++)
	  if (real)
	    bmgs_paste(rbuf + m * bc->nrecv[i][d], bc->recvsize[i][d],
		       a2 + m * ng2, bc->size2, bc->recvstart[i][d]);
	  else
	    bmgs_pastez((const double_complex*)(rbuf +
						m * bc->nrecv[i][d]),
			bc->recvsize[i][d],
			(double_complex*)(a2 + m * ng2),
			bc->size2, bc->recvstart[i][d]);
	rbuf += bc->nrecv[i][d] * nin;
      }
  
  // This does not work on the ibm with gcc!  We do a blocking send instead.
  for (int d = 0; d < 2; d++)
    if (sendreq[d] != 0)
      MPI_Wait(&sendreq[d], MPI_STATUS_IGNORE);
#endif // PARALLEL
}
