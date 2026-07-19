/* pyvertexenum.c - ctypes-friendly shim around lrslib for vertex enumeration.
 *
 * Original lrslib Copyright: David Avis 2003,2011 avis@cs.mcgill.ca
 * R interface written by Robert Robere, October 2012.
 * Python (ctypes) interface derived from the R interface, 2026.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lrslib.h"

#ifdef _WIN32
#define VE_EXPORT __declspec(dllexport)
#define VE_NULL_DEVICE "NUL"
#else
#define VE_EXPORT
#define VE_NULL_DEVICE "/dev/null"
#endif

/* status codes returned by ve_enumerate */
#define VE_OK 0
#define VE_ERR_INIT 1        /* lrs_init / lrs_alloc_dat / lrs_alloc_dic failed */
#define VE_ERR_NOMEM 2       /* out of memory in the shim itself */
#define VE_ERR_NOBASIS 3     /* lrs_getfirstbasis failed (e.g. infeasible system) */
#define VE_ERR_BADINPUT 4    /* bad dimensions or NULL pointer */

VE_EXPORT const char *
ve_lrs_version (void)
{
  return TITLE VERSION;
}

/* Enumerate the vertices (and extreme rays) of { x : Ax <= b }.
 *
 * Inputs are exact rationals: entry (i,j) of A is A_num[i*n+j] / A_den[i*n+j]
 * (row-major), entry i of b is b_num[i] / b_den[i].  All denominators must be
 * nonzero and positive.
 *
 * On success (*out_data) points to a malloc'd row-major (*out_rows) x (n+1)
 * matrix of doubles.  In each row, column 0 is 1.0 for a vertex and 0.0 for an
 * extreme ray; columns 1..n are the coordinates.  The caller must release the
 * matrix with ve_free().  On failure nothing has to be released.
 */
VE_EXPORT int
ve_enumerate (long m, long n,
              const long *A_num, const long *A_den,
              const long *b_num, const long *b_den,
              long *out_rows, double **out_data)
{
  lrs_dic *P = NULL;      /* structure for holding current dictionary and indices */
  lrs_dat *Q = NULL;      /* structure for holding static problem data            */
  lrs_mp_matrix Lin;      /* holds input linearities if any are found             */
  lrs_mp_matrix matrix = NULL;   /* accumulated output rows                       */
  lrs_mp_matrix tmpMatrix;
  long col;               /* output column index for dictionary                   */
  long prune = FALSE;     /* if TRUE, getnextbasis will prune tree and backtrack  */
  long *num = NULL;
  long *den = NULL;
  long i, j;
  long numVertices = 0;
  long allocated = 8;
  long oldAllocated;
  double *result = NULL;
  double tmpDoub;
  int status = VE_OK;
  FILE *devnull = NULL;

  if (m < 1 || n < 1 || A_num == NULL || A_den == NULL || b_num == NULL
      || b_den == NULL || out_rows == NULL || out_data == NULL)
    return VE_ERR_BADINPUT;
  for (i = 0; i < m; i++)
    {
      if (b_den[i] == 0)
        return VE_ERR_BADINPUT;
      for (j = 0; j < n; j++)
        if (A_den[i * n + j] == 0)
          return VE_ERR_BADINPUT;
    }

  *out_rows = 0;
  *out_data = NULL;

  if (!lrs_init ("vertexenum"))
    return VE_ERR_INIT;

  /* lrs_init leaves the global lrs_ofp/lrs_ifp NULL; point them at the null
     device so any stray fprintf inside lrslib cannot dereference NULL */
  devnull = fopen (VE_NULL_DEVICE, "w");
  if (devnull != NULL)
    {
      lrs_ifp = devnull;
      lrs_ofp = devnull;
    }

  Q = lrs_alloc_dat ("LRS globals");
  if (Q == NULL)
    {
      status = VE_ERR_INIT;
      goto done;
    }
  Q->m = m;
  Q->n = n + 1;

  P = lrs_alloc_dic (Q);
  if (P == NULL)
    {
      status = VE_ERR_INIT;
      goto done;
    }

  /* copy the input: lrs H-representation row i is  b_i + (-A_i) x >= 0 */
  num = (long *) malloc ((n + 1) * sizeof (long));
  den = (long *) malloc ((n + 1) * sizeof (long));
  if (num == NULL || den == NULL)
    {
      status = VE_ERR_NOMEM;
      goto done;
    }
  for (i = 0; i < m; i++)
    {
      num[0] = b_num[i];
      den[0] = b_den[i];
      for (j = 0; j < n; j++)
        {
          num[j + 1] = -A_num[i * n + j];
          den[j + 1] = A_den[i * n + j];
        }
      lrs_set_row (P, Q, i + 1, num, den, GE);
    }

  if (!lrs_getfirstbasis (&P, Q, &Lin, TRUE))
    {
      status = VE_ERR_NOBASIS;
      goto done;
    }

  matrix = lrs_alloc_mp_matrix (allocated, Q->n);
  if (matrix == NULL)
    {
      status = VE_ERR_NOMEM;
      goto done;
    }

  /* reverse search: collect one output row per vertex/ray */
  do
    {
      for (col = 0; col <= P->d; col++)
        {
          if (lrs_getsolution (P, Q, matrix[numVertices], col))
            {
              numVertices++;
              if (numVertices == allocated)
                {
                  oldAllocated = allocated;
                  allocated *= 2;
                  tmpMatrix = matrix;
                  matrix = lrs_alloc_mp_matrix (allocated, Q->n);
                  if (matrix == NULL)
                    {
                      matrix = tmpMatrix;
                      allocated = oldAllocated;
                      status = VE_ERR_NOMEM;
                      goto done;
                    }
                  for (i = 0; i < numVertices; i++)
                    for (j = 0; j < Q->n; j++)
                      copy (matrix[i][j], tmpMatrix[i][j]);
                  lrs_clear_mp_matrix (tmpMatrix, oldAllocated, Q->n);
                }
            }
        }
    }
  while (!Q->lponly && lrs_getnextbasis (&P, Q, prune));

  /* convert to doubles: row i of the result is
     [is_vertex, x_1, ..., x_n] with is_vertex 1.0 for vertices, 0.0 for rays */
  if (numVertices > 0)
    {
      result = (double *) malloc ((size_t) numVertices * Q->n * sizeof (double));
      if (result == NULL)
        {
          status = VE_ERR_NOMEM;
          goto done;
        }
      for (i = 0; i < numVertices; i++)
        {
          for (j = 0; j < Q->n; j++)
            {
              if (zero (matrix[i][0]))
                mptodouble (matrix[i][j], &tmpDoub);   /* extreme ray */
              else
                rattodouble (matrix[i][j], matrix[i][0], &tmpDoub);  /* vertex: matrix[i][0] is the denominator */
              result[i * Q->n + j] = tmpDoub;
            }
        }
    }

  *out_rows = numVertices;
  *out_data = result;

done:
  free (num);
  free (den);
  if (matrix != NULL)
    lrs_clear_mp_matrix (matrix, allocated, Q->n);
  if (P != NULL)
    lrs_free_dic (P, Q);
  if (Q != NULL)
    lrs_free_dat (Q);
  if (devnull != NULL)
    {
      lrs_ifp = NULL;
      lrs_ofp = NULL;
      fclose (devnull);
    }
  return status;
}

VE_EXPORT void
ve_free (double *data)
{
  free (data);
}
