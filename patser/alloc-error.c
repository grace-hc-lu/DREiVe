/* Copyright 1990, 1992, 1993, 1994, 1995, 1996, 1997 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   Dept. of Molecular, Cellular, and Developmental Biology
 *   University of Colorado
 *   Campus Box 347
 *   Boulder, CO  80309-0347
 *
 *   hertz@colorado.edu
 */


#ifndef VMS
#include <stdio.h>
#include <stdlib.h>
#else
#include stdio
#include stdlib
#endif


/* Gerald Z. Hertz
 * Source for the following seven functions:
 *     malloc_error()
 *     calloc_error()
 *     realloc_error()
 *     recalloc_error()
 *     free_error()
 *
 *     bug_report()
 *     check_array()
 *
 * The first 3 functions are simply "malloc()", "calloc()", and "realloc()"
 * with an error message.  In addition, "realloc_error()" uses "malloc()"
 * if the pointer is NUll in order to ensure portability.
 *
 * recalloc_error() combines calloc() with realloc().
 * free_error() is free() with an error message.
 *
 * bug_report() is a general function for reporting program bugs.
 * check_array() makes sure memory is not molested, but
 *               is only defined in "alloc-error-debug.c".
 */


/* WARNING: "calloc_error()" WILL NOT  INITIALIZE
 * ARRAYS TO 0'S SINCE IT REALLY USES "malloc()". */


/* Malloc with an error message and an exit if unsuccessful. */
void *malloc_error(
     unsigned size,    /* Size of the element. */
     char *variable,   /* Name of the variable. */
     char *function)   /* Name of the calling function. */
{
  void *array;


  array = (void *)malloc(size);
  if (array == (void *)NULL)
    {
      fprintf(stderr, "Cannot malloc space for \"%s\" in the ", variable);
      fprintf(stderr, "\"%s\" function.\n", function);
      exit(1);
    }
  else
    return(array);
}

/* Calloc with an error message and an exit if unsuccessful.
 * Also, the element count is an "int" rather than an "unsigned". */
void *calloc_error(
     int elt_count,        /* Number of elements in the array. */
     unsigned elt_size,    /* Size of each element. */
     char *variable,       /* Name of the variable. */
     char *function)       /* Name of the calling function. */
{
  void *array;


  array = (void *)malloc((unsigned)elt_count * elt_size);

  if (array == (void *)NULL)
    {
      fprintf(stderr, "Cannot calloc space for \"%s\" in the ", variable);
      fprintf(stderr, "\"%s\" function.\n", function);
      exit(1);
    }
  else
    return(array);
}

/* Realloc with an error message and an exit if unsuccessful.  Also,
 * checks whether "array" is equal to NUll, if so, uses "malloc()" to
 * avoid problems with compilers that cannot realloc NULL pointers. */
void *realloc_error(
     void *array,      /* The array whose size is to be modified. */
     unsigned size,    /* Size of the array. */
     char *variable,   /* Name of the variable. */
     char *function)   /* Name of the calling function. */
{
  if (array == (void *)NULL) array = (void *)malloc(size);
  else array = (void *)realloc(array, size);

  if (array == (void *)NULL)
    {
      fprintf(stderr, "Cannot realloc space for \"%s\" in the ", variable);
      fprintf(stderr, "\"%s\" function.\n", function);
      exit(1);
    }
  else
    return(array);
}

/* Realloc an array with an error message and an exit if unsuccessful.
 * Also, checks whether "array" is equal to NUll, if so, uses "malloc()"
 * to avoid problems with compilers that cannot realloc NULL pointers. */
void *recalloc_error(
     void *array,          /* The array whose size is to be modified. */
     int elt_count,        /* Number of elements in the array. */
     unsigned elt_size,    /* Size of each element. */
     char *variable,       /* Name of the variable. */
     char *function)       /* Name of the calling function. */
{
  if (array == (void *)NULL)
    array = (void *)malloc((unsigned)elt_count * elt_size);
  else array = (void *)realloc(array, (unsigned)elt_count * elt_size);

  if (array == (void *)NULL)
    {
      fprintf(stderr, "Cannot recalloc space for \"%s\" in the ", variable);
      fprintf(stderr, "\"%s\" function.\n", function);
      exit(1);
    }
  else return(array);
}

/* Free dynamic memory, but do nothing if handed a NULL pointer. */
void free_error(
     void *array,
     char *variable,   /* Name of the variable. */
     char *function)   /* Name of the calling function. */
{
  if (array != (void *)NULL) free(array);
}


/* A general function for reporting program bugs. */
void bug_report(
     char *function)
{
  fprintf(stderr, "PROGRAM BUG: %s\n", function);
  exit(1);
}


/* NULL function.
 * When called outside of "alloc-error-debug.c",  alloc_function = "check". */
void *check_array(
     char *alloc_function,  /* Name of the allocating function. */
     void *array,           /* Array being checked. */
     char *variable,        /* Name of the variable. */
     char *function)        /* Name of the calling function. */
{
  return((int *)NULL);
}
