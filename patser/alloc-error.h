/* Copyright 1997 Gerald Z. Hertz
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


/* Declare housekeeping functions defined in "alloc_error.c" and
 * "alloc-error-debug.c". */


/* Malloc with an error message and an exit if unsuccessful. */
extern void *malloc_error(
     unsigned size,       /* Size of the array. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */

/* Calloc with an error message and an exit if unsuccessful.
 * Also, the element count is an "int" rather than an "unsigned". */
extern void *calloc_error(
     int elt_count,       /* Number of elements in the array. */
     unsigned elt_size,   /* Size of each element. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */

/* Realloc with an error message and an exit if unsuccessful.  Also,
 * checks whether "array" is equal to NUll, if so, uses "malloc()" to
 * avoid problems with compilers that cannot realloc NULL pointers. */
extern void *realloc_error(
     void *array,         /* The array whose size is to be modified. */
     unsigned size,       /* Size of the array. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */

/* Realloc an array with an error message and an exit if unsuccessful.
 * Also, checks whether "array" is equal to NUll, if so, uses "malloc()"
 * to avoid problems with compilers that cannot realloc NULL pointers. */
extern void *recalloc_error(
     void *array,         /* The array whose size is to be modified. */
     int elt_count,       /* Number of elements in the array. */
     unsigned elt_size,   /* Size of each element. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */

/* Free dynamic memory, but do nothing if handed a NULL pointer. */
extern void free_error(
     void *array,         /* The array being freed. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */

/* A general function for reporting program bugs. */
extern void bug_report(
     char *function);     /* Function containing the bug. */



/* The following function is only defined in "alloc-error-debug.c".
 * In "alloc-error.c", the function is defined to do nothing. */

/* Check an array to make sure the magic ints and char are all consistent
 * with unmolested memory.  Returns a void pointer to the gross memory.
 * When called outside of "alloc-error-debug.c",  alloc_function = "check". */
extern void *check_array(
     char *alloc_function,/* Name of allocating function.*/
     void *array,         /* Array being checked. */
     char *variable,      /* Name of the variable. */
     char *function);     /* Name of the calling function. */
