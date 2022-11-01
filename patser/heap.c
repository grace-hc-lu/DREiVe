/* Copyright 1997, 1999, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */

#ifndef OPTIONS
#include <string.h>
#include "alloc-error.h"
#define HEAP_CHUNK_SIZE 32  /* Size unit for allocating memory for tie heap. */
#else
#include "options.h"
#endif
#include "heap.h"


/* Functions for manipulating priority queues with a heap.
 * HEAP CONTAINS THE ITEMS WITH THE GREATEST VALUES;
 * THEREFORE, THE ELEMENT[1] CONTAINS THE MINIMUM VALUE IN THE HEAP.
 * The [0] first item in the heap array is not used, but must contain a
 * sentinel value (-INFINITY) if the function "heap_up()" is to be used
 * [functions calling "heap_up()":
 *  heap_insert(), tieheap_add(), tieheap_to_heap()]. */

/* Also contains functions for manipulating "tie heaps", my invention
   storing additional elements when all the elements tied for the
   lowest score are also saved.  The priority queue portion contains
   all the elements down to the second lowest score, which is at index
   position [1].  However, the elements having the lowest score are at
   the end of the heap array. */

#if 0
/* Requires a comparing function.  The return value of the comparing
   function is the negative of the value used by qsort() and
   bsearch(); thus, the comparing fuction will achieve the same result
   with hsort() and qsort(), even though hsort() sorts in descending order. */

/* The comparing function.  The input is pointers to the elements.
 * Returns -1 if elt2 < elt1;
 *          0 if elt2 = elt1;
 *          1 if elt2 > elt1. */
int compar(const void *elt1, const void *elt2);


/* Functions for manipulating heaps: */

/* Sort an unordered array in descending order using the heap sort.
 * First index of the array is 0.  Assume (*compar)() returns:
 * -1 if elt2 < elt1; 0 if elt2 = elt1; 1 if elt2 > elt1. */
void hsort(
     void *array,            /* The unordered array. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Order an unordered array into a heap. */
void heap_initialize(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Sort a heap array in descending order.
 * That is, turn a heap array into an ordered array with the top
 * of the priority queue at the bottom of the list. */
void heap_sort(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Replace first item in the heap, if "new_elt" is greater than first element.
 * "new_elt" points to whichever element is less and, thus, off the heap. */
void heap_replace(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Assign current element to appropriate position in heap by moving down.
 * Assumes that everthing below the starting position is a proper heap. */
void heap_down(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Insert an element at the end of the heap.
 * Returns the new count of the number of elements in the heap.
 * This function does not make sure the array is sufficiently large.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
int heap_insert(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Assign current element to appropriate position in the heap by moving up.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
void heap_up(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/


/* Functions for manipulating tie heaps: a heap containing the
   second lowest score, followed by a list of the lowest scores. */

/* Order an unordered array into a tie heap. */
void tieheap_initialize(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end).
				Its value is determined by this function. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Add a new element to the tie heap if its score is greater than or
   equal to the current lowest score.  Assumes the size of the tie
   heap (count) is >= its target size (count_target).  Delete elements
   and readjust the tie heap as necessary.  Returns the memory address
   of the heap in case additional memory has be reallocated.

   If "new_elt" replaces an existing element, its final value will
   equal the element being replaced rather than its original value.
   Thus, deleted values will be "new_elt" and "tiheap[new count + 1]"
   through "tieheap[old count]". */
void *tieheap_add(
     void *new_elt,          /* The new element to be inserted. */
     void *tieheap,          /* The array containing the tie heap. */
     size_t *count,          /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t *count_alloc,    /* Number of elements allocated for tie heap. */
     size_t count_target,    /* Target for number of elements in tie heap. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Convert a standard heap into a tie heap. */
void heap_to_tieheap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Convert a tie heap into a standard heap. */
void tieheap_to_heap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t heap_count,      /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t count_target,    /* Target for number of elements in tie heap. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/
#endif



/* Sort an unordered array in descending order using the heap sort.
 * First index of the array is 0.  Assume (*compar)() returns:
 * -1 if elt2 < elt1; 0 if elt2 = elt1; 1 if elt2 > elt1. */
void hsort(
     void *array,            /* The unordered array. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  /* Create a heap array with an intial index of 1. */
  void *heap = (char *)array - size;

  heap_initialize(heap, count, size, compar);
  heap_sort(heap, count, size, compar);
}

/* Order an unordered array into a heap. */
void heap_initialize(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  int position;               /* Position within the heap. */

  for (position = count / 2; position >= 1; --position)
    heap_down(position, heap, count, size, compar);
}

/* Sort a heap array in descending order.
 * That is, turn a heap array into an ordered array with the top
 * of the priority queue at the bottom of the list. */
void heap_sort(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *elt_1;                         /* Copy of first element of heap. */
  void *heap_1 = (char *)heap + size;  /* Pointer to first element of heap. */
  void *heap_c;                        /* Pointer to last element of heap. */


  /* Allocate memory for "elt_1". */
  elt_1 = malloc_error(size, "elt_1", "heap_sort()");

  while (count > 1)
    {
      heap_c = (char *)heap + count * size;

      memcpy(elt_1, heap_1, size);
      memcpy(heap_1, heap_c, size);
      memcpy(heap_c, elt_1, size);

      heap_down(1, heap, --count, size, compar);
    }

  /* Free memory for "elt_1". */
  free_error(elt_1, "elt_1", "heap_sort()");
}

/* Replace first item in the heap, if "new_elt" is greater than first element.
 * "new_elt" points to whichever element is less and, thus, off the heap. */
void heap_replace(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *elt_1;                         /* Copy of first element of heap. */
  void *heap_1 = (char *)heap + size;  /* Pointer to first element of heap. */


  if ((*compar)(heap_1, new_elt) > 0)
    {
      /* Allocate memory for "elt_1". */
      elt_1 = malloc_error(size, "elt_1", "heap_replace()");

      memcpy(elt_1, heap_1, size);
      memcpy(heap_1, new_elt, size);
      memcpy(new_elt, elt_1, size);
      heap_down(1, heap, count, size, compar);

      /* Free memory for "elt_1". */
      free_error(elt_1, "elt_1", "heap_replace()");
    }
}

/* Assign current element to appropriate position in heap by moving down.
 * Assumes that everthing below the starting position is a proper heap. */
void heap_down(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *elt;                    /* The element to be placed. */
  int child_idx = position * 2; /* Child index of current position in heap. */
  void *parent;                 /* Pointer to the current position. */
  void *child_1;                /* Pointer to the younger or smaller child . */
  void *child_2;                /* Pointer to the older child. */



  /* Determine current "parent".
   * Allocate memory for "elt" and determine its value. */
  parent = (char *)heap + position * size;
  elt = malloc_error(size, "elt", "heap_down()");
  memcpy(elt, parent, size);

  /* Work down the heap to find the appropriate position in the heap. */
  while (child_idx < count)
    {
      /* Determine pointers to the children. */
      child_1 = (char *)heap + child_idx * size;
      child_2 = (char *)child_1 + size;

      /* Determine which child has the smaller value. */
      if ((*compar)(child_1, child_2) < 0)
	{
	  ++child_idx;
	  child_1 = child_2;
	}

      /* Determine whether current position should be
       * switched with its smallest child. */
      if ((*compar)(elt, child_1) < 0)
	{
	  /* Replace the parent with the child. */
	  memcpy(parent, child_1, size);

	  /* Precede down to the next level. */
	  child_idx *= 2;
	  parent = child_1;
	}
      else break;
    }

  /* Check whether the child equals count. */
  if (child_idx == count)
    {
      /* Determine pointer to the child. */
      child_1 = (char *)heap + child_idx * size;

      /* Determine whether current position should be
       * switched with its child. */
      if ((*compar)(elt, child_1) < 0)
	{
	  /* Replace the parent with the child. */
	  memcpy(parent, child_1, size);
	  parent = child_1;
	}
    }

  /* Place the element back into the heap. */
  memcpy(parent, elt, size);
  free_error(elt, "elt", "heap_down()");
}

/* Insert an element at the end of the heap.
 * Returns the new count of the number of elements in the heap.
 * This function does not make sure the array is sufficiently large.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
int heap_insert(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  /* Copy the "new_elt" to the end of the heap. */
  memcpy((char *)heap + (++count) * size, new_elt, size);

  /* Assign the "new_elt" to its proper position within the heap. */
  heap_up(count, heap, size, compar);

  /* Return the new number of elements in the heap. */
  return(count);
}

/* Assign current element to appropriate position in the heap by moving up.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
void heap_up(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *elt;                      /* The element to be placed. */
  int parent_idx;                 /* Index to parent of current position. */
  void *parent;                   /* Pointer to parent of current position. */
  void *position_pointer;         /* Pointer to the current position. */


  /* Initialize the "parent_idx", "position_pointer", and "parent" pointer. */
  parent_idx = position / 2;
  position_pointer = (char *)heap + position * size;
  parent = (char *)heap + parent_idx * size;

  /* Allocate memory for "elt" and assign value of the element to be placed. */
  elt = malloc_error(size, "elt", "heap_up()");
  memcpy(elt, position_pointer, size);

  /* Work up the heap to find the appropriate position in the heap. */
  while ((*compar)(parent, elt) < 0)
    {
      /* Copy the parent to the current position. */
      memcpy(position_pointer, parent, size);

      /* Precede up to the next level of the heap. */
      position_pointer = parent;
      parent_idx = parent_idx / 2;
      parent = (char *)heap + parent_idx * size;
    }

  /* Place the element back into the heap. */
  memcpy(position_pointer, elt, size);
  free_error(elt, "elt", "heap_up()");
}


/* Functions for manipulating tie heaps: a heap containing the
   second lowest score, followed by a list of the lowest scores. */

/* Order an unordered array into a tie heap. */
void tieheap_initialize(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end).
				Its value is determined by this function. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  if (count == 1) *heap_count = 0;

  else
    {
      /* Order the tie heap into a heap. */
      heap_initialize(tieheap, count, size, compar);

      /* Create a tie heap from a standard heap. */
      *heap_count = count;
      heap_to_tieheap(tieheap, count, heap_count, size, compar);
    }
}

/* Add a new element to the tie heap if its score is greater than or
   equal to the current lowest score.  Assumes the size of the tie
   heap (count) is >= its target size (count_target).  Delete elements
   and readjust the tie heap as necessary.  Returns the memory address
   of the heap in case additional memory has be reallocated.

   If "new_elt" replaces an existing element, its final value will
   equal the element being replaced rather than its original value.
   Thus, deleted values will be "new_elt" and "tiheap[new count + 1]"
   through "tieheap[old count]". */
void *tieheap_add(
     void *new_elt,          /* The new element to be inserted. */
     void *tieheap,          /* The array containing the tie heap. */
     size_t *count,          /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t *count_alloc,    /* Number of elements allocated for tie heap. */
     size_t count_target,    /* Target for number of elements in tie heap. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *tmp_elt; /* Temporary location for lowest element being overwritten. */

  /* Pointer to last element of the heap portion of the tie heap. */
  void *heap_end;

  /* -1 if new_elt < end of tieheap;
      0 if new_elt = end of tie heap;
      1 if new_elt > end of tie heap */
  int status = (*compar)((char *)tieheap + *count * size, new_elt);



  /* If the new element is tied with the current lowest score, insert it
     at the end of the tie heap after ensuring there is enough memory. */
  if (status == 0)
    {
      /* Increment size of tie heap; make sure enough memory is allocated. */
      if (++(*count) == *count_alloc)
	tieheap = recalloc_error(tieheap, *count_alloc += HEAP_CHUNK_SIZE,
				 size, "tieheap", "tieheap_add() 1");

      /* Copy the "new_elt" to the end of the tieheap. */
      memcpy((char *)tieheap + *count * size, new_elt, size);
    }


  /* If new element is greater than current lowest score, insert into heap. */
  else if (status == 1)
    {
      /* Increment the size of the standard heap;
	 determine if the current lowest scores are not being dropped. */
      if (++(*heap_count) < count_target)
	{
	  /* Increment tie-heap size; make sure enough memory is allocated. */
	  if (++(*count) == *count_alloc)
	    tieheap = recalloc_error(tieheap, *count_alloc += HEAP_CHUNK_SIZE,
				     size, "tieheap", "tieheap_add() 2");

	  /* Shift lowest score from the end of heap to end of tie heap. */
	  heap_end = (char *)tieheap + *heap_count * size;
	  memcpy((char *)tieheap + *count * size, heap_end, size);

	  /* Copy the new element to the end of the standard heap. */
	  memcpy(heap_end, new_elt, size);

	  /* Update the standard heap. */
	  heap_up(*heap_count, tieheap, size, compar);
	}

      /* Drop the current lowest scores. */
      else if (*heap_count == count_target)
	{
	  /* Isolate the end of the standard heap.
	     Drop low scores (except the first) by shortening count. */
	  heap_end = (char *)tieheap + *heap_count * size;
	  *count = count_target;

	  /* Copy the new element to the end of the standard heap.
	   * Copy the old lowest element that is being overwritten by
	   * new_elt to new_elt. */
	  tmp_elt = malloc_error(size, "tmp_elt", "tieheap_add()");
	  memcpy(tmp_elt, heap_end, size);
	  memcpy(heap_end, new_elt, size);
	  memcpy(new_elt, tmp_elt, size);
	  free_error(tmp_elt, "tmp_elt", "tieheap_add()");

	  /* Update the standard heap. */
	  heap_up(count_target, tieheap, size, compar);

	  /* Create a tie heap from the standard heap. */
	  heap_to_tieheap(tieheap, *count, heap_count, size, compar);
	}

      else bug_report("tieheap_add()");
    }

  return(tieheap);
}

/* Convert a standard heap into a tie heap.  Can also be used on a
   partial tie heap in which some members of the standard heap belong
   in the terminal array. */
void heap_to_tieheap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  void *elt_1;                         /* Copy of the first element of heap. */
  void *heap_1;             /* Pointer to first element of heap. */
  void *heap_end;           /* Pointer to last element of the standard heap. */
  void *tieheap_end;        /* Pointer to last element of the tie heap. */


  if ((count <= 1) || (*heap_count == 0)) *heap_count = 0;

  else
    {
      /* Isolate the first and last elements of the tie heap. */
      heap_1 = (char *)tieheap + size;
      tieheap_end = (char *)tieheap + size * count;
      
      if ((*heap_count == count) || ((*compar)(heap_1, tieheap_end) == 0))
	{
	  /* Allocate memory for "elt_1". */
	  elt_1 = malloc_error(size, "elt_1", "heap_to_tieheap()");

	  do
	    {
	      if (*heap_count == 1)
		{
		  *heap_count = 0;
		  break;
		}
	      else
		{
		  /* Exchange tieheap[1] with last element of standard heap. */
		  heap_end = (char *)tieheap + *heap_count * size;
		  memcpy(elt_1, heap_1, size);
		  memcpy(heap_1, heap_end, size);
		  memcpy(heap_end, elt_1, size);

		  /* Shorten and reinitialize standard heap portion. */
		  --(*heap_count);
		  heap_down(1, tieheap, *heap_count, size, compar);
		}
	    }
	  while ((*compar)(heap_1, tieheap_end) == 0);

	  /* Free memory for "elt_1". */
	  free_error(elt_1, "elt_1", "heap_to_tieheap()");
	}
    }
}

/* Convert a tie heap into a standard heap. */
void tieheap_to_heap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t heap_count,      /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t count_target,    /* Target for number of elements in tie heap. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2)) /* Comparing function*/
{
  int i;


  /* Order an unordered array into a heap. */
  if (count < count_target) heap_initialize(tieheap, count, size, compar);

  /* Convert a tie heap into a standard heap. */
  else
    {
      for (i = heap_count + 1; i <= count; ++i)
	heap_up(i, tieheap, size, compar);
    }
}
