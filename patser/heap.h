/* Copyright 1999, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


#include <stddef.h>

/* Functions for manipulating priority queues with a heap.
 * HEAP CONTAINS THE ITEMS WITH THE GREATEST VALUES;
 * THEREFORE, THE ELEMENT[1] CONTAINS THE MINIMUM VALUE IN THE HEAP.
 * The [0] first item in the heap array is not used, but must contain a
 * sentinel value (-INFINITY) if the function "heap_up()" is to be used
 * [functions calling "heap_up()": heap_insert(), tieheap_add()]. */

#if 0
/* Requires a comparing function.  The return value of the comparing
   function is the negative of the value used by qsort() and
   bsearch(); thus, the comparing fuction can be used by both heap_sort()
   and bsearch() even though heap_sort() sorts in descending order. */

/* The comparing function.  The input is pointers to the elements.
 * Returns -1 if elt2 < elt1;
 *          0 if elt2 = elt1;
 *          1 if elt2 > elt1. */
int compar(const void *elt1, const void *elt2);
#endif


/* Functions for manipulating heaps: */

/* Sort an unordered array in descending order using the heap sort.
 * First index of the array is 0.  Assume (*compar)() returns:
 * -1 if elt2 < elt1; 0 if elt2 = elt1; 1 if elt2 > elt1. */
extern void hsort(
     void *array,            /* The unordered array. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Order an unordered array into a heap. */
extern void heap_initialize(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Sort a heap array in descending order.
 * That is, turn a heap array into an ordered array with the top
 * of the priority queue at the bottom of the list. */
extern void heap_sort(
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Replace first item in the heap, if "new_elt" is greater than first element.
 * "new_elt" points to whichever element is less and, thus, off the heap. */
extern void heap_replace(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Assign current element to appropriate position in heap by moving down. */
extern void heap_down(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Insert an element at the end of the heap.
 * Returns the new count of the number of elements in the heap.
 * This function does not make sure the array is sufficiently large.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
extern int heap_insert(
     void *new_elt,          /* The new element to be inserted. */
     void *heap,             /* The array containing the heap. */
     size_t count,           /* The number of elements in the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Assign current element to appropriate position in the heap by moving up.
 * The [0] item of the heap array must contain a sentinel value (-INFINITY). */
extern void heap_up(
     int position,           /* The current position in the heap. */
     void *heap,             /* The array containing the heap. */
     size_t size,            /* The size of one of the elements of the heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/


/* Functions for manipulating tie heaps: a heap containing the
   second lowest score, followed by a list of the lowest scores. */

/* Order an unordered array into a tie heap. */
extern void tieheap_initialize(
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
extern void *tieheap_add(
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
extern void heap_to_tieheap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t *heap_count,     /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/

/* Convert a tie heap into a standard heap. */
extern void tieheap_to_heap(
     void *tieheap,          /* The array containing the tie heap. */
     size_t count,           /* The number of elements in the tie heap. */
     size_t heap_count,      /* The number of elements in the heap portion of
				tie heap (excludes low scores at the end). */
     size_t count_target,    /* Target for number of elements in tie heap. */
     size_t size,            /* Size of one of the elements of the tie heap. */
     int (*compar)(const void *elt1, const void *elt2));/* Comparing function*/
