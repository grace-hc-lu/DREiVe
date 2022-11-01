/* Copyright 1990, 1995, 1996, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


/* This file contains the "#define" statements. */

/* #define TEST 1             When defined, the functions for determining
			    * exact and analytical cutoffs are included. */
#define HIGH_CUTOFF 1      /* When defined, the cutoff score is such that the
			    * marginal probability is less than or equal to
			    * the target.  When NOT defined, the cutoff score
			    * is such that the marginal probability is
			    * greater than or equal to the target. */
#define NEW_FUDGE 1        /* When defined, the small sample size correction
			    * is Fudge * a priori probability.  When NOT
			    * defined, the small sample size correction is
			    * simply Fudge. */


#define HEAP_CHUNK_SIZE 32  /* Size unit for allocating memory for tie heap. */
#define BIG_CHUNK_SIZE 512   /* The maximum size for file names. */
#define CHUNK_SIZE 512       /* The size unit for allocating dynamic memory. */
#define NAME_LENGTH 20       /* Presumed maximum length of a sequence name. */

#define BASE 2.71828182846   /* The base for taking logarithms. */
#define INFINITY 1.0e100     /* A large double. */
#define INT_INF 2147483647   /* The largest 4 byte positive integer. */
#define ERROR 1.0e-6         /* Leeway factor when determining whether
			      * two scores or two column sums are the same. */
#define PI 3.14159265359
#define SQRT2PI 2.50662827463 /* (2 * PI)^{0.5} */

#define LETTER char           /* Storage type of the letters of a sequence. */

#define YES 1
#define NO 0


/* Round a floating point number to the closest integer.
 * Positive numbers: >= 0.5 is rounded up.
 * Negative numbers: <= -0.5 is rounded down. */
#define ROUND_TO_INT(n) (((n) > 0.0) ? ((int)((n) + 0.5)) : ((int)((n) - 0.5)))

/* Sum numbers given their logarithms.  Subtract the larger logarithm
 * from each logarithm to avoid adding 0's or infinities within the
 * precision of the machine.  Returns the logarithm of the sum. */
#define SUM_LN(n1, n2) (((n1) >= (n2)) ? \
			((n1) + log(exp((n2) - (n1)) + 1.0)) : \
			((n2) + log(exp((n1) - (n2)) + 1.0)))

/* Subtract the second number from the first given their logarithms.
* Subtracts the logarithm of the first number from each logarithm to
* avoid 0's or infinities within the precision of the machine.
* Returns the logarithm of the difference.
* If the second number is >= the first number, returns -INFINITY. */
#define DIFF_LN(n1, n2) (((n1) >= (n2) + ERROR) ? \
			 ((n1) + log(1.0 - exp((n2) - (n1)))) : \
			 (((n1) > (n2)) ? ((n1) + log((n1) - (n2))) : \
			  -INFINITY))


/* Structure for holding information about top scores when "Top" > 0. */
typedef struct struct_PAIR
{
  int             position; /* The position within the current sequence.
			     * Negative integers means complementary strand. */
  double          score;    /* The score of the L-mer. */
} PAIR;

/* Structure for holding the heap when the exact cutoff score is determined.
 *      1) the score of the corresponding L-mer; and
 *      2) the probability of the corresponding L-mer. */
struct struct_HEAP
{
  double        score;
  double        ln_prob;
};
typedef struct struct_HEAP HEAP;
