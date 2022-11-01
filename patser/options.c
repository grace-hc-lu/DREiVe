/* Copyright 1990, 1994, 1995, 1996, 2000, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


/* Declaration of global variables:
  * Input files.
  * Limiting the scores to print.
  * Alphabet options.
  * File pointers.
  * Matrix file.
  * Determined by file containing the list of sequence files.
  * Determined by the contents of the sequence files.
  */


#include <stdio.h>
#include "definitions.h"     /* File containing the "#define" statements. */


/* Testing flag indicating calculations of cutoff scores are being tested. */
char Test = 0;        /* 0: no testing.
		       * 1: determine cutoff score analytically.
		       * 2: determine cutoff score analtyically and exactly. */


/* Input files. */
char *Mat_file = "matrix";                 /* The name of the matrix file. */
char *Seq_file = (char *)NULL;             /* The name of the sequence file. */


/* Adjusting or restricting which scores are printed. */
char Thresh = '\0';    /* Options for having threshold scores.
			* = '\0': no thresholds;
			* = 'l': lower threshold;
			* = 'u': upper threshold;
			* = 'b': both lower and upper thresholds. */
double L_thresh;       /* Lower threshold for printing scores, inclusive. */
double U_thresh;       /* Upper threshold for printing scores, exclusive. */
int Top = 0;	  /* Print indicated number of top scores for each sequence. */
char Print_order = 0; /* Flag determining how the top scores are printed.
		       *  0: print in order of increasing position number.
		       *  1: print in order of decreasing score. */
char Error = 1;       /* Flag determining how unrecognized symbols are treated.
		       *  0: treat as an error.
		       *  1: treat as a discontinuity, but print a warning.
		       *  2: treat as a discontinuity, and print NO warning. */
double Max_ln_p_value = -INFINITY; /* Upper threshold for printing p-values. */
char Auto_cutoff = NO;/* YES: determine a cutoff score assuming the
		       * ln(cutoff p-value) = -adjusted_information_content. */
double Equal_err = ERROR;/* Small difference for considering 2 scores equal. */


#if 0
/* Alphabet options that are defined in "alpha.c". */
char *Alpha_file = "alphabet";             /* The name of the alphabet file. */
char Comp_status = NO; /* Flag for whether to scan complementary strand. */
char Ascii = YES;      /* Flag for whether alphabet characters are ASCII or
			* integer.  YES: ASCII; NO: integer. */
char Case_sensitive = 1; /* Flag for whether letter alphabet is case sensitive.
			  * 0: case sensitive; 1: not case sensitive;
			  * 2: not case sensitive, but indicate the locations
			  *    of lowercase letters. */
char Comp_flag = 0;    /* Flag indicating whether alphabet is complementary and
			* whether any letters are their own complement.
			*     0: no complements; 1: no own; 2: some own */
int A_size = 0;  /* Alphabet size---i.e., number of letters. */
int *A;          /* Translation between an integral index and the corresponding
		  * letter/integer.  Convention for nucleic acid alphabets:
		  * if "n" is the index of a letter, its complement will have
		  * the index (n+1), if the letters are not the same. */
int *A_comp;     /* Translation between the index of a letter and
		  * the index of its complement.*/
double *P;       /* The prior frequencies of the letters: obtained from the
		  * normalization data. */
#endif


/* File pointers. */
FILE *M_fp;        /* Pointer to file containing the matrix of integers. */
FILE *S_fp;        /* Pointer to file containing the list of sequence files. */


/* Matrix file. */
char Sum_mat = YES;  /* Flag indicating whether matrix file contains an
		      * alignment matrix (YES) or a weight matrix (NO). */
double Fudge = 1.0;  /* Correction added to the elements of the alignment
		      * matrix before converting to a weight matrix. */
char Vertical = NO;  /* Flag indicating whether matrix is vertical (YES)---i.e.
		      * rows correspond to positions---or horizontal (NO). */
double **Matrix;     /* The matrix after converting to a weight matrix. */
double **Align_mat;  /* The alignment matrix. */
int L = 0;           /* The width of the matrix. */
int L_1;             /* The width of the matrix minus 1. */
char *Letter_array;  /* Flags for making sure each letter occurs only
		      * once in the matrix file. */
int *Letter_index;   /* The letter index of the corresponding row or
		      * column of the matrix. */
char Print = 0;    /* Flag indicating whether to print the specificity matirx.
		     * 0: do not print; 1: print a horizontal matrix;
		     * 2: print a vertical matrix. */

/* Information contents. */
double Inf;         /* Information content of the alignment matrix. */
double Inf_adj;     /* Adjusted information content of the alignment matrix. */
double Inf_fudge;   /* Information content of fudged alignment matrix. */


/* Determined by file containing the list of sequence files. */
char File[BIG_CHUNK_SIZE]; /* The name of the current sequence file. */
FILE *Fp;                  /* Pointer to the current sequence file. */
char Circle;    /* Flag indicating whether the current sequence is circular. */


/* Determined by the contents of a sequence file. */
LETTER *Sequence = (LETTER *)NULL;  /* The sequence of the current sequence. */
int Seq_length;                     /* The length of the current sequence. */


/* Variables for determining the array of marginal probabilities
 * for integer scores. */
double *Marginal_prob = (double *)NULL;    /* ln(p-value of integer scores). */
double Max_int_range = 1.0e4;   /* The range, starting from 0, for
				 * approximating a score with an integer. */
double Min_score = 0.0;         /* Minimum score for calculating
				 * marginal probabilities. */
double Multiple = 1.0;/* For multiplying scores before converting to integer.*/
double Max_exact_score;   /* Maximum score of the exact weight matrix. */
double Min_exact_score;   /* Minimum score of the exact weight matrix. */
double Min_pseudo_score;  /* Minimum score of weight matrix after
			   * assigning values to the -INFINITY elements. */
char Minus_inf = NO;   /* Flag indicating whether matrix contains -INFINITY. */
int Max_int_score;        /* Maximum score of integral weight matrix. */
int Min_int_score = 0;    /* Minimum ingegral score - 1 for
			   * estimating probabilities.*/
double *ArrayA;          /* Holds final ln(probabilities) of integer scores. */
double *Ln_P;             /* ln(P[i] */
int Skip = 1;  /* Number of positions to skip in the marginal probability array
		* when calculating marginal probabilities analytically. */
