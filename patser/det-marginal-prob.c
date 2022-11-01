/* Copyright 1996, 2000, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


#include "options.h"
/* #define PRINTDISTRIBUTION 1    Print final distribution when defined. */
#define PRINTRANGE 1           /* Print the range of scores and p-values. */
#define T_SIZE A_size      /* The total number of symbols in the matrix. */
#define LN_PROB_SUM 0.0    /* ln(sum of all possible probabilities)
			    * Can be > 0 in the presence of gaps. */


/* Functions for determining the exact p-values of the approximate matrix.
 * The following external variables are used:
 *
 *    T_SIZE                Multiple              Max_int_score
 *    P            	    Max_exact_score 	  Min_int_score
 *    Marginal_prob	    Min_exact_score 	  ArrayA
 *    Max_int_range	    Min_pseudo_score	  Ln_P
 *    Min_score    	    Minus_inf             Equal_err
 */


/* Properties of the integral weight matrix. */
static int *STmax_score;   /* Maximum score for each column. */
static int *STmin_score;   /* Minimum (non -INFINITY) score for each column. */
static int STmax_int_range;/* The range: Max_int_score - Min_int_score. */

static double *STarrayA = (double *)NULL;   /* Array for holding
					     * intermediate results. */
static double *STarrayB;   /* Array for holding current column. */
static double *STarrayC;   /* Array for holding intermediate results. */
static int *STarray_index; /* Index to STarrayB[],
			    * the array holding the current column. */


/* Determine the p-values of the approximate matrix.
 * Place p-values in Marginal_prob[];
 * Multiple is the conversion factor from exact scores to integral scores. */
void det_marginal_prob(weight_mat, width)
     double **weight_mat;           /* The weight matrix. */
     int width;                     /* Width of the matrix. */
{
  void determine_variables();
  void appr_mat_space();
  void determine_distribution();
  void determine_marginal_prob();


  /* Make sure the integer range has not been set too large. */
  if (Max_int_range > (double)(INT_INF - width - 3))
    {
      fprintf(stderr, "The integer range of %g is too large,\n",Max_int_range);
      fprintf(stderr, "it must be less than or equal to %d - %d - 3\n",
	      INT_INF, width);
      exit(1);
    }

  /* Reset Marginal_prob[0] to the beginning of the allocated memory. */
  Marginal_prob += Min_int_score;

  /* Allocate space for STmax_score[] and STmin_score[]; determine Minus_inf,
   * Max_exact_score, Min_exact_score, Min_pseudo_score, Multiple,
   * STmax_score[], STmin_score[], Max_int_score, Min_int_score,
   * and STmax_int_range. */
  determine_variables(weight_mat, width);

  /* Allocate space for "STarrayA[]", "STarrayB[]", "STarrayC[]",
   * and "STarray_index[]". */
  appr_mat_space();

  /* Determine the logarithm of the distribution of the integral scores of
   * the "weight_mat[][]" and place in ArrayA[]. */
  determine_distribution(weight_mat, width);

  /* Determine the logarithm of the p-values and place in Marginal_prob[]. */
  determine_marginal_prob();


  /* Free arrays except for STarrayA. */
  free_error(STmax_score, "STmax_score", "exact_appr_cutoff()");
  free_error(STmin_score, "STmin_score", "exact_appr_cutoff()");
  free_error(STarrayB, "STarrayB", "exact_appr_cutoff()");
  free_error(STarrayC, "STarrayC", "exact_appr_cutoff()");
  free_error(STarray_index, "STarray_index", "exact_appr_cutoff()");
}

/* Allocate space for STmax_score[] and STmin_score[]; determine Minus_inf,
 * Max_exact_score, Min_exact_score, Min_pseudo_score, Multiple, STmax_score[],
 * STmin_score[], Max_int_score, Min_int_score, and STmax_int_range. */
void determine_variables(weight_mat, width)
     double **weight_mat;     /* The weight matrix. */
     int width;               /* Width of the matrix. */
{
  int i, j;
  char max_min_equality;  /* YES: maximum and minimum (non -INFINITY)
			   *      scores are equal. */
  double *min_score;      /* Value of column element with lowest score. */
  double *max_score;      /* Value of column element with largest score. */
  int min_int_score;      /* Minimum integral score greater than -INFINITY. */
  double min_prob_score;  /* Minimum score for estimating probabilities. */
  double adj_minus_infinity(); /* Assign values to the -INFINITY elements.
				* Determine Min_pseudo_score. */


  /* Allocate space for "STmax_score[]" and "STmin_score[]". */
  STmax_score = (int *)calloc_error(width, sizeof(int),
				    "STmax_score", "determine_variables()");
  STmin_score = (int *)calloc_error(width, sizeof(int),
				    "STmin_score", "determine_variables()");
  max_score = (double *)calloc_error(width, sizeof(double),
				     "max_score", "determine_variables()");
  min_score = (double *)calloc_error(width, sizeof(double),
				     "min_score", "determine_variables()");

  /* Determine "Max_exact_score" and "Min_exact_score". */
  for (j = 0, Max_exact_score = Min_exact_score = 0.0; j < width; ++j)
    {
      for (i = 0, max_score[j] = -INFINITY, min_score[j] = INFINITY;
	   i < T_SIZE; ++i)
	{
	  if (weight_mat[i][j] > -INFINITY)
	    {
	      if (weight_mat[i][j] < min_score[j])
		min_score[j] = weight_mat[i][j];
	      if (weight_mat[i][j] > max_score[j])
		max_score[j] = weight_mat[i][j];
	    }
	  else Minus_inf = YES;
	}
      if (min_score[j] < INFINITY)
	{
	  Max_exact_score += max_score[j];
	  Min_exact_score += min_score[j];
	}
    }

#ifdef PRINTRANGE
  /* Print the maximum score, the minimum score, and the range of scores. */
  printf("                                    maximum score: %8.3f\n",
	 Max_exact_score);
  if (Minus_inf == NO)
    printf("                                    minimum score: %8.3f\n",
	   Min_exact_score);
  else
    {
      printf("                                    minimum score: -INFINITY\n");
      printf("                        minimum score > -INFINITY: %8.3f\n",
	     Min_exact_score);
    }
  printf("            range of scores: %8.3f - %8.3f = %8.3f\n",
	 Max_exact_score, Min_exact_score, Max_exact_score - Min_exact_score);
  printf("\n");
#endif


  /* Determine the minimum score for calculating p-values. */
  if (Min_score >= Min_exact_score) min_prob_score = Min_score;
  else min_prob_score = Min_exact_score;

  /* Make sure the minimum score for calculating p-values is
   * not greater than the maximum possible score. */
  if (min_prob_score > Max_exact_score)
    {
      fprintf(stderr, "The minimum score for calculating p-values is greater");
      fprintf(stderr, "\nthan the maximum possible score.\n");
      exit(1);
    }
  /* Adjustment for the special case when the minimum score for
   * calculating p-values is equal to the maximum possible score. */
  else if (min_prob_score == Max_exact_score) max_min_equality = YES;
  else max_min_equality = NO;

  /* Determine "Multiple". */
  if (max_min_equality == YES) Multiple = Max_int_range / Max_exact_score;
  else Multiple = Max_int_range / (Max_exact_score - min_prob_score);

  /* Determine "STmax_score[]", "STmin_score[]", and Max_int_score. */
  for (j = 0, min_int_score = Max_int_score = 0; j < width; ++j)
    {
      if (min_score[j] < INFINITY)
	{
	  STmax_score[j] = ROUND_TO_INT(Multiple * max_score[j]);
	  STmin_score[j] = ROUND_TO_INT(Multiple * min_score[j]);

	  Max_int_score += STmax_score[j];
	  min_int_score += STmin_score[j];
	}
      else STmax_score[j] = STmin_score[j] = 0;
    }

  /* Determine "Min_int_score". */
  Min_int_score = ROUND_TO_INT(Multiple * min_prob_score);
  if ((Min_exact_score >= Min_score) && (min_int_score < Min_int_score))
    Min_int_score = min_int_score;
  --Min_int_score;
  if ((max_min_equality == YES) || (Min_int_score >= Max_int_score))
    Min_int_score = Max_int_score - 1;

  /* Determine "STmax_int_range". */
  STmax_int_range = Max_int_score - Min_int_score;
  if (STmax_int_range < 1) bug_report("determine_variables()");

#ifdef TEST
  /* Assign values to the -INFINITY elements of the weight matrix.
   * Returns Min_exact_score if there are no -INFINITY scores; otherwise,
   * returns a lower score determined by the values assigned
   * to the -INFINITY elements of the weight matrix. */
  Min_pseudo_score = adj_minus_infinity(weight_mat, max_score,
					min_score, width);
#endif


  free_error(max_score, "max_score", "determine_variables()");
  free_error(min_score, "min_score", "determine_variables()");
}

/* Assign values to the -INFINITY elements of the weight matrix.
 * Returns Min_exact_score if there are no -INFINITY scores; otherwise,
 * returns a lower score determined by the values assigned
 * to the -INFINITY elements of the weight matrix. */
double adj_minus_infinity(weight_mat, max_score, min_score, width)
     double **weight_mat;    /* The weight matrix. */
     double *max_score;      /* Value of column element with largest score. */
     double *min_score;      /* Value of column element with lowest score. */
     int width;              /* Width of the matrix. */
{
  int i, j;
  double min_pseudo_score;
  double inf_offset;         /* Number subtracted from maximum scores to
			      * establish the effective -INFINITY score. */


  if (Minus_inf == NO) min_pseudo_score = Min_exact_score;

  else
    {
      inf_offset = Max_exact_score - Min_exact_score + 1.1 / Multiple;

      for (j = 0, min_pseudo_score = 0.0; j < width; ++j)
	{
	  for (i = 0; i < T_SIZE; ++i)
	    {
	      if (weight_mat[i][j] == -INFINITY)
		{
		  weight_mat[i][j] = max_score[j] - inf_offset;
		  min_score[j] = weight_mat[i][j];
		}
	    }
	  min_pseudo_score += min_score[j];
	}
    }

  return(min_pseudo_score);
}

/* Allocate space for "STarrayA[]", "STarrayB[]", "STarrayC[]",
 * and "STarray_index[]". */
void appr_mat_space()
{
  int i;
  int max_score_1 = STmax_int_range + 1;  /* Maximum size of arrays. */


  STarrayA = (double *)recalloc_error(STarrayA, max_score_1, sizeof(double),
				    "STarrayA", "appr_mat_space()");

  STarrayB = (double *)calloc_error(max_score_1, sizeof(double),
				    "STarrayB", "appr_mat_space()");

  STarrayC = (double *)calloc_error(max_score_1, sizeof(double),
				    "STarrayC", "appr_mat_space()");

  STarray_index = (int *)calloc_error(T_SIZE, sizeof(int),
				      "STarray_index", "appr_mat_space()");

  for (i = 0; i < max_score_1; ++i)
    {
      STarrayA[i] = -INFINITY;
      STarrayC[i] = -INFINITY;
      STarrayB[i] = 0.0;
    }
}


/* Determine the logarithm of the distribution of the integral scores of
 * the "weight_mat[][]" and place in ArrayA[]. */
void determine_distribution(weight_mat, width)
     double **weight_mat;     /* The weight matrix. */
{
  int i, j;
  double ln_prob_AB;       /* ln(column_probability in STarrayB[] multiplied
			    *    by previous_probability in STarrayA[]) */
  double ln_prob_C;        /* The ln(probability currently in STarrayC[]). */
  double *temp_array;      /* Temporary pointer when exchanging pointers to
			    * "STarrayA" and "STarrayC". */
  int max_score;           /* Maximum score prior to current multiplication. */
  int new_max_score;       /* Maximum score after current multiplication. */
  int min_score;           /* Minimum score prior to current multiplication. */
  int new_min_score;       /* Minimum score after current multiplication. */
  int diff_sub;      /* Number of different substitution scores (<= A_size). */

  /* Pointers to multiplication arrays so maximum index = maximum score. */
  double *arrayA;       /* Pointer to STarrayA: product of previous columns. */
  double *arrayB;       /* Pointer to STarrayB: current column. */
  double *arrayC;       /* Pointer to STarrayC: product of current column
			 *                      with previous columns. */
  /* Minimum index in each of the multiplication arrays. */
  int min_idx_A;        /* arrayA */
  int min_idx_B;        /* arrayB */
  int min_idx_C;        /* arrayC */

  /* The score index for each multiplication array. */
  int idx_A;            /* arrayA */
  int idx_B;            /* arrayB */
  int idx_C;            /* arrayC */


  /* Initialize "arrayA", "max_score", and "min_score" for
   * the first matrix column.*/
  max_score = STmax_score[0];
  min_score = STmin_score[0];
  min_idx_A = -(STmax_int_range - max_score);
  if (min_score > min_idx_A) min_idx_A = min_score - 1;
  arrayA = -min_idx_A + STarrayA;
  for (i = 0; i < T_SIZE; ++i)
    {
      if (weight_mat[i][0] <= -INFINITY) idx_A = min_idx_A;
      else
	{
	  idx_A = ROUND_TO_INT(Multiple * weight_mat[i][0]);
	  if (idx_A < min_idx_A) idx_A = min_idx_A;
	}

      if (arrayA[idx_A] > -INFINITY)
	arrayA[idx_A] = SUM_LN(Ln_P[i], arrayA[idx_A]);
      else arrayA[idx_A] = Ln_P[i];
    }

  for (j = 1; j < width; ++j)
    {
      /* Insert the probabilities for the current column into "arrayB[]".
       * Set the "STarray_index[]" to indicate once which scores occur
       * in the current column. */
      min_idx_B = STmax_score[j] - STmax_int_range;
      arrayB = -min_idx_B + STarrayB;
      for (i = 0, diff_sub = 0; i < T_SIZE; ++i)
	{
	  if ((weight_mat[i][j] <= -INFINITY) ||
	     ((idx_B = ROUND_TO_INT(Multiple * weight_mat[i][j])) < min_idx_B))
	    idx_B = min_idx_B;

	  if (arrayB[idx_B] == 0.0)
	    {
	      arrayB[idx_B] = P[i];
	      STarray_index[diff_sub++] = idx_B;
	    }
	  else arrayB[idx_B] += P[i];
	}
      /* Convert the probabilities in "arrayB[]" to logarithms. */
      for (i = 0; i < diff_sub; ++i)
	arrayB[STarray_index[i]] = log(arrayB[STarray_index[i]]);


      /* Multiply "arrayA[]" by "arrayB[]". */
      new_max_score = max_score + STmax_score[j];
      new_min_score = min_score + STmin_score[j];
      min_idx_C = -(STmax_int_range - new_max_score);
      if (new_min_score > min_idx_C) min_idx_C = new_min_score - 1;
      arrayC = -min_idx_C + STarrayC;
      if (arrayA[min_idx_A] > -INFINITY)
	{
	  arrayC[min_idx_C] = arrayA[min_idx_A] + LN_PROB_SUM;
	  arrayA[min_idx_A] = -INFINITY;
	}
      for (idx_A = max_score; idx_A > min_idx_A; --idx_A)
	{
	  if (arrayA[idx_A] > -INFINITY)
	    {
	      for (i = 0; i < diff_sub; ++i)
		{
		  idx_B = STarray_index[i];

		  if ((idx_C = idx_A + idx_B) < min_idx_C) idx_C = min_idx_C;

		  ln_prob_AB = arrayA[idx_A] + arrayB[idx_B];
		  ln_prob_C = arrayC[idx_C];

		  if (ln_prob_C <= -INFINITY) arrayC[idx_C] = ln_prob_AB;
		  else arrayC[idx_C] = SUM_LN(ln_prob_AB, ln_prob_C);
		}

	      /* Zero the current element of "arrayA[]". */
	      arrayA[idx_A] = -INFINITY;
	    }
	}


      /* Zero "arrayB[]". */
      for (i = 0; i < diff_sub; ++i) arrayB[STarray_index[i]] = 0.0;

      /* Exchange pointers to "STarrayA[]" and "STarrayC[]". */
      temp_array = STarrayA;
      STarrayA = STarrayC;
      STarrayC = temp_array;
      arrayA = arrayC;
      min_idx_A = min_idx_C;

      /* Update "max_score". */
      max_score = new_max_score;
      min_score = new_min_score;
    }

  /* Transfer information on arrayA into external variables.
   * Max_int_score: the maximum index determined in determine_variables().
   * Min_int_score: the minimum index might change from value set
   *                in determine_variables() due to rounding error.
   *      ArrayA[]: holds final ln(probabilities) of integer scores. */
  Min_int_score = min_idx_A;    
  ArrayA = arrayA;
}


/* Determine the logarithm of the p-values and place in Marginal_prob[]. */
void determine_marginal_prob()
{
  int i, j, k;
  int max_score;    /* The greater of the maximum integral score determined
		     * from the maximum exact score and the integral matrix. */
  double delta_ln_p;/* The change in the ln(p-value) between
		     * two observed integral scores. */


  /* Allocate space for "Marginal_prob[]". */
  max_score = ROUND_TO_INT(Multiple * Max_exact_score);
  if (Max_int_score > max_score) max_score = Max_int_score;
  Marginal_prob = (double *)recalloc_error(Marginal_prob,
				 max_score - Min_int_score + 1, sizeof(double),
				 "Marginal_prob", "determine_marginal_prob()");
  Marginal_prob -= Min_int_score;


  /* Determine the ln(p-value) of the Max_int_score. */
  Marginal_prob[Max_int_score] = ArrayA[Max_int_score];

  /* Copy ln(p-values) for scores greater than Max_int_score. */
  for (i = Max_int_score + 1; i <= max_score; ++i)
    Marginal_prob[i] = Marginal_prob[Max_int_score];


  /* Determine the main part of Marginal_prob[]. */
  for (i = Max_int_score - 1, j = Max_int_score;
       i > Min_int_score; --i)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  Marginal_prob[i] = SUM_LN(ArrayA[i], Marginal_prob[j]);

	  /* Interpolate ln(p-values) for unobserved scores. */
	  delta_ln_p = (Marginal_prob[i] - Marginal_prob[j]) / (double)(j - i);
	  for (k = j--; j > i; k = j--)
	    Marginal_prob[j] = Marginal_prob[k] + delta_ln_p;
	}
    }

  /* Copy ln(p-values) for scores less than the lowest observed score j. */
  for (k = j - 1; k > Min_int_score; --k) Marginal_prob[k] = Marginal_prob[j];

  /* Determine the ln(p-value) for the Min_int_score. */
  Marginal_prob[i] = SUM_LN(ArrayA[i], Marginal_prob[j]);


#ifdef TEST
  /* Print the ln(sum of the distribution)---better be close to ln(1) = 0. */
  printf("    ln(the sum of the distribution): %g\n",
	 Marginal_prob[Min_int_score]);

  /* Print the maximum and minimum scores. */
  printf("  Maximum score (true/approximated): %8.3f/%-8.3f\n",
	 Max_exact_score, (double)Max_int_score / Multiple);
  printf("                      Minimum score: %8.3f\n", Min_exact_score);
  printf("\n");
#endif

#ifdef PRINTRANGE
  /* Print minimum score for numerically estimating p-values. */
  printf("           minimum score for calculating p-values: %8.3f\n",
	 (double)(Min_int_score + 1) / Multiple);

  /* Print the maximum and minimum numerically calculated p-values. */
  printf("       maximum ln(numerically calculated p-value): %8.3f\n",
	 Marginal_prob[Min_int_score + 1]);
  printf("       minimum ln(numerically calculated p-value): %8.3f\n",
	 Marginal_prob[Max_int_score]);

  printf("\n");
#endif


#ifdef PRINTDISTRIBUTION
  printf("a priori letter probabilities:\n");
  for (i = 0; i < T_SIZE; ++i) printf("%5.3f\n", P[i]);
  printf("\n");

  printf("Multiple: %.15g\n", Multiple);
  printf("\n");

  printf("integral score | probability | p-value | ln(p-value)\n");
  for (i = max_score; i > Max_int_score; i--)
    printf("%11d:          0 %10g %11f\n",
	   i, exp(Marginal_prob[i]), Marginal_prob[i]);
  for (i = Max_int_score; i >= Min_int_score; i--)
    printf("%11d: %10g %10g %11f\n",
	   i, exp(ArrayA[i]), exp(Marginal_prob[i]), Marginal_prob[i]);
#endif
}

/* Determine and print the cutoff score given the logarithm
 * of the target p-value.  Returns the cutoff score. */
double cutoff_approx(ln_prob, note)
     double ln_prob;    /* The target p-value. */
     char *note;        /* Note describing the source of the cutoff p-value. */
{
  int i, j, k;
  double cutoff_score;     /* Numerically calculated cutoff score. */
  double ln_prob_score;    /* ln(probability * score) */
  double avg;              /* Average score above the cutoff
			    * [initially holds ln(average)]. */


  /* Determine the cutoff score and the average score. */
  for (i = Max_int_score, j = Max_int_score + 1,
       k = Max_int_score - Min_int_score, avg = -INFINITY;
       (k > 0) && (Marginal_prob[i] <= ln_prob); --i, --k)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  avg = SUM_LN(avg, ln_prob_score);
	}
    }
#ifndef HIGH_CUTOFF
  if (Marginal_prob[j] < ln_prob)
    {
      for ( ; (k > 0) && (ArrayA[i] <= -INFINITY); --i, --k);

      if (k > 0)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  avg = SUM_LN(avg, ln_prob_score);
	}
    }
#endif
  if (j <= Max_int_score)
    avg = (exp(avg - Marginal_prob[j]) + (double)Min_int_score) / Multiple;
  else avg = 0.0;


  /* Print the target cutoff p-value. */
  printf("ln(cutoff p-value)%s: %8.3f\n", note, ln_prob);

  /* Print the cutoff score and the achieved p-value. */
  if (j > Max_int_score)
    {
      cutoff_score = (double)Max_int_score / Multiple;
      printf("      numerically calculated cutoff score is the ");
      printf("maximum score: %8.3f\n", cutoff_score);

      printf("       minimum ln(numerically calculated p-value): %8.3f\n",
	     Marginal_prob[Max_int_score]);
    }
  else if ((i == Min_int_score) &&
	   (ArrayA[Min_int_score] > -INFINITY) &&
	   (Marginal_prob[j] < ln_prob))
    {
      /* Information for the minimum score for numerically
       * estimating p-values is now printed automatically above. */
      cutoff_score = -INFINITY;
      printf("The maximum calculated p-value is less than the cutoff\n");
      printf("p-value%s.\n", note);
      printf("Try decreasing the minimum score for calculating p-values\n");
      printf("using the \"-M\" command line option.\n");
    }
  else
    {
      cutoff_score = (double)j / Multiple;
      printf("              numerically calculated cutoff score: %8.3f\n",
	     cutoff_score);

      printf("        ln(numerically calculated cutoff p-value): %8.3f\n",
	     Marginal_prob[j]);
    }

  /* Print the average score above the cutoff. */
  printf("average score above numerically calculated cutoff: %8.3f\n", avg);
  printf("\n");

  return(cutoff_score);
}


/* Determine and print the cutoff score given the average score of
 * the L-mers above the cutoff. */
double cutoff_approx_avg(avg_score)
     double avg_score;       /* The target average score. */
{
  int i, j, k;
  double ln_int_avg_score;/* ln(target average score converted to "integer") */
  double ln_prob_score;       /* ln(probability * score) */
  double ln_prob_score_sum;   /* ln(SUM ln_prob_score). */
  double ln_int_avg;          /* ln(average score above the current cutoff). */
  double avg;                 /* Average score above the cutoff. */
  double cutoff_score;        /* The cutoff score. */


  /* Convert the targeted average score to its "integral" representation
   * and take its logarithm. */
  if (avg_score * Multiple < (double)(Min_int_score + 1))
    {
      printf("minimum numerically calculated cutoff score: ");
      printf("%8.3f\n", (double)(Min_int_score + 1) / Multiple);
      printf("    targeted average score above the cutoff: %8.3f\n",
	     avg_score);
      printf("\n");

      return(-INFINITY);
    }
  ln_int_avg_score = log(avg_score * Multiple - (double)Min_int_score);

  /* Determine the integral cutoff score and the average score. */
  for (i = Max_int_score, j = Max_int_score + 1,
       k = Max_int_score - Min_int_score,
       ln_int_avg = INFINITY, ln_prob_score_sum = -INFINITY;
       (k > 0) && (ln_int_avg > ln_int_avg_score); --i, --k)
    {
      if (ArrayA[i] > -INFINITY)
	{
	  j = i;
	  ln_prob_score = ArrayA[i] + log((double)k);
	  ln_prob_score_sum = SUM_LN(ln_prob_score_sum, ln_prob_score);

	  ln_int_avg = ln_prob_score_sum - Marginal_prob[j];
	}
    }
  avg = (exp(ln_int_avg) + (double)Min_int_score) / Multiple;

  if ((k == 0) && (ArrayA[Min_int_score] > -INFINITY) &&
      (ln_int_avg > ln_int_avg_score))
    {
      printf("THE MINIMUM SCORE FOR CALCULATING PROBABILITIES\n");
      printf("IS TOO HIGH TO FIND THE TARGETED AVERAGE SCORE.\n");
    }

  /* Determine the cutoff score. */
  cutoff_score = (double)j / Multiple;

  /* Print the target average score, the cutoff score, the p-value,
   * and the calculated average score. */
  printf("  targeted average score above the cutoff: %8.3f\n", avg_score);
  printf("      numerically calculated cutoff score: %8.3f\n", cutoff_score);
  printf("                   ln(calculated p-value): %8.3f\n",
	 Marginal_prob[j]);
  printf("calculated average score above the cutoff: %8.3f\n", avg);
  printf("\n");

  return(cutoff_score);
}

/* Extract the p-value from Marginal_prob[] for the given score.
 * Print the p-value. */
void extract_marginal_prob(score)
     double score;
{
  int int_score = ROUND_TO_INT(Multiple * score);  /* Integral score. */
  double ln_prob;                                  /* ln(p-value) */


  if (Max_int_range > 0.0)
    {
      if ((score >= -INFINITY / 2.0) &&
	  (score >= Min_score - Equal_err) &&
#ifdef TEST
	  (score >= Min_exact_score - 1.0 / Multiple) &&
#endif
	  (score >= Min_exact_score - Equal_err))
	{
	  if (int_score > Min_int_score) ln_prob = Marginal_prob[int_score];
	  else ln_prob = Marginal_prob[Min_int_score + 1];

	  printf("  ln(p-value)= %7.2f", ln_prob);
	}

      else if (Min_score <= Min_exact_score)
	{
	  ln_prob = 0.0;
	  printf(" ln(p-value)= %5.2f", ln_prob);
	}

#if 0
      else printf("  ln(p-value)> %7.2f", Marginal_prob[Min_int_score + 1]);
#endif
    }
}
