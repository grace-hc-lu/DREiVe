/* Copyright 1990, 1994, 1996, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


#include "options.h"
#include "heap.h"


/* Determine the score of the L-mers of the current sequence.
 * The following external variables are used:
 *
 *    Comp_status           A_comp                File
 *    Thresh     	    S_fp  		  Fp
 *    L_thresh   	    Matrix		  Circle
 *    U_thresh   	    L     		  Sequence
 *    Top        	    L_1   		  Seq_length
 *    A_size                Equal_err
 */


/* Function for determing whether a score should be printed.  Equals either
 * "print_threshold()" or "save_top()".
 * RETURN VALUES
 *     YES: continue processing the current sequence.
 *      NO: stop the processing of the current sequence. */
static char (*print_option)();

static int STmax_idx; /* Maximum sequence index at which an L-mer can begin. */

/* Variables for printing a list of top scores. */
static int STarray_size;      /* Size of the array holding the information
			       * about top scores when "Top" > 0. */
static int STnum_top;         /* The number of top scoring L-mers in
			       * the current sequence. */
static int STnum_max_top;     /* The number of top scoring L-mers,
			       * excluding the minimum top score. */
static PAIR *STtop_array;     /* Array holding the information about
			       * top scores when "Top" > 0. */


/* Initialize information for "process_sequence()". */
void init_variables()
{
  char print_threshold();  /* Print score if within threshold limits. */
  char save_top();         /* Save score if tops for current sequence. */


  /* Initialize "(*print_option)()". */
  if (Top == 0) print_option = print_threshold;
  else print_option = save_top;

  if (Top > 0)
    {
      STarray_size = Top + 1;
      STtop_array = (PAIR *)calloc_error(STarray_size, sizeof(PAIR),
					 "STtop_array", "init_variables()");
      STtop_array[0].score = -INFINITY;
    }
}



/* Determine the score of the L-mers of the current sequence. */
void process_sequence()
{
  int i, j;
  int circle_size;    /* Number of letters within circular wrap around. */

  int score_l_mer();       /* Score the current L-mer. */
  void print_top_scores(); /* Print top scores as stored in "STtop_array[]". */


  /* Initialize the number of top scores. */
  STnum_top = 0;

  /* Determine the size of the circular wraparound. */
  if (Circle == YES)
    {
      if (Seq_length >= L) circle_size = L - 1;
      else circle_size = Seq_length;
    }
  else circle_size = 0;

  /* Add the circular wraparound to the Sequence[]. */
  Sequence = (LETTER *)recalloc_error((char *)Sequence, Seq_length+circle_size,
			     sizeof(LETTER), "Sequence", "process_sequence()");
  for (i = 0, j = Seq_length; i < circle_size; ++i, ++j)
    Sequence[j] = Sequence[i];


  /* Process each L-mer of the Sequence[]. */
  for (i = 0, STmax_idx = Seq_length + circle_size - L; i <= STmax_idx; ++i)
    i = score_l_mer(i);                 /* Score the current L-mer. */

  /* Print results if only top scores for each sequence are being printed. */
  if (Top > 0) print_top_scores();

  if (Fp != S_fp) fclose(Fp);
  return;
}


/* Score the current L-mer and its complement, if desired.
 * Print the score, depending on the threshold conditions
 * and whether only the top scores are being printed.
 * Returns the index preceding the next L-mer. */
int score_l_mer(idx)
     int idx;              /* The starting index for the current L-mer. */
{
  int i, j;
  LETTER *l_mer;           /* The L-mer currently being scored. */
  int letter;              /* The current letter of the current L-mer. */
  double score;            /* Score of current L-mer or its complement. */


  /* Determine the score of the current L-mer. */
  for (i = 0, l_mer = Sequence + idx, score = 0; i < L; ++i)
    {
      letter = (int)(l_mer[i]);
      if (letter >= A_size) return(idx + i);

      score += Matrix[letter][i];
    }

  /* Print score depending on threshold conditions. */
  if ((*print_option)(1, idx, score) == NO) return(STmax_idx);


  if (Comp_status == YES)
    {
      /* Determine the score of the current L-mer's complement. */
      for (i = 0, j = L_1, score = 0; i < L; ++i, --j)
	{
	  letter = A_comp[(int)(l_mer[j])];
	  score += Matrix[letter][i];
	}

      /* Print score depending on threshold conditions. */
      if ((*print_option)(-1, idx, score) == NO) return(STmax_idx);
    }

  return(idx);
}


/* "(*print_option)()" is the function for determing whether a score
 * should be printed.  Equals either
 *     print_threshold(): Top == 0
 *     save_top(): Top > 0.
 * RETURN VALUES
 *     YES: continue processing the current sequence.
 *      NO: stop the processing of the current sequence.
 */

/* Determine whether the score of the current L-mer should be printed.
 * YES: continue processing the current sequence (always returns YES). */
char print_threshold(comp_strand, idx, score)
     int comp_strand;   /* Flag for whether L-mer is from complementary strand.
			 * 1: not from complement; -1: from complement. */
     int idx;           /* The starting index for the current L-mer. */
     double score;      /* Score of current L-mer or its complement. */
{
  int status = NO;    /* Status of whether score is within threshold limits. */
  int position;         /* Position of the score within the sequence. */
  void print_score();   /* Print a score. */


  /* Determine whether the score is within the threshold limits. */
  switch (Thresh)
    {
    case '\0':
      status = YES;
      break;
    case 'l':
      if (score >= L_thresh)
	status = YES;
      break;
    case 'u':
      if (score < U_thresh)
	status = YES;
      break;
    case 'b':
      if ((score >= L_thresh) && (score < U_thresh))
	status = YES;
      break;
    default:
      fprintf(stderr, "PROGRAM BUG: print_option\n");
    }

  /* Print the results if "status" equals YES. */
  if (status == YES)
    {
      position = idx + 1;
      if (comp_strand < 0) position = -position;

      print_score(File, position, score);
    }

  return(YES);
}

/* Determine whether the score of the current L-mer is a top score
 * and should be printed.  RETURN VALUES
 *     YES: continue processing the current sequence.
 *      NO: stop the processing of the current sequence
 *          if upper threshold has been exceeded. */
char save_top(comp_strand, idx, score)
     int comp_strand;   /* Flag for whether L-mer is from complementary strand.
			 * 1: not from complement; -1: from complement. */
     int idx;           /* The starting index for the current L-mer. */
     double score;      /* Score of current L-mer or its complement. */
{
  PAIR pair;  /* New position and score information for inserting into heap. */

  /* Compare 2 PAIR structures bases on the score member (lower < higher).
     Return values causes the heap to store the highest scores.
     Returns -1 if pair2 < pair1;
            0 if pair2 = pair1;
            1 if pair2 > pair1. */
  int compar_score(const void *pair1, const void *pair2);



  /* Return to calling function if score is less than the lower threshold. */
  if (((Thresh == 'l') || (Thresh == 'b')) && (score < L_thresh))
    return(YES);

  /* Quit processing the sequence if greater than the upper threshold. */
  if (((Thresh == 'u') || (Thresh == 'b')) && (score >= U_thresh))
    {
      STnum_top = 0;
      return(NO);
    }


  /* The first "Top" L-mers are automatically the top scoring L-mers so far. */
  if (STnum_top < Top)
    {
      ++STnum_top;
      STtop_array[STnum_top].position = comp_strand * (idx + 1);
      STtop_array[STnum_top].score = score;

      if (STnum_top == Top)
	tieheap_initialize(STtop_array, STnum_top, &STnum_max_top,
			   sizeof(PAIR), compar_score);
    }

  /* Add a new L-mer to the tieheap if its score is >= the current
     lowest score.  Delete elements and readjust the tieheap as necessary. */
  else if (score >= STtop_array[STnum_top].score - Equal_err)
    {
      pair.position = comp_strand * (idx + 1);
      pair.score = score;

      STtop_array = tieheap_add(&pair, STtop_array, &STnum_top, &STnum_max_top,
			      &STarray_size, Top, sizeof(PAIR), compar_score);
    }


  return(YES);
}

/* Compare 2 PAIR structures bases on the score member (lower < higher).
   Return values causes the heap to store the highest scores.
   Returns -1 if pair2 < pair1;
            0 if pair2 = pair1;
            1 if pair2 > pair1. */
int compar_score(const void *pair1, const void *pair2)
{
  double score1 = ((PAIR *)pair1)->score;
  double score2 = ((PAIR *)pair2)->score;


  if (fabs(score2 - score1) <= Equal_err) return(0);
  else if (score2 < score1) return(-1);
  else if (score2 > score1) return(1);

  else bug_report("compar_score()");
}


/* Compare 2 PAIR structures based on the score member.
   The primary comparison is based on the score (lower < higher). The secondary
   comparison is based on the absolute value of the position (higher < lower).
   The tertiary comparison is based on the sign (negative < positive).
   Return values causes the heap to store the highest scores.
   Returns -1 if pair2 < pair1;
            0 if pair2 = pair1;
            1 if pair2 > pair1. */
int compar_score_2(const void *pair1, const void *pair2)
{
  double score1 = ((PAIR *)pair1)->score;
  double score2 = ((PAIR *)pair2)->score;
  int position1;
  int position2;


  if (fabs(score2 - score1) <= Equal_err)
    {
      position1 = ((PAIR *)pair1)->position;
      position2 = ((PAIR *)pair2)->position;

      /* If absolute value of the position members differ. */
      if (abs(position1) < abs(position2)) return(-1);
      else if (abs(position1) > abs(position2)) return(1);

      /* If absolute value of the position members are the same. */
      else if (position1 == position2) return(0);
      else if (position2 < position1) return(-1);
      else if (position2 > position1) return(1);
    }

  else if (score2 < score1) return(-1);
  else if (score2 > score1) return(1);

  bug_report("compar_score_2()");
}

/* Print the top scores as stored in "STtop_array[]". */
void print_top_scores()
{
  int i;
  int position;         /* The position of the current L-mer. */
  double score;         /* Score of current L-mer or its complement. */
  void print_score();   /* Print a score. */

  /* Compare 2 PAIR structures bases on the position member.  The
     primary comparison is bases on the absolute value of the position
     (lower < higher), the secondary comparison is based on the sign
     (positive < negative).  Return values causes hsort() to sort in
     ascending order.
     Returns -1 if pair1 < pair2;
              0 if pair1 = pair2;
              1 if pair1 > pair2. */
  int compar_position(const void *pair1, const void *pair2);

/* Compare 2 PAIR structures based on the score member.
   The primary comparison is based on the score (lower < higher). The secondary
   comparison is based on the absolute value of the position (higher < lower).
   The tertiary comparison is based on the sign (negative < positive).
   Return values causes the heap to store the highest scores.
   Returns -1 if pair2 < pair1;
            0 if pair2 = pair1;
            1 if pair2 > pair1. */
  int compar_score_2(const void *pair1, const void *pair2);


  /* Sort the STtop_array according to increasing position number. */
  if (Print_order == 0)
    hsort(STtop_array + 1, STnum_top, sizeof(PAIR), compar_position);

  /* Sort the STtop_array according to decreasing score. */
  else if (Print_order == 1)
    hsort(STtop_array + 1, STnum_top, sizeof(PAIR), compar_score_2);

  else bug_report("print_top_scores()");
  

  for (i = 1; i <= STnum_top; ++i)
    {
      position = STtop_array[i].position;
      score = STtop_array[i].score;

      print_score(File, position, score);
    }
}

/* Compare 2 PAIR structures bases on the position member.  The
   primary comparison is bases on the absolute value of the position
   (lower < higher), the secondary comparison is based on the sign
   (positive < negative).  Return values causes hsort() to sort in
   ascending order.
   Returns -1 if pair1 < pair2;
            0 if pair1 = pair2;
            1 if pair1 > pair2. */
int compar_position(const void *pair1, const void *pair2)
{
  int position1 = ((PAIR *)pair1)->position;
  int position2 = ((PAIR *)pair2)->position;


  /* If absolute value of the position members differ. */
  if (abs(position1) < abs(position2)) return(-1);
  else if (abs(position1) > abs(position2)) return(1);

  /* If absolute value of the position members are the same. */
  else if (position1 == position2) return(0);
  else if (position1 > position2) return(-1);
  else if (position1 < position2) return(1);

  else bug_report("compar_position()");
}


/* Print a score. */
void print_score(seq_name, position, score)
     char *seq_name;           /* Name of the sequence being scored. */
     int position;             /* Position being scored.  A negative position
				* number means the score is on the
				* complementary strand. */
     double score;             /* The score. */
{
  void extract_marginal_prob(); /* Extract and print marginal
				 * probability from Marginal_prob[]. */


  /* L-mer is not on the complementary strand. */
  if (position > 0)
    printf("%*s  position= %6d  score= ",
	   NAME_LENGTH, seq_name, position);

  /* L-mer is on the complementary strand. */
  else
    printf("%*s  position= %6dC score= ",
	   NAME_LENGTH, seq_name, abs(position));

  if ((score < -INFINITY / 2.0) ||
#ifdef TEST
      (score < Min_exact_score - 1.0 / Multiple) ||
#endif
      (score < Min_exact_score - Equal_err))
    printf("-INFINITY");

  else printf("%6.2f", score);

  extract_marginal_prob(score);
  printf("\n");
}
