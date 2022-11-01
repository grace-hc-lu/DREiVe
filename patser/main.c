/* Copyright 1990, 1996, 2000 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


/* Gerald Z. Hertz Fri Jul  6 18:18 1990
 * This program scores the L-mers of the indicated sequences
 * against the indicated summary matrix. */

#include "options.h"

main(argc, argv)
     int argc;
     char *argv[];
{
  int i;
  double cutoff_score;    /* Cutoff score based on a ln(p-value). */

  void command_line();    /* Read the command line arguments. */
  void process_matrix();  /* Convert summary matrix to a specificity matrix. */
  void read_matrix();     /* Read the specificity matrix. */
  void init_read_seq();   /* Initialize information for "read_sequence()". */
  void init_variables();  /* Initialize information for "process_sequence()".*/
  int read_sequence();    /* Read a sequence; determine whether circular. */
  void process_sequence();/* Determine score of L-mers of the sequence. */

  void det_marginal_prob(); /* Approximate marginal probabilities of scores. */
  double cutoff_approx();   /* Determine cutoff score from Marginal_prob[]. */
  void cutoff_exact();      /* Determine the exact cutoff score. */
  void det_avg_var();       /* Determine avg, var, max, and min scores. */
  void cutoff_analytical(); /* Determine the cutoff score analytically. */

  /* Determine the cutoff scores by trying to get the average score
   * above the cutoff to equal the information content. */
  double cutoff_approx_avg(); /* Determine cutoff score from Marginal_prob[].*/
  void heap_sort();              /* Sort the heap in descending order. */
  double cutoff_exact_avg();     /* Determine the exact cutoff score. */
  double cutoff_analytical_avg();/* Determine the cutoff score analytically. */

  /* Analytical calculations at high standard deviation. */
  void func_avg_score_1();       /* avg score as a function of cutoff. */
  double func_ln_prob_1();       /* Marginal probability of a score. */
  /* Analytical calculations at low standard deviation: uses gamma functions */
  void func_avg_score_2();       /* avg score as a function of cutoff. */
  double func_ln_prob_2();       /* Marginal probability of a score. */
  /* Analytical calculations at low standard deviation: uses gamma functions
   * and an explicit maximum score. */
  void func_avg_score_3();       /* avg score as a function of cutoff. */
  double func_ln_prob_3();       /* Marginal probability of a score. */
  void det_marginal_prob_analytical();  /* Determine marginal probabilities
					 * analytically. */


  /* Read the command line arguments.  The default settings are contained in
   * the "options.c" file.  See "directions" for more details. */
  command_line(argc, argv);

  /* Determine Ln_P[] = ln(P[]). */
  Ln_P = (double *)calloc_error(A_size, sizeof(double), "Ln_P", "main()");
  for (i = 0; i < A_size; ++i) Ln_P[i] = log(P[i]);

  if (Sum_mat == YES)
    /* Convert the summary matrix to a specificity matrix. */
    process_matrix();
  else
    /* Read the weight matrix. */
    read_matrix();

  /* Create "Marginal_prob[]", array containing marginal probabilities
   * for scores after converting the weight matrix to integers. */
  if (Max_int_range > 0.0)
    {
      det_marginal_prob(Matrix, L);

      if (Sum_mat == YES)
	{
	  /* Print cutoff score based on sample size
	   * adjusted information content. */
	  cutoff_score = cutoff_approx(-Inf_adj,
			" based on sample size adjusted information content");
	  if (Auto_cutoff == YES)
	    {
	      L_thresh = cutoff_score;
	      if (Thresh == '\0') Thresh = 'l';
	      else if (Thresh == 'u') Thresh = 'b';
	      else bug_report("main() 1");
	    }

#ifdef TEST
	  /* Determine cutoff scores based on a target marginal probability of
	   * negative the information content of the alignment matrix after
	   * adding fudge factors. */
	  if (Test >= 1) cutoff_approx(-Inf_fudge,
			" based on pseudo-count adjusted information content");
	  if (Test == 2) cutoff_exact(Matrix, -Inf_fudge);
	  if (Test >= 1)
	    {
	      det_avg_var(Align_mat, Matrix);
	      cutoff_analytical(-Inf_fudge);
	    }

	  /* Determine cutoff scores based on a target average score equal to
	   * the information content of the alignment matrix after
	   * adding fudge factors. */
	  if (Test >= 1) cutoff_approx_avg(Inf_fudge);
	  if (Test == 2)
	    {
	      heap_sort();
	      cutoff_exact_avg(Inf_fudge);
	    }
	  if (Test >= 1)
	    {
	      printf("Analytical formula for high standard deviation.\n");
	      cutoff_analytical_avg(Inf_fudge,func_avg_score_1,func_ln_prob_1);

	      printf("Gamma analytical formula for low standard deviation.\n");
	      cutoff_analytical_avg(Inf_fudge,func_avg_score_2,func_ln_prob_2);

	      printf("Gamma analytical formula for low standard deviation ");
	      printf("and explicit maximum score.\n");
	      cutoff_analytical_avg(Inf_fudge,func_avg_score_3,func_ln_prob_3);

	      /* Determine marginal probabilities analytically. */
	      det_marginal_prob_analytical();
	    }
#endif
	}

      if (Max_ln_p_value > -INFINITY)
	{
	  L_thresh = cutoff_approx(Max_ln_p_value, " set by the user");
	  if (L_thresh <= -INFINITY)
	    {
	      fprintf(stderr,
		 "Try decreasing the minimum score for calculating p-values\n");
	      fprintf(stderr, "using the \"-M\" command line option.\n");
	      exit(1);
	    }
	  if (Thresh == '\0') Thresh = 'l';
	  else if (Thresh == 'u') Thresh = 'b';
	  else bug_report("main() 2");
	}
    }
  else Min_exact_score = -INFINITY / 2.0;

  /* Initialize information for "read_sequence()". */
  init_read_seq();

  /* Initialize information for "process_sequence()". */
  init_variables();

  /* Read a sequence and determine whether it is circular.
   * Returns YES, if a sequence exists. */
  while (read_sequence() == YES)
    /* Determine the score of the L-mers of the current sequence. */
    process_sequence();

  exit(0);
}
