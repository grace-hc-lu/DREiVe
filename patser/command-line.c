/* Copyright 1990, 1994, 1995, 1996, 1997, 2000, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


#include "options.h"
#include "parse-line.h"

#define MESSAGE "Use the \"-h\" option for detailed directions"

/* This file contains functions for reading the command line.
 * The following external variables are used from "options.c":
 *
 *    Mat_file            L_thresh            M_fp
 *    Seq_file    	  U_thresh 	      S_fp
 *    Comp_status 	  Top      	      Sum_mat
 *    Thresh      	  Error    	      Vertical
 *    Fudge               Case_sensitive      Max_int_range
 *    Test                Min_score           Skip
 *    Max_ln_p_value      Auto_cutoff         Equal_err
 */


/* Functions specific to this particular program and
 * needed in the OPTION vector. */
extern int pl_Lower();   /* Determine lower threshhold and mark "Thresh". */
extern int pl_Upper();   /* Determine upper threshhold and mark "Thresh". */
extern int pl_Print_top(); /* Determine number of top scores to print for each
			      sequence and assign to "Top".
			      Assume "Top" equals 1 if not explicitly given. */


/* Vector holding the possible command line options. */
static OPTION ST_options[] =
{
  "h", 	 OPT, pl_Help,      NULL,            "Print directions",
  "m", 	 OPT, pl_String,    &Mat_file,       "Name of matrix file---default name is \"matrix\"",
  "w", 	 OPT, pl_NBool,     &Sum_mat,        "Matrix is a weight matrix",
  "v", 	 OPT, pl_Bool,      &Vertical,       "Vertical matrix---rows correspond to positions",
  "f", 	 OPT, pl_String,    &Seq_file,       "Name of sequence file---default: standard input",
  "b", 	 OPT, pl_Nn_Double, &Fudge,          "Correction added to the elements of the alignment matrix (default: 1)",
  "a", 	 OPT, pl_Alpha_af,  NULL,            "Name of ascii alphabet file---default name is \"alphabet\"",
  "i", 	 OPT, pl_Alpha_if,  NULL,            "Name of integer alphabet file",
  "A", 	 OPT, pl_Alpha_ac,  NULL,            "Ascii alphabet information",
  "CS",  OPT, pl_NBool,     &Case_sensitive, "Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)",
  "CM",  OPT, pl_Int_2,     &Case_sensitive, "Ascii alphabet is case insensitive, but mark the location of lowercase letters",
  "c", 	 OPT, pl_Bool,      &Comp_status,    "Score the complementary strand",
  "ls",  OPT, pl_Lower,     NULL,            "Lower-threshold score, inclusive (formerly the -l option)",
  "li",  OPT, pl_Bool,      &Auto_cutoff,    "Determine lower-threshold score from adjusted information content",
  "lp",  OPT, pl_Np_Double, &Max_ln_p_value, "Determine lower-threshold score from a maximum ln(p-value)",
  "u",   OPT, pl_Upper,     NULL,            "Upper-threshold score, exclusive",
  "t", 	 OPT, pl_Print_top, &Top,            "Print only the top scores",
  "ds",  OPT, pl_Bool,      &Print_order,    "Print top scores in order of decreasing score (default: print in order of position)",
  "e", 	 OPT, pl_Nn_Double, &Equal_err,      "Small difference for considering 2 scores equal (default: 0.000001)",
  "d0",  OPT, pl_NBool,     &Error,          "Treat unrecognized characters as errors",
  "d1",  OPT, pl_Bool,      &Error,          "Treat unrecognized characters as discontinuities, but print warning (the default)",
  "d2",  OPT, pl_Int_2,     &Error,          "Treat unrecognized characters as discontinuities, and print NO warning",
  "p", 	 OPT, pl_Int_2,     &Print,          "Vertically print the weight matrix",
  "R", 	 OPT, pl_Nn_Double, &Max_int_range,  "Set the range for approximating a weight matrix with integers (default: 10000)",
  "M", 	 OPT, pl_Double,    &Min_score,      "Set the minimum score for calculating the p-value of scores (default: 0)",
#ifdef TEST
  "Ts",  OPT, pl_P_Int,     &Skip,      "Number of positions to skip in the marginal probability array when calculating p-values analytically (default: 1)",
  "T1",  OPT, pl_Bool,      &Test,      "Determine the cutoff score analytically based on the information content after adding the small sample size correction",
  "T2",  OPT, pl_Int_2,     &Test,      "Determine the cutoff score exactly and analytically based on the information content after adding the small sample size correction",
#endif
  END
};


/* Vector indicating mutually exclusive options. */
static char *ST_exclusive[] =
{
  "a", "i", "A", END,
  "CS", "CM", "i", END,
  "w", "b", END,
  "ls", "li", "lp", END,
  "d0", "d1", "d2", END,
#ifdef TEST
  "w", "T1", END,
  "w", "T2", END,
#endif
  END
};


/* Read the command line options.  Print to standard output the PID
 * and the options chosen. */
void command_line(argc, argv)
     int argc;
     char *argv[];
{
  void adjust_alphabet();
  void check_variables();
  void print_options();

  /* Read the command line. */
  if (parse_line(ST_options, ST_exclusive, argc, argv) == NO)
    usage(ST_options, MESSAGE, argv);

  /* Read alphabet from indicated file, if not already read from command line.
   * Make sure the alphabet is complementary if both strands are being used.
   * Convert to uppercase letters if (Case_sensitive != 0).
   * Make sure none of the characters in the alphabet occur more than once.
   * Convert the alphabet normalizations to frequencies. */
  adjust_alphabet();

  /* Check the validity of the command line variables, and determine the
   * variables that are dependent on the command line variables. */
  check_variables();

  /* Print the options chosen. */
  print_options(argc, argv);
}


/* Command specific functions.  Returns the index of the next
 * command line variable.  Return values <= 0 indicate an error.
 * The following is its format: /
int function(variable, argc, argv, i, k)
     char *variable;     / Address of the variable to be updated. /
     int argc;           / Number of variables on the command line. /
     char *argv[];       / The table of command line strings. /
     int i;              / The index to the current command line string. /
     int k;              / Index to current position of the current string. /
 */


/* Functions specific to this particular program. /
int pl_Lower();        / Determine lower threshhold and mark "Thresh". /
int pl_Upper();        / Determine upper threshhold and mark "Thresh". /
int pl_Print_top();    / Determine number of top scores to print for each
			 sequence and assign to "Top".
			 Assume "Top" equals 1 if not explicitly given. /
 */

/* Determine the lower threshhold and mark "Thresh". */
int pl_Lower(variable, argc, argv, i, k)
     char *variable;     /* Address of the variable to be updated. */
     int argc;           /* Number of variables on the command line. */
     char *argv[];       /* The table of command line strings. */
     int i;              /* The index to the current command line string. */
     int k;              /* Index to current position of the current string. */
{
  if (Thresh == '\0') Thresh = 'l';
  else if (Thresh == 'u') Thresh = 'b';
  else bug_report("pl_Lower()");

  return(pl_Double(&L_thresh, argc, argv, i, k));
}


/* Determine the upper threshhold and mark "Thresh". */
int pl_Upper(variable, argc, argv, i, k)
     char *variable;     /* Address of the variable to be updated. */
     int argc;           /* Number of variables on the command line. */
     char *argv[];       /* The table of command line strings. */
     int i;              /* The index to the current command line string. */
     int k;              /* Index to current position of the current string. */
{
  if (Thresh == '\0') Thresh = 'u';
  else if (Thresh == 'l') Thresh = 'b';
  else bug_report("pl_Upper()");

  return(pl_Double(&U_thresh, argc, argv, i, k));
}


/* Determine number of top scores to print for each sequence and
   assign to "Top".  Assume "Top" equals 1 if not explicitly given;
   otherwise, must be a positive number. */
int pl_Print_top(
     void *variable,     /* Address of the variable to be updated. */
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int i,              /* The index to the current command line string. */
     int k)              /* Index to current position of the current string. */
{
  char c = '\0';

  if ((i >= argc) || ((isdigit(argv[i][k]) == NO) && (k == 0)))
    {
      *(int *)variable = 1;
      return(i);
    }

  else if ((sscanf(&argv[i][k], "%d%c", (int *)variable, &c) != 1)
	   || (c != '\0') || (*((int *)variable) <= 0))
    {
      fprintf(stderr, "\"%s\" (from item %d on the command line) ",
	      &argv[i][k], i);
      fprintf(stderr, "is not a positive integer.\n\n");
      return(NO);
    }

  else
    return(i + 1);
}


/* Check the validity of the command line variables, and determine the
 * variables that are dependent on the command line variables. */
void check_variables()
{
  int status = YES;       /* Validity of default values in "options.c" file. */


  /* Open the matrix file. */
  if ((M_fp = fopen(Mat_file, "r")) == (FILE *)NULL)
    {
      fprintf(stderr, "Matrix file \"%s\" cannot be opened.\n", Mat_file);
      status = NO;
    }

  /* Open the file containing the names of the sequence files. */
  if (Seq_file == (char *)NULL) S_fp = stdin;
  else if ((S_fp = fopen(Seq_file, "r")) == (FILE *)NULL)
    {
      fprintf(stderr, "File of sequences \"%s\" cannot be opened.\n",Seq_file);
      status = NO;
    }

  /* Make sure the upper threshold is greater than the lower threshold. */
  if ((Thresh == 'b') && (U_thresh <= L_thresh))
    {
      fprintf(stderr, "Upper threshold for printing scores must be ");
      fprintf(stderr, "greater than lower threshold.\n");
      status = NO;
    }


  /* Exit if any of the values are not valid. */
  if (status == NO) exit(1);
}


/* Print the options chosen. */
void print_options(argc, argv)
     int argc;
     char *argv[];
{
  int i;
  void print_alpha();


  /* Echo the command line. */
  printf("COMMAND LINE:");
  for (i = 0; i < argc; ++i) printf(" %s", argv[i]);
  printf("\n\n");

  printf("File containing the matrix: %s\n", Mat_file);

  if (Seq_file == (char *)NULL)
    printf("Sequence information from the standard input\n");
  else printf("File containing the sequence information: %s\n", Seq_file);

  /* Print: type of matrix and Fudge factor (if applicable). */
  printf("Type of matrix: %s\n", (Sum_mat == YES) ? "alignment" : "weight");
  if (Sum_mat == YES)
    {
      printf("Total pseudo-counts added to the elements ");
      printf("of the alignment matrix: %.3g\n", Fudge);
    }
  printf("\n");

  /* Print: the range for approximating a weight matrix with integers, and
   * the minimum score for calculating the p-values of scores. */
  if (Max_int_range <= 0.0)
    printf("Do not numerically approximate the p-values of the scores.\n");
  else
    {
      printf("Range for approximating a weight matrix with integers: %d\n",
	     ROUND_TO_INT(Max_int_range));
      printf("Minimum score for calculating p-values: %.3g\n", Min_score);
    }

  /* Print: status of complementary strand and unrecognized characters. */
  if (Comp_status == YES)
    printf("Also score the complementary strands\n");

  if (Error == 0)
    printf("Treat unrecognized symbols as errors.\n");
  else if ((Error == 1) || (Error == 2))
    printf("Treat unrecognized symbols as discontinuities in the sequence.\n");
  else bug_report("print_options() 1");

  /* Print the adjustments and restrictions on which scores are printed. */
  if (Top == 0)
    {
      if ((Thresh == '\0') &&
	  (Max_ln_p_value <= -INFINITY) && (Auto_cutoff == NO))
	printf("Print ALL the scores\n");
      else if (Thresh == 'l')
	printf("Print scores greater than or equal to %.2f\n", L_thresh);
      else if (Thresh == 'u')
	printf("Print scores less than %.2f\n", U_thresh);
      else if (Thresh == 'b')
	{
	  printf("Print scores greater than or equal to ");
	  printf("%.2f and less than %.2f\n", L_thresh, U_thresh);
	}
    }
  else
    {
      if (Thresh == '\0')
	printf("Print %d top score(s) for each sequence\n", Top);
      else if (Thresh == 'l')
	{
	  printf("Print %d top score(s) for each sequence, if greater ", Top);
	  printf("than or equal to %.2f\n", L_thresh);
	}
      else if (Thresh == 'u')
	{
	  printf("Print %d top score(s) for each sequence, if less ", Top);
	  printf("than %.2f\n", U_thresh);
	}
      else if (Thresh == 'b')
	{
	  printf("Print %d top score(s) for each sequence, if between ", Top);
	  printf("%.2f and %.2f\n", L_thresh, U_thresh);
	}
      else bug_report("print_options() 2");
    }
  if (Max_ln_p_value > -INFINITY)
    {
      printf("Print scores having a ln(p-value) less than or equal ");
      printf("to %.2f\n", Max_ln_p_value);
    }
  else if (Auto_cutoff == YES)
    {
      printf("Print scores having a ln(p-value) less than or equal\n");
      printf("to the negative sample-size adjusted information content\n");
    }
  printf("\n");

  /* Print the alphabet information. */
  print_alpha();
  printf("\n");
}
