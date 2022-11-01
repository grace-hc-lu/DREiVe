/* Copyright 1990, 1994, 1995, 1996 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   Dept. of Molecular, Cellular, and Developmental Biology
 *   University of Colorado
 *   Campus Box 347
 *   Boulder, CO  80309-0347
 *
 *   hertz@boulder.colorado.edu
 */


#include "options.h"


/* Convert the alignment matrix to a weight matrix.
////////////////////////////////////////////////////////////////////////////
// MORE INFORMATION ON WEIGHT MATRICES                                    //
// Weight matrices are used to evaluate other sequences according to      //
// how closesly they match the summarized sequences.                      //
////////////////////////////////////////////////////////////////////////////
// The following procedure is used to convert the alignment matrix into a //
// weight matrix:                                                         //
****************************************************************************
** add (P[i] * Fudge) to each matrix-element;                             **
** divide each element by ((total # of sequences) + Fudge);               **
** normalize above ratio to the prior frequency of corresponding letter;  **
** take the logarithm, to the base BASE, of each fraction.                **
***************************************************************************/

/* Global external variables used:
    Mat_file            P                   L_1
    Ascii       	M_fp    	    Letter_array
    Comp_flag   	Vertical	    Letter_index
    A_size      	Matrix              Align_mat
*/


/* Read the alignment matrix and derive the weight matrix. */
void process_matrix()
{
  int i;
  double **read_h_align_mat(); /* Read the horizontal matrix of numbers. */
  double **read_v_align_mat(); /* Read the vertical matrix of numbers. */
  double **weight_matrix();  /* Determine the weight matrix. */
  double sum_column();       /* Determine the sum of the matrix columns. */
  double inf_mat();    /* Determine information content of alignment matrix.*/
  double adjust_inf(); /* Determine adjusted information of alignment matrix.*/
  void fudge_align_mat();    /* Add (P[i] * Fudge) to each matrix-element. */
  void print_h_align_mat();  /* Print horizontal alignment matrix and width. */
  void print_v_align_mat();  /* Print vertical alignment matrix and width. */
  void print_h_matrix();     /* Print the horizontal weight matrix. */
  void print_v_matrix();     /* Print the vertical weight matrix. */


  /* Reads the contents of the matrix file into a dynamic matrix, determine "L"
   * and "L_1", make sure the matrix is not empty, and close matrix file. */
  if (Vertical == NO)
    Align_mat = read_h_align_mat();        /* Read horizontal matrix. */
  else
    Align_mat = read_v_align_mat();        /* Read a vertical matrix. */

#if 0
  /* Determine the sum of the matrix columns and make sure all
   * columns have the same sum. */
  column_sum = sum_column(Align_mat);
#endif

  /* Print the alignment matrix and its width. */
  if (Vertical == NO) print_h_align_mat(Align_mat);
  else print_v_align_mat(Align_mat);

  /* Determine the information content of the alignment matrix. */
  Inf = inf_mat(Align_mat);
  /* Determine the adjusted information content of the alignment matrix. */
  Inf_adj = adjust_inf(Align_mat, Inf);

  /* Add (P[i] * Fudge) to each matrix-element. */
  fudge_align_mat(Align_mat);
  /* Determine the information content of fudged alignment matrix. */
  Inf_fudge = inf_mat(Align_mat);

  /* Determine the weight matrix and save it in the "Matrix[][]". */
  Matrix = weight_matrix(Align_mat);

#if 0
  /* Print the alignment matrix and its width. */
  if (Vertical == NO) print_h_align_mat(Align_mat);
  else print_v_align_mat(Align_mat);

  /* Free the memory allocated for the alignment matrix. */
  for (i = 0; i < A_size; ++i)
    free_error((char *)Align_mat[i], "Align_mat[i]", "process_matrix()");
  free_error((char *)Align_mat, "Align_mat", "process_matrix()");
#endif


  /* Print the weight matrix. */
  if (Print == 1) print_h_matrix(Matrix, L, Mat_file);
  else if (Print == 2) print_v_matrix(Matrix, L, Mat_file);

  /* Print the information contents. */
  printf("Information content (base e): %7.3f\n", Inf);
  printf("Sample size adjusted information content\n     (information ");
  printf("content minus the average information\n     expected from an ");
  printf("arbitrary alignment of random sequences): %7.3f\n", Inf_adj);
  printf("Information content after adding pseudo-counts: %7.3f\n", Inf_fudge);
  printf("\n");
}


/* Read a horizontal matrix.
 * Reads the contents of the matrix file into a dynamic matrix, determine "L"
 * and "L_1", make sure the matrix is not empty, and close matrix file. */
double **read_h_align_mat()
{
  int i;
  void read_first_row();  /* Read the first row of the matrix. */
  void read_h_row();      /* Read a row of the matrix, except for first row. */
  double **align_mat;     /* The alignment matrix in the matrix file. */
  long find_line();       /* Find the beginning of a data line. */


  /* Allocate space for "Letter_array". */
  Letter_array = (char *)calloc_error(A_size, sizeof(char),
				      "Letter_array", "read_h_align_mat");
  for (i = 0; i < A_size; ++i) Letter_array[i] = NO;

  /* Allocate space for "align_mat[][]". */
  align_mat = (double **)calloc_error(A_size, sizeof(double *),
				      "align_mat", "read_h_align_mat");

  /* Read the first row of the matrix and determinine L. */
  read_first_row(align_mat);

  /* Determine the rest of the rows of the matrix. */
  for (i = 2; i <= A_size; ++i) read_h_row(align_mat, i);

  /* Make sure there are no illegal characters at the end of the file. */
  if (find_line(M_fp, Mat_file) != -1L)
    {
      fprintf(stderr, "Extra characters at the end of the matrix file.\n");
      exit(1);
    }

  L_1 = L - 1;        /* Determine "L_1". */

  /* Close the matrix file and free space for "Letter_array[]". */
  fclose(M_fp);
  free_error(Letter_array, "Letter_array", "read_h_align_mat()");

  return(align_mat);
}

/* Read a vertical matrix.
 * Reads the contents of the matrix file into a dynamic matrix, determine "L"
 * and "L_1", make sure the matrix is not empty, and close matrix file. */
double **read_v_align_mat()
{
  int i;
  void read_letter_row();    /* Read a row of letters. */
  char read_v_row();         /* Read a row of integers from vertical matrix. */
  double **align_mat;          /* The alignment matrix in the matrix file. */


  /* Allocate space for "Letter_array" and "Letter_index". */
  Letter_array = calloc_error(A_size, sizeof(char),
			      "Letter_array", "read_v_align_mat");
  for (i = 0; i < A_size; ++i) Letter_array[i] = NO;
  Letter_index = (int *)calloc_error(A_size, sizeof(int),
			      "Letter_index", "read_v_align_mat");

  /* Allocate space for "align_mat[][]". */
  align_mat = (double **)calloc_error(A_size, sizeof(double *),
				      "align_mat", "read_v_align_mat");
  for (i = 0; i < A_size; ++i) align_mat[i] = (double *)NULL;

  /* Read the letters in the first row of the matrix. */
  read_letter_row();

  /* Read the rows of integers and increment "L". */
  L = 0;
  while (read_v_row(align_mat) == YES) ++L;

  if (L == 0)
    {
      fprintf(stderr, "The first row of the alignment matrix ");
      fprintf(stderr, "contains no integers.\n\n");
      exit(1);
    }

  /* Truncate the rows of "align_mat[][]". */
  for (i = 0; i < A_size; ++i)
    align_mat[i] = (double *)recalloc_error((char *)align_mat[i],
			L, sizeof(double), "align_mat[]","read_v_align_mat()");

  L_1 = L - 1;        /* Determine "L_1". */

  /* Close matrix file and free space for "Letter_array" and "Letter_index". */
  fclose(M_fp);
  free_error(Letter_array, "Letter_array", "read_v_align_mat()");
  free_error((char *)Letter_index, "Letter_index", "read_v_align_mat()");

  return(align_mat);
}


/* Read the first row of the matrix and determinine L. */
void read_first_row(align_mat)
     double **align_mat;
{
  int letter;                /* Letter index corresponding to current row. */
  int array_size= CHUNK_SIZE;/* The size of the array holding the first row. */
  int read_letter();         /* Read the letter index corresponding to row.
			      * Comment lines beginning with ';', '%', or '#'
			      * are ignored along with the '|'. */
  double read_pos_double();  /* Read a non-negative double from the matrix
			      * file.  The '|' character is ignored. */


  /* Determine letter index of first row. */
  if ((letter = read_letter()) == EOF)
    {
      fprintf(stderr, "Cannot find the letter corresponding to the ");
      fprintf(stderr, "first row of the alignment matrix.\n");
      exit(1);
    }

  /* Allocate initial space for the first row. */
  align_mat[letter] = (double *)calloc_error(array_size, sizeof(double),
					     "align_mat[]", "read_first_row");

  /* Read the first row of the alignment matrix. */
  while ((align_mat[letter][L] = read_pos_double(M_fp, 1)) >= 0.0)
    {
      ++L;

      if (L >= array_size)
	{
	  array_size = array_size + CHUNK_SIZE;
	  align_mat[letter] =
	    (double *)recalloc_error((char *)align_mat[letter],
		  array_size, sizeof(double), "align_mat[]", "read_first_row");
	}
    }

  /* Truncate the row, if the matrix is not empty. */
  if (L == 0)
    {
      fprintf(stderr, "The first row of the alignment matrix ");
      fprintf(stderr, "contains no non-negative numbers.\n\n");
      exit(1);
    }
  else
    align_mat[letter] = (double *)recalloc_error((char *)align_mat[letter],
			   L, sizeof(double), "align_mat[]", "read_first_row");
}

/* Read a row of the horizontal matrix, except for first row. */
void read_h_row(align_mat, row)
     double **align_mat;  /* The alignment matrix. */
     int row;             /* The row currently being read. */
{
  int i;
  int letter;
  int read_letter();         /* Read the letter index corresponding to row.
			      * Comment lines beginning with ';', '%', or '#'
			      * are ignored along with the '|'. */
  double read_pos_double();  /* Read a non-negative double from the matrix
			      * file.  The '|' character is ignored. */


  /* Determine letter index of the current row. */
  if ((letter = read_letter()) == EOF)
    {
      fprintf(stderr, "Cannot find the letter corresponding to ");
      fprintf(stderr, "row %d of the alignment matrix.\n", row);
      exit(1);
    }

  /* Allocate space for the current row. */
  align_mat[letter] = (double *)calloc_error(L, sizeof(double),
					     "align_mat[]", "read_h_row");

  for (i = 0; i < L; ++i)
    {
      if ((align_mat[letter][i] = read_pos_double(M_fp, row)) < 0.0)
	{
	  fprintf(stderr, "Row %d of the alignment matrix has ", row);
	  fprintf(stderr, "less elements than row 1.\n");
	  exit(1);
	}
    }

  /* Make sure there are no extra numbers or other inappripriate
   * characters at the end of the line. */
  if (read_pos_double(M_fp, row) >= 0.0)
    {
      fprintf(stderr, "Row %d of the alignment matrix has ", row);
      fprintf(stderr, "more elements than row 1.\n");
      exit(1);
    }
}

/* Read the letter index corresponding to the current row.  Comment lines
 * beginning with ';', '%', or '#' are ignored along with the '|'.
 * Returns the letter index or EOF, if a letter is not found. */
int read_letter()
{
  int letter;              /* The letter corresponding to the current row. */
  int letter_index;        /* The index of the current letter. */
  int translate_letter();  /* Translate letter into integral representation. */
  long find_line();        /* Find the beginning of a data line. */


  if (Ascii == YES)
    {
      while (1)
	{
	  if (find_line(M_fp, Mat_file) == -1L) return(EOF);
	  letter = getc(M_fp);
	  if (letter == '|');
	  else
	    {
	      if ((Case_sensitive != 0) && (islower(letter) != NO))
		letter = toupper(letter);
	      letter_index = translate_letter(letter);
	      break;
	    }
	}
    }
  else
    {
      while (1)
	{
	  if (find_line(M_fp, Mat_file) == -1L) return(EOF);
	  if (fscanf(M_fp, "%d", &letter) == 1)
	    {
	      letter_index = translate_letter(letter);
	      break;
	    }
	  else if (getc(M_fp) == '|');
	  else return(EOF);
	}
    }

  /* Make sure the letter does not occur on more than one row of the matrix. */
  if (Letter_array[letter_index] == NO)
    {
      Letter_array[letter_index] = YES;
      return(letter_index);
    }
  else
    {
      if (Ascii == YES) fprintf(stderr, "The letter \"%c\" ", letter);
      else fprintf(stderr, "The integer letter \"%d\" occurs ", letter);

      fprintf(stderr, "on more than one row of the alignment matrix.\n");
      exit(1);
    }
}

/* Translate a letter into its integral representation. */
int translate_letter(letter)
     int letter;
{
  int i;

  /* Find the index of the "letter". */
  for (i = 0; i < A_size; ++i)
    {
      if (letter == A[i])
	return(i);
    }

  /* Error message if the "letter" does not match
   * any of the letters in the alphabet. */
  if (Ascii == YES)
    fprintf(stderr, "The letter \"%c\" ", letter);
  else
    fprintf(stderr, "The integer symbol \"%d\" ", letter);
  fprintf(stderr, "in the matrix file is not a member of the alphabet.\n");
  exit(1);
}


/* Read a non-negative double from the matrix information.  The '|' character
 * is ignored.  Returns the non-negative number or -1, if a non-negative
 * number is not found by the end of the current line.
 * If a non-negative number is not found, will either scan to the end of the 
 * line or print an error message and "exit(1)" if an inappropriate character
 * is found (i.e., not spaces, '|', or comments */
double read_pos_double(fp, row)
     FILE *fp;
     int row;     /* The row currently being read. */
{
  double number;     /* The return value of this function. */
  int letter;        /* A letter in the matrix file. */


  while ((letter = getc(fp)) != EOF)
    {
      /* Return if end of line; pass over any other whitespace or '|'. */
      if (letter == '\n') return(-1.0);
      else if (isspace(letter) != NO);
      else if (letter == '|');

      /* Read the integer. */
      else if ((isdigit(letter) != NO) || (letter == '.') || (letter == '+')) 
	{
	  if ((ungetc(letter, fp) == EOF) || (fscanf(fp, "%lf", &number) != 1))
	    {
	      fprintf(stderr, "Problems reading a non-negative number from ");
	      fprintf(stderr, "the matrix information.\n");
	      exit(1);
	    }
	  else return(number);
	}

      /* Treat a comment as if it is the end of the line. */
      else if ((letter == ';') || (letter == '#') || (letter == '%'))
	{
	  while ((letter = getc(fp)) != EOF) if (letter == '\n') return(-1.0);
	}

      /* An error if the character is inappropriate. */
      else
	{
	  fprintf(stderr, "The character \"%c\" occurs inappropriat", letter);
	  fprintf(stderr, "ely in row %d of the matrix information.\n", row);
	  exit(1);
	}
    }

  /* Return -1 if no more letters. */
  if (feof(fp) == NO)
    {
      fprintf(stderr, "An error occurred while reading the matrix file.\n");
      exit(1);
    }
  else return(-1.0);
}

/* Read the letters in the first row.  Fill in "Letter_index[]"
 * with the corresponding letter indices. */
void read_letter_row()     /* Read a row of letters. */
{
  int letter;              /* The letter just read from the matrix file. */
  int letter_index;        /* The index of the letter just read. */
  int i;
  void scan_end();         /* Scan to the end of the row. */
  int translate_letter();  /* Translate letter into integral representation. */
  long find_line();        /* Find the beginning of a data line. */


  /* Find the first data line in the matrix file. */
  if (find_line(M_fp, Mat_file) == -1L)
    {
      fprintf(stderr, "The matrix file does not contain a matrix.\n");
      exit(1);
    }

  /* Read letters when the alphabet contains ASCII characters. */
  if (Ascii == YES)
    {
      for (i = 0; i < A_size; ++i)
	{
	  letter = getc(M_fp);
	  if ((letter == EOF) || (letter == '\n') ||
	      (letter == ';') || (letter == '%') || (letter == '#'))
	    {
	      fprintf(stderr, "The vertical matrix file is missing ");
	      fprintf(stderr, "letters in its first row.\n");
	      exit(1);
	    }
	  else if ((isspace(letter) != NO) || (letter == '|'))
	    --i;
	  else
	    {
	      if ((Case_sensitive != 0) && (islower(letter) != NO))
		letter = toupper(letter);
	      letter_index = translate_letter(letter);

	      /* Make sure letter does not occur in more than one column. */
	      if (Letter_array[letter_index] == NO)
		{
		  Letter_array[letter_index] = YES;
		  Letter_index[i] = letter_index;
		}
	      else
		{
		  fprintf(stderr, "The letter \"%c\" occurs in more ", letter);
		  fprintf(stderr,"than one column of the alignment matrix.\n");
		  exit(1);
		}
	    }
	}
    }

  /* Read letters when the alphabet contains integer characters. */
  else
    {
      for (i = 0; i < A_size; ++i)
	{
	  letter = getc(M_fp);
	  if ((letter == EOF) || (letter == '\n') ||
	      (letter == ';') || (letter == '%') || (letter == '#'))
	    {
	      fprintf(stderr, "The vertical matrix file is missing ");
	      fprintf(stderr, "letters in its first row.\n");
	      exit(1);
	    }
	  else if ((isspace(letter) != NO) || (letter == '|'))
	    --i;
	  else
	    {
	      if (ungetc(letter, M_fp) == EOF)
		{
		  fprintf(stderr, "\"ungetc\" error\n");
		  exit(1);
		}

	      if (fscanf(M_fp, "%d", &letter) == 1)
		letter_index = translate_letter(letter);
	      else
		{
		  letter = getc(M_fp);

		  fprintf(stderr, "The character \"%c\" is inapp", letter);
		  fprintf(stderr, "ropriate for the first line\nof a ");
		  fprintf(stderr, "vertical matrix that has an integer ");
		  fprintf(stderr, "alphabet.\n");
		  exit(1);
		}
		
	      /* Make sure letter does not occur in more than one column. */
	      if (Letter_array[letter_index] == NO)
		{
		  Letter_array[letter_index] = YES;
		  Letter_index[i] = letter_index;
		}
	      else
		{
		  fprintf(stderr, "The integer letter \"%d\" ", letter);
		  fprintf(stderr, "occurs in more than one column of the ");
		  fprintf(stderr, "alignment matrix.\n");
		  exit(1);
		}
	    }
	}
    }
  scan_end("the letter row");
}

/* Read a row of integers from the vertical matrix.
 * Returns YES if the row is found; otherwise, NO. */
char read_v_row(align_mat)
     double **align_mat;
{
  int i;
  static int array_size = 0; /* The size of the arrays holding the rows. */
  double read_pos_double();  /* Read a non-negative double from the matrix
			      * file.  The '|' character is ignored. */
  long find_line();          /* Find the beginning of a data line. */
  

  if (find_line(M_fp, Mat_file) == -1L) return(NO);

  /* Make sure the alignment matrix has enough space in each row. */
  if (L >= array_size)
    {
      array_size = array_size + CHUNK_SIZE;
      for (i = 0; i < A_size; ++i)
	align_mat[i] = (double *)recalloc_error((char *)align_mat[i],
		      array_size, sizeof(double), "align_mat[]", "read_v_row");
    }

  /* Read a row of the alignment matrix. */
  for (i = 0; i < A_size; ++i)
    {
      if ((align_mat[Letter_index[i]][L] = read_pos_double(M_fp, L + 1)) < 0.0)
	{
	  if (i == 0)
	    {
	      --i;
	      if (find_line(M_fp, Mat_file) == -1L)
		return(NO);
	    }
	  else
	    {
	      fprintf(stderr, "Row %d of the alignment matrix contains ", L+1);
	      fprintf(stderr, "too few non-negative numbers.\n");
	      exit(1);
	    }
	}
    }

  /* Make sure there are no extra non-negative doubles or other inappripriate
   * characters at the end of the line. */
  if (read_pos_double(M_fp, L + 1) >= 0.0)
    {
      fprintf(stderr, "Row %d of the matrix contains more ", L + 1);
      fprintf(stderr, "elements than the number of columns.\n");
      exit(1);
    }

  return(YES);
}

/* Scan to the end of the row.  Ignore '|', spaces, and comments. */
void scan_end(row)
     char *row;    /* Message fitting into the error message below
		    * to describe which row is being scanned. */
{
  int letter;

  letter = getc(M_fp);
  while ((letter != EOF) && (letter != '\n'))
    {
      if ((isspace(letter) != NO) || (letter == '|'))
	letter = getc(M_fp);
      else if ((letter == ';') || (letter == '#') || (letter == '%'))
	while (((letter = getc(M_fp)) != EOF) && (letter != '\n'));
      else
	{
	  fprintf(stderr, "The character \"%c\" at the end of ", letter);
	  fprintf(stderr, "%s of the matrix is not a valid character.\n", row);
	  exit(1);
	}
    }
}


#if 0
/* Determine the sum of the matrix columns and make sure all
 * columns have the same sum. */
double sum_column(align_mat)
     double **align_mat; /* The alignment matrix. */
{
  int i, j;
  double column_sum;        /* The sum of one of the matrix colums. */
  double temp_sum;          /* A temporary sum of one of the matrix colums. */
  int status = NO;        /* NO: no errors.  YES: some errors. */


  /* Determine the sum of the first column. */
  for (i = 0, column_sum = 0.0; i < A_size; ++i) column_sum += align_mat[i][0];

  /* Compare the sum of the first column with the sums of the other columns.*/
  for (j = 1; j < L; ++j)
    {
      /* Determine the sum of the column. */
      for (i = 0, temp_sum = 0.0; i < A_size; ++i) temp_sum += align_mat[i][j];

      /* Compare the sum of first column with the sum of current column. */
      if (fabs(temp_sum - column_sum) / column_sum > ERROR)
        {
          status = YES;

          fprintf(stderr, "The number of sequences represented at position ");
	  fprintf(stderr, "%d of the alignment matrix\ndoes not match the ",
		  j + 1);
          fprintf(stderr, "number of sequences at the first position.\n\n");
        }
    }

  if (status == YES) exit(1);
  return(column_sum);
}
#endif


/* Determine the information content of the alignment matrix.*/
double inf_mat(align_mat)
     double **align_mat;   /* The alignment matrix. */
{
  int i, j;
  double col_inf;          /* The information content of the current column. */
  double col_sum;          /* The sum of the current column. */
  double inf;              /* The information content. */


  /* Manipulate each position of each row. */
  for (j = 0, inf = 0.0; j < L; ++j)
    {
      for (i = 0, col_inf = col_sum = 0.0; i < A_size; ++i)
	{
	  if (align_mat[i][j] > 0.0)
	    {
	      col_inf += align_mat[i][j] * log(align_mat[i][j] / P[i]);
	      col_sum += align_mat[i][j];
	    }
	}
      if (col_sum > 0.0) col_inf = col_inf / col_sum - log(col_sum);
      inf += col_inf;
    }

  return(inf);
}
	  
/* Determine the adjusted information content of the alignment matrix. */
double adjust_inf(align_mat, inf)
     double **align_mat;   /* The alignment matrix. */
     double inf;           /* The information content. */
{
  int i, j, k, l;
  double inf_adj;   /* Adjusted information content of the alignment matrix. */
  double col_sum;          /* Sum of the current column. */
  int col_sum_int;         /* Column sum rounded to the closes int. */
  int max_col_sum;         /* The maximum rounded column sum. */
  double *ln_fact;         /* A list of ln(n!). */
  int fact_num;            /* The number factorials determined. */
  double avg_inf();        /* Determine the average information content. */

  /* Array holding background information in case columns
   * contain differing sums. */
  double *background = (double *)NULL;


  fact_num = 1;
  ln_fact = (double *)calloc_error(fact_num + 1, sizeof(double),
				   "ln_fact", "adjust_inf()");
  ln_fact[0] = ln_fact[1] = 0.0;

  for (j = 0, max_col_sum = -1, inf_adj = inf; j < L; ++j)
    {
      /* Determine the column sum. */
      for (i = 0, col_sum = 0.0; i < A_size; ++i) col_sum += align_mat[i][j];
      col_sum_int = (int)(col_sum + 0.5);

      /* Allocate additional space for ln_fact[]. */
      if (col_sum_int > fact_num) 
	{
	  ln_fact = (double *)recalloc_error((char *)ln_fact,
		   col_sum_int + 1, sizeof(double), "ln_fact", "adjust_inf()");
	  for (k = fact_num + 1, l = fact_num; k <= col_sum_int; l = k++)
	    ln_fact[k] = log((double)k) + ln_fact[l];
	  fact_num = col_sum_int;
	}

      /* Allocate additional space for background[]. */
      if (col_sum_int > max_col_sum) 
	{
	  background = (double *)recalloc_error((char *)background,
		col_sum_int + 1, sizeof(double), "background", "adjust_inf()");
	  for (k = max_col_sum + 1; k <= col_sum_int; ++k)background[k] = -1.0;
	  max_col_sum = col_sum_int;
	}

      if (background[col_sum_int] < 0.0)
	background[col_sum_int] = avg_inf(col_sum_int, ln_fact);

      inf_adj -= background[col_sum_int];
    }

  free_error((char *)ln_fact, "ln_fact", "adjust_inf()");
  free_error((char *)background, "background", "adjust_inf()");
  return(inf_adj);
}

/* Determine the average information content. */
double avg_inf(n, ln_fact)
     int n;                     /* The number of sequences in the alignment. */
     double *ln_fact;           /* A list of ln(n!). */
{
  int i, j, k;
  double dj, dk;          /* (double)j and (double)k */
  double dn = (double)n;  /* (double)n */
  double avg;             /* The average information content. */
  double ln_combi;        /* ln[N!/n!(N-n)!] for components of distribution. */
  double prob;            /* Binomial probability for observed information. */
  double ln_p;            /* ln(P[]). */
  double ln_1_p;          /* ln(1 - P[]). */


  for (i = 0, avg = 0.0; i < A_size; ++i)
    {
      if (P[i] <= 0.0) avg += 0.0;

      else if (P[i] >= 1.0) avg += dn * log(dn / P[i]);

      else
	{
	  ln_p = log(P[i]);
	  ln_1_p = log(1.0 - P[i]);

	  for (j = 1, k = n - 1; j <= n; ++j, --k)
	    {
	      dj = (double)j;
	      dk = (double)k;
	      ln_combi = ln_fact[n] - ln_fact[j] - ln_fact[k];
	      prob = exp(ln_combi + dj * ln_p + dk * ln_1_p);

	      avg += prob * dj * log(dj / P[i]);
	    }
	}
    }

  if (n > 0) avg = avg / dn - log(dn);
  else avg = 0.0;

  return(avg);
}


/* Add (P[i] * Fudge) to each matrix-element. */
void fudge_align_mat(align_mat)
     double **align_mat;   /* The alignment matrix. */
{
  int i, j;


  /* Manipulate each position of each row. */
  for (i = 0; i < A_size; ++i)
    {
#ifdef NEW_FUDGE
      for (j = 0; j < L; ++j) align_mat[i][j] += P[i] * Fudge;
#else
      for (j = 0; j < L; ++j) align_mat[i][j] += Fudge;
#endif
    }
}


/* Determine the weight matrix by
 * 1) normalizing the fudged elements of the alignment matrix
 *    by the column sum and P[]; and
 * 2) taking the logarithm of the normalized elements. */
double **weight_matrix(align_mat)
     double **align_mat;       /* The alignment matrix. */
{
  int i, j;
  double col_sum;              /* The sum of the current column. */
  double **matrix;             /* The weight matrix. */


  /* Allocate space for "matrix[][]". */
  matrix = (double **)calloc_error(A_size, sizeof(double *),
				  "matrix", "weight_matrix");
  for (i = 0; i < A_size; ++i)
    matrix[i] = (double *)calloc_error(L, sizeof(double),
				       "matrix[]", "weight_matrix");

  /* Manipulate each position of each row. */
  for (j = 0; j < L; ++j)
    {
      /* Determine the column sum. */
      for (i = 0, col_sum = 0.0; i < A_size; ++i) col_sum += align_mat[i][j];
#ifndef NEW_FUDGE
      col_sum -= (double)(A_size - 1) * Fudge;
#endif

      /* Determine the weights. */
      for (i = 0; i < A_size; ++i)
	{
	  if (align_mat[i][j] == 0.0) matrix[i][j] = -INFINITY;
	  else matrix[i][j] = log(align_mat[i][j] / col_sum / P[i]);
	}
    }

  return(matrix);
}


/* Print the horizontal alignment matrix and its width. */
void print_h_align_mat(align_mat)
     double **align_mat;   /* The alignment matrix. */
{
  int i, j;


  /* Print the width of the matrix. */
  printf("width of the alignment matrix: %d\n", L);

  /* If alphabet does not have complements, or if some letter are their own
   * complements, print rows in the order they appear in the matrix. */
  if (Comp_flag != 1)
    {
      for (i = 0; i < A_size; ++i)
	{
	  /* Print the current row. */
	  (Ascii == YES) ? printf("%c |", A[i]) : printf("%3d |", A[i]);
	  for (j = 0; j < L; ++j) printf(" %3g", align_mat[i][j]);
	  printf("\n");
	}
    }

  /* If alphabet has complements, and no letter is its own complement,
   * rearrange rows so that the alphabet is symmetrical. */
  else
    {
      for (i = 0; i < A_size; i = i + 2)
	{
	  /* Print the current row. */
	  (Ascii == YES) ? printf("%c |", A[i]) : printf("%3d |", A[i]);
	  for (j = 0; j < L; ++j) printf(" %3g", align_mat[i][j]);
	  printf("\n");
	}
      for (i = A_size - 1; i > 0; i = i - 2)
	{
	  /* Print the current row. */
	  (Ascii == YES) ? printf("%c |", A[i]) : printf("%3d |", A[i]);
	  for (j = 0; j < L; ++j) printf(" %3g", align_mat[i][j]);
	  printf("\n");
	}
    }

  printf("\n");
}

/* Print the vertical alignment matrix and its width. */
void print_v_align_mat(align_mat)
     double **align_mat;   /* The alignment matrix. */
{
  int i, j;


  /* Print the width of the matrix. */
  printf("width of the alignment matrix: %d\n", L);

  /* If alphabet does not have complements, or if some letter are their own
   * complements, print rows in the order they appear in the matrix. */
  if (Comp_flag != 1)
    {
      /* Print the alphabet headings. */
      (Ascii == YES) ? printf("   %c ", A[0]) : printf("%3d", A[0]);
      for (i = 1; i < A_size; ++i)
	(Ascii == YES) ? printf("    %c ", A[i]) : printf(" %5d", A[i]);
      printf("\n");

      /* Print the columns. */
      for (j = 0; j < L; ++j)
	{
	  printf("%5g", align_mat[0][j]);
	  for (i = 1; i < A_size; ++i) printf(" %5g", align_mat[i][j]);
	  printf("\n");
	}
    }

  /* If alphabet has complements, and no letter is its own complement,
   * rearrange rows so that the alphabet is symmetrical. */
  else
    {
      /* Print the alphabet headings. */
      (Ascii == YES) ? printf("   %c ", A[0]) : printf("%3d", A[0]);
      for (i = 2; i < A_size; i = i + 2)
	(Ascii == YES) ? printf("    %c ", A[i]) : printf(" %5d", A[i]);
      for (i = A_size - 1; i > 0; i = i - 2)
	(Ascii == YES) ? printf("    %c ", A[i]) : printf(" %5d", A[i]);
      printf("\n");

      /* Print the columns. */
      for (j = 0; j < L; ++j)
	{
	  printf("%5g", align_mat[0][j]);
	  for (i = 2; i < A_size; i = i + 2)
	    printf(" %5g", align_mat[i][j]);
	  for (i = A_size - 1; i > 0; i = i - 2)
	    printf(" %5g", align_mat[i][j]);
	  printf("\n");
	}
    }

  printf("\n");
}
