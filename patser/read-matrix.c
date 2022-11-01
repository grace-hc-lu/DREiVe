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


/* Read a weight matrix.  The following external variables are used:
    Mat_file            M_fp                L_1
    Ascii       	Vertical 	    Letter_array
    Comp_flag   	Matrix   	    Letter_index
    A_size      	A                   L
*/


/* Reads the contents of the matrix file into a dynamic matrix, determine
 * "L" and "L_1", make sure the matrix is not empty, make sure the number
 * of elements is divisible by "A_size, and close matrix file. */
void read_matrix()
{
  int i;
  void read_h_sp_mat();     /* Read the horizontal weight matrix. */
  void read_v_sp_mat();     /* Read the vertical weight matrix. */
  void print_h_matrix();    /* Print the horizontal weight matrix. */
  void print_v_matrix();    /* Print the vertical weight matrix. */


  /* Allocate space for "Letter_array". */
  Letter_array = calloc_error(A_size, sizeof(char),
				"Letter_array", "read_matrix");
  for (i = 0; i < A_size; ++i) Letter_array[i] = NO;

  /* Allocate space for "Matrix[][]". */
  Matrix = (double **)calloc_error(A_size, sizeof(double *),
				    "Matrix", "read_matrix");
  for (i = 0; i < A_size; ++i) Matrix[i] = (double *)NULL;

  /* Read the weight matrix. */
  if (Vertical == NO)
    {
      read_h_sp_mat();                    /* Read horizontal weight matrix. */
      print_h_matrix(Matrix, L, Mat_file);/* Print horizontal weight matrix. */
    }
  else
    {
      read_v_sp_mat();                      /* Read vertical weight matrix. */
      print_v_matrix(Matrix, L, Mat_file);  /* Print vertical weight matrix. */
    }

  /* Determine "L_1". */
  L_1 = L - 1;

  /* Close the matrix file and free space for "Letter_array[]". */
  fclose(M_fp);
  free_error(Letter_array, "Letter_array", "read_matrix()");
}


/* Read the horizontal weight matrix.
 * Reads the contents of the matrix file into a dynamic matrix,
 * determine "L" and make sure the matrix is not empty. */
void read_h_sp_mat()
{
  int i;
  void read_first_sp_row(); /* Read the first row of the weight matrix. */
  void read_sp_row();       /* Read a row of the weight matrix,
			     * except for first row. */
  long find_line();         /* Find the beginning of a data line. */


  /* Read the first row of the matrix and determinine L. */
  read_first_sp_row();

  /* Determine the rest of the rows of the matrix. */
  for (i = 2; i <= A_size; ++i)
    read_sp_row(i);

  /* Make sure there are no illegal characters at the end of the file. */
  if (find_line(M_fp, Mat_file) != -1L)
    {
      fprintf(stderr, "Extra characters at the end of the matrix file.\n");
      exit(1);
    }
}


/* Read the vertical weight matrix.
 * Reads the contents of the matrix file into a dynamic matrix,
 * determine "L" and make sure the matrix is not empty. */
void read_v_sp_mat()
{
  int i;
  void read_letter_row();       /* Read a row of letters. */
  char read_double_row();       /* Read a row of integers. */


  /* Allocate space for "Letter_index". */
  Letter_index = (int *)calloc_error(A_size, sizeof(int),
			      "Letter_index", "read_v_sp_matrix");

  /* Read the letters in first row of matrix and intialize "Letter_index[]". */
  read_letter_row();

  /* Read the rows of doubles and increment "L". */
  L = 0;
  while (read_double_row() == YES) ++L;

  /* Make sure the matrix is not empty. */
  if (L == 0)
    {
      fprintf(stderr, "The first row of the weight matrix ");
      fprintf(stderr, "contains no numbers.\n\n");
      exit(1);
    }

  /* Truncate the rows of "Matrix[][]". */
  for (i = 0; i < A_size; ++i)
    Matrix[i] = (double *)realloc_error((char *)Matrix[i],
					(unsigned)L * sizeof(double),
					"Matrix[]", "read_v_sp_matrix()");

  /* Free space for "Letter_index". */
  free_error((char *)Letter_index, "Letter_index", "read_v_sp_mat()");
}

/* Read a row of doubles.  Returns YES if row is found; otherwise, NO. */
char read_double_row()
{
  double number;
  static int array_size = 0; /* The size of the arrays holding the rows. */
  int i;
  int read_double();         /* Read a double precision number from the matrix
			      * file.  The '|' character is ignored.  Returns
			      * YES: double found; NO: double not found. */
  long find_line();          /* Find the beginning of a data line. */
  

  if (find_line(M_fp, Mat_file) == -1L) return(NO);

  /* Make sure the summary matrix has enough space in each row. */
  if (L >= array_size)
    {
      array_size = array_size + CHUNK_SIZE;
      for (i = 0; i < A_size; ++i)
	Matrix[i] = (double *)recalloc_error((char *)Matrix[i],
		    array_size, sizeof(double), "Matrix[]", "read_double_row");
    }

  /* Read a row of the summary matrix. */
  for (i = 0; i < A_size; ++i)
    {
      if (read_double(M_fp, &number, L + 1) == YES)
	Matrix[Letter_index[i]][L] = number;
      else
	{
	  fprintf(stderr, "Row %d of the weight matrix ", L + 1);
	  fprintf(stderr, "contains too few numbers.\n");
	  exit(1);
	}
    }

  /* Make sure there are no extra numbers or other inappripriate
   * characters at the end of the line. */
  if (read_double(M_fp, &number, L + 1) == YES)
    {
      fprintf(stderr, "Row %d of the matrix contains more ", L + 1);
      fprintf(stderr, "elements than the number of columns.\n");
      exit(1);
    }

  return(YES);
}


/* Read the first row of the weight matrix and determinine L. */
void read_first_sp_row()
{
  int letter;                /* Letter index corresponding to current row. */
  double number;             /* The current matrix element. */
  int array_size= CHUNK_SIZE;/* The size of the array holding the first row. */
  int read_letter();         /* Read the letter index corresponding to row.
			      * Comment lines beginning with ';', '%', or '#'
			      * are ignored along with the '|'. */
  int read_double();         /* Read a double precision number from the matrix
			      * file.  The '|' character is ignored.  Returns
			      * YES: double found; NO: double not found. */


  /* Determine letter index of first row. */
  if ((letter = read_letter()) == EOF)
    {
      fprintf(stderr, "Cannot find the letter corresponding to the ");
      fprintf(stderr, "first row of the weight matrix.\n");
      exit(1);
    }

  /* Allocate initial space for the first row. */
  Matrix[letter] = (double *)calloc_error(array_size, sizeof(double),
					  "Matrix[]", "read_first_sp_row");

  /* Read the first row of the summary matrix. */
  while (read_double(M_fp, &number, 1) == YES)
    {
      Matrix[letter][L] = number;
      ++L;

      if (L >= array_size)
	{
	  array_size = array_size + CHUNK_SIZE;
	  Matrix[letter] = (double *)realloc_error((char *)Matrix[letter],
				     (unsigned)array_size * sizeof(double),
				     "Matrix[]", "read_first_sp_row");
	}
    }

  /* Truncate the row, if the matrix is not empty. */
  if (L == 0)
    {
      fprintf(stderr, "The first row of the weight matrix ");
      fprintf(stderr, "contains no numbers.\n\n");
      exit(1);
    }
  else
    Matrix[letter] = (double *)realloc_error((char *)Matrix[letter],
				(unsigned)L * sizeof(double),
				"Matrix[]", "read_first_sp_row");
}

/* Read a row of the weight matrix, except for first row. */
void read_sp_row(row)
     int row;           /* The row currently being read. */
{
  int letter;                /* Letter index corresponding to current row. */
  double number;             /* The current matrix element. */
  int i;
  int read_letter();         /* Read the letter index corresponding to row.
			      * Comment lines beginning with ';', '%', or '#'
			      * are ignored along with the '|'. */
  int read_double();         /* Read a double precision number from the matrix
			      * file.  The '|' character is ignored.  Returns
			      * YES: double found; NO: double not found. */


  /* Determine letter index of the current row. */
  if ((letter = read_letter()) == EOF)
    {
      fprintf(stderr, "Cannot find the letter corresponding to ");
      fprintf(stderr, "row %d of the weight matrix.\n", row);
      exit(1);
    }

  /* Allocate space for the current row. */
  Matrix[letter] = (double *)calloc_error(L, sizeof(double),
					  "Matrix[]", "read_sp_row");

  for (i = 0; i < L; ++i)
    {
      if (read_double(M_fp, &number, row) == YES)
	Matrix[letter][i] = number;
      else
	{
	  fprintf(stderr, "Row %d of the weight matrix has ", row);
	  fprintf(stderr, "less elements than row 1.\n");
	  exit(1);
	}
    }

  /* Make sure there are no extra numbers at the end of the line. */
  if (read_double(M_fp, &number, row) == YES)
    {
      fprintf(stderr, "Row %d of the weight matrix has ", row);
      fprintf(stderr, "more elements than row 1.\n");
      exit(1);
    }
}


/* Read a double precision number from the matrix file.  The '|' character
 * is ignored.  Returns YES: double found; NO: double not found. */
int read_double(fp, number, row)
     FILE *fp;
     double *number;    /* The address of the number to be determined. */
     int row;           /* The number of the current matrix row. */
{
  int letter;     /* A letter in the matrix file. */

  while (1)
    {
      letter = getc(fp);

      /* Check for the end of the line and
       * pass over any other whitespace or '|'. */
      if (letter == '\n')
	return(NO);
      else if (isspace(letter) != NO);
      else if (letter == '|');

      /* Treat a comment as if it is the end of the line. */
      else if ((letter == ';') || (letter == '#') || (letter == '%'))
	{
	  while (((letter = getc(fp)) != '\n') && (letter != EOF));
	  return(NO);
	}

      /* Read the number. */
      else if ((ungetc(letter, fp) == EOF)||(fscanf(M_fp, "%lf", number) != 1))
	{
	  fprintf(stderr, "Problems reading a number from row %d ", row);
	  fprintf(stderr, "of the matrix file.\n");
	  exit(1);
	}
      else
	return(YES);
    }
}


/* Print the horizontal matrix and its width. */
void print_h_matrix(matrix, width, mat_file)
     double **matrix;        /* The matrix in the matrix file. */
     int width;              /* Width of the matrix. */
     char *mat_file;         /* Name of the matrix file. */
{
  int i, j;
  void print_double();       /* Print double in the %f or %e formats to
			      * fit the designated space. */


  /* Print the width of the matrix. */
  printf("width of the matrix in file \"%s\": %d\n", mat_file, width);

  /* If alphabet does not have complements, or if some letter are their own
   * complements, print rows in the order they appear in the matrix. */
  if (Comp_flag != 1)
    {
      for (i = 0; i < A_size; ++i)
	{
	  /* Print the current row. */
	  (Ascii == YES) ? printf("%c |", A[i]) : printf("%3d |", A[i]);
	  for (j = 0; j < width; ++j)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 7, 2);
	    }
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
	  for (j = 0; j < width; ++j)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 7, 2);
	    }
	  printf("\n");
	}
      for (i = A_size - 1; i > 0; i = i - 2)
	{
	  /* Print the current row. */
	  (Ascii == YES) ? printf("%c |", A[i]) : printf("%3d |", A[i]);
	  for (j = 0; j < width; ++j)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 7, 2);
	    }
	  printf("\n");
	}
    }

  printf("\n");
}



/* Print the vertical matrix and its width. */
void print_v_matrix(matrix, width, mat_file)
     double **matrix;        /* The matrix in the matrix file. */
     int width;              /* Width of the matrix. */
     char *mat_file;         /* Name of the matrix file. */
{
  int i, j;
  void print_double();       /* Print double in the %f or %e formats to
			      * fit the designated space. */


  /* Print the width of the matrix. */
  printf("width of the matrix in file \"%s\": %d\n", mat_file, width);



  /* If alphabet does not have complements, or if some letter are their own
   * complements, print rows in the order they appear in the matrix. */
  if (Comp_flag != 1)
    {
      /* Print the alphabet headings. */
      (Ascii == YES) ? printf("     %c  ", A[0]) : printf("%8d", A[0]);
      for (i = 1; i < A_size; ++i)
	(Ascii == YES) ? printf("      %c  ", A[i]) : printf(" %8d", A[i]);
      printf("\n");

      /* Print the columns. */
      for (j = 0; j < width; ++j)
	{
	  print_double(stdout, matrix[0][j], 8, 2);
	  for (i = 1; i < A_size; ++i)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 8, 2);
	    }
	  printf("\n");
	}
    }

  /* If alphabet has complements, and no letter is its own complement,
   * rearrange rows so that the alphabet is symmetrical. */
  else
    {
      /* Print the alphabet headings. */
      (Ascii == YES) ? printf("     %c  ", A[0]) : printf("%8d", A[0]);
      for (i = 2; i < A_size; i = i + 2)
	(Ascii == YES) ? printf("      %c  ", A[i]) : printf(" %8d", A[i]);
      for (i = A_size - 1; i > 0; i = i - 2)
	(Ascii == YES) ? printf("      %c  ", A[i]) : printf(" %8d", A[i]);
      printf("\n");

      /* Print the columns. */
      for (j = 0; j < width; ++j)
	{
	  print_double(stdout, matrix[0][j], 8, 2);
	  for (i = 2; i < A_size; i = i + 2)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 8, 2);
	    }
	  for (i = A_size - 1; i > 0; i = i - 2)
	    {
	      printf(" ");
	      print_double(stdout, matrix[i][j], 8, 2);
	    }
	  printf("\n");
	}
    }

  printf("\n");
}

/* Print a double using either the %f or %e formats.
 * If the number can fit in the specified width use the %f format.
 * If the |number| is less than 0.1 or too large to fit, use the %e format. */
void print_double(fp, x, width, precision)
     FILE *fp;
     double x;       /* The number being printed. */
     int width;      /* The space for the printed number (>= 7). */
     int precision;  /* The precision for the printed number
		      * in the %f format (>= 1 and <= width - 3). */
{
  double logten;             /* log_10 of numbers absolute value. */
  int adj;                   /* Precision adjustment when number < 0. */
  int adj_prec;              /* Precision adjusted for the %e format. */
  double log_low_cutoff;     /* Cutoff logarithm for the %f format. */


  if (precision > width - 3)
    {
      fprintf(stderr, "\"print_double()\" has been passed a ");
      fprintf(stderr, "precision greater than width minus 3.\n");
      exit(1);
    }
  else if (precision < 0)
    {
      fprintf(stderr, "\"print_double()\" has been passed a ");
      fprintf(stderr, "precision less than 0.\n");
      exit(1);
    }

  if (width < 7)
    {
      fprintf(stderr, "\"print_double()\" has been passed a width ");
      fprintf(stderr, "less than 7.\n");
      exit(1);
    }

  if (x == 0.0)
    {
      fprintf(fp, "%*.*f", width, precision, x);
      return;
    }
  else if (x > 0.0) adj = 0;
  else adj = 1;

  if (precision < 1) log_low_cutoff = 0.0;
  else log_low_cutoff = -1.0;

  logten = log10(fabs(x));

  if (logten >= 100.0)
    {
      adj_prec = width - 7 - adj;
      if (adj_prec < 0) adj_prec = 0;
      fprintf(fp, "%*.*e", width, adj_prec, x);
    }
  else if (logten >= (double)(width - precision - 1 - adj))
    fprintf(fp, "%*.*e", width, width - 6 - adj, x);
  else if (logten >= log_low_cutoff) fprintf(fp, "%*.*f", width, precision, x);
  else if (logten >= -99.0) fprintf(fp, "%*.*e", width, width - 6 - adj, x);
  else
    {
      adj_prec = width - 7 - adj;
      if (adj_prec < 0) adj_prec = 0;
      fprintf(fp, "%*.*e", width, adj_prec, x);
    }
}
