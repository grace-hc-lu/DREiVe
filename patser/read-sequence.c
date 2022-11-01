/* Copyright 1991, 1994, 1996 Gerald Z. Hertz
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

/* If defined, print the sequence: USED FOR DEBUGGING PURPOSES ONLY.
#define PRINT_SEQUENCE 1
 */

/* Functions for reading a sequence from a sequence file.
 * The following external variables are used:
 *
 *    Error                 A                     Circle
 *    Ascii     	    S_fp		  Sequence
 *    A_size    	    File		  Seq_length
 *              	    Fp                    Case_sensitive
 */


/* Function for reading a letter from the sequence.  Equals either
 * "read_ascii_letter()" or "read_int_letter()".  Returns EOF: no more
 * letters; NO: letter is unrecognized; YES: letter is recognized. */
static int (*read_seq_letter)();
static int STseq_array = 0;   /* Amount of space allocated for "Sequence[]". */

/* Initialize information for "read_sequence()". */
void init_read_seq()
{
  int read_ascii_letter(); /* Read an ASCII letter from a sequence. */
  int read_int_letter();   /* Read an integer letter from a sequence. */


   /* Initialize "(*read_seq_letter)()". */
  if (Ascii == YES) read_seq_letter = read_ascii_letter;
  else read_seq_letter = read_int_letter;
}



/* Open sequence file, determine whether it is circular, and read the sequence.
 * Returns YES, if sequence exists; returns NO, if sequence does not exist.*/
int read_sequence()
{
  int letter;            /* Index of the letter just read. */
  int status;            /* Return status of (*read_seq_letter)(). */
  int open_sequence();   /* Open sequence file, determine whether circular. */
#ifdef PRINT_SEQUENCE
  int i;
#endif


  /* Open a sequence file and determine whether the sequence is circular.
   * Returns NO, if there are not any more sequences listed. */
  if (open_sequence() == NO)
    {
      if (S_fp != stdin) fclose(S_fp);
      return(NO);
    }

  /* Initialize the line for marking lowercase letters. */
  if ((Ascii == YES) && (Case_sensitive == 2))
    printf("%20s lower case letters:", File);

  /* Read the letters in the sequence file. */
  for (Seq_length=0; (status=(*read_seq_letter)(&letter)) != EOF; ++Seq_length)
    {
      /* Make sure that "Sequence[]" is big enough to hold the sequence. */
      if (Seq_length >= STseq_array)
	{
	  STseq_array += CHUNK_SIZE;
	  Sequence = (LETTER *)recalloc_error((char *)Sequence, STseq_array,
				sizeof(LETTER), "Sequence", "read_sequence()");
	}

      /* Insert the letter into the "Sequence[]". */
      if (status == YES) Sequence[Seq_length] = (LETTER)letter;
      else if ((Error == 1) || (Error == 2))
	Sequence[Seq_length] = (LETTER)(A_size + 1);
      else exit(1);
    }
  if (Fp != S_fp) fclose(Fp);

  /* End the line for marking lowercase letters. */
  if ((Ascii == YES) && (Case_sensitive == 2)) printf("\n");

  /* Truncate the sequence to the appropriate length. */
  STseq_array = Seq_length;
  Sequence = (LETTER *)recalloc_error((char *)Sequence, STseq_array,
				sizeof(LETTER), "Sequence", "read_sequence()");

#ifdef PRINT_SEQUENCE
  /* Print the current sequence: FOR DEBUGGING PURPOSES ONLY. */
  printf("%s: ", File);
  for (i = 0; i < Seq_length; ++i)
    {
      if (Ascii == YES) printf("%c", A[(int)Sequence[i]]);
      else printf(" %d", A[(int)Sequence[i]]);
    }
  printf("\n");
#endif

  return(YES);
}

/* Read the name of a sequence and determine whether it is circular.
 * Determine whether or not the sequence is in a separate file.
 * Determine the pointer to the file containing the sequence.
 * Returns YES, if sequence exists; returns NO, if sequence does not exist. */
int open_sequence()
{
  int find_line_2();       /* Find the beginning of a data line. */


  /* Find the beginning of a data line. */
  if (find_line_2(S_fp) == EOF) return(NO);

  /* Initialize the flag that indicates whether sequence is circular. */
  Circle = NO;
  
  /* Read the initial string in the "line" of the file of sequence files. */
  if (fscanf(S_fp, "%s", File) == EOF) return(NO);
  
  /* Determine whether the sequence is circular. */
  if (!strcmp(File, "-c"))
    {
      Circle = YES;
      if (fscanf(S_fp, "%s", File) == EOF)
	{
	  fprintf(stderr, "The list of files has the wrong format.\n");
	  fprintf(stderr, "It contains a \"-c\" that is not followed ");
	  fprintf(stderr, "by another character.\n");
	  exit(1);
	}
    }
  

  /* Determine whether the sequence is in a separate file,
   * or with the list of sequence names. */
  if (find_line_2(S_fp) == '\\')
    {
      if (getc(S_fp) != (int)'\\')
	{
	  fprintf(stderr, "A \"getc\" error while reading the ");
	  fprintf(stderr, "sequence information.\n");
	  exit(1);
	}
      Fp = S_fp;
      return(YES);
    }

  /* Open sequence file, if sequence is contained in separate file. */
  else
    {
      if ((Fp = fopen(File, "r")) == (FILE *)NULL)
	{
	  fprintf(stderr, "Sequence file \"%s\" cannot be opened.\n", File);
	  exit(1);
	}
      return(YES);
    }
}


/* Function for reading an ASCII letter from the sequence.  Returns EOF: no
 * more letters; NO: letter is unrecognized; YES: letter is recognized. */
int read_ascii_letter(letter_index)
     int *letter_index;       /* Index of the current letter. */
{
  int letter;                 /* The ASCII code for the current letter. */
  int i;


  while (((letter = getc(Fp)) != EOF) && (letter != '\\'))
    {
      /* Ignore the letter, if it is whitespace, a digit, '-', '.', or '/'. */
      if (isspace(letter) != NO);
      else if (isdigit(letter) != NO);
      else if (letter == '-');
      else if (letter == '.');
      else if (letter == '/');

      /* Scan over any comments. */
      else if ((letter == ';') || (letter == '#') || (letter == '%'))
	while (((letter = getc(Fp)) != '\n') && (letter != EOF));

      /* Determine the index of the letter. */
      else if (isalpha(letter) != NO)
	{
	  if ((Case_sensitive != 0) && (islower(letter) != NO))
	    {
	      letter = toupper(letter);
	      if (Case_sensitive == 2) printf(" %d", Seq_length + 1);
	    }
	  for (i = 0; i < A_size; ++i)
	    {
	      if (letter == A[i])
		{
		  *letter_index = i;
		  return(YES);
		}
	    }
	  /* Error message if the letter is not recognized. */
	  if ((Error == 0) || (Error == 1))
	    {
	      fprintf(stderr, "The following letter in sequence ");
	      fprintf(stderr, "\"%s\" is not recognized: %c\n", File, letter);
	    }
	  return(NO);
	}

      else
	{
	  fprintf(stderr, "Ascii code \"%d\" is an invalid ", letter);
	  fprintf(stderr, "character in sequence \"%s\".\n", File);
	  exit(1);
	}
    }

  /* Return EOF if no more letters. */
  if (letter == '\\') return(EOF);
  else if (feof(Fp) == NO)
    {
      fprintf(stderr, "An error occurred while reading sequence ");
      fprintf(stderr, "file \"%s\".\n", File );
      exit(1);
    }
  else return(EOF);
}

/* Function for reading an integer letter from the sequence.  Returns EOF: no
 * more letters; NO: letter is unrecognized; YES: letter is recognized. */
int read_int_letter(letter_index)
     int *letter_index;       /* Index of the current letter. */
{
  int letter;                 /* The integer code for the current letter. */
  int status;                 /* Status of whether an integer was read. */
  int i;


  while ((status = fscanf(Fp, "%d", &letter)) != EOF)
    {
      /* Determine the index of the letter. */
      if (status == 1)
	{
	  for (i = 0; i < A_size; ++i)
	    {
	      if (letter == A[i])
		{
		  *letter_index = i;
		  return(YES);
		}
	    }
	  /* Error message if the integer letter is not recognized. */
	  if ((Error == 0) || (Error == 1))
	    {
	      fprintf(stderr, "The following integer in sequence ");
	      fprintf(stderr, "\"%s\" is not recognized: %d\n", File, letter);
	    }
	  return(NO);
	}

      /* What to do if the next characters do not form an integer. */
      else if (status == 0)
	{
	  letter = getc(Fp);
	  if (letter == '\\') break;

	  /* Ignore a letter, if it is whitespace or a '/'. */
	  if (isspace(letter) != NO);
	  else if (letter == '/');

	  /* Scan over any comments. */
	  else if ((letter == ';') || (letter == '#') || (letter == '%'))
	    while (((letter = getc(Fp)) != '\n') && (letter != EOF));

	  else
	    {
	      fprintf(stderr, "The letter \"%c\" (ASCII character ", letter);
	      fprintf(stderr, "%d) in sequence\n\"%s\" ", letter, File);
	      fprintf(stderr, "is not a recognized character.\n");
	      exit(1);
	    }
	}

      else bug_report("read_int_letter()");
    }

  /* Return EOF if no more letters. */
  if (letter == '\\') return(EOF);
  else if (feof(Fp) == NO)
    {
      fprintf(stderr, "An error occurred while reading sequence ");
      fprintf(stderr, "file \"%s\".\n", File );
      exit(1);
    }
  else return(EOF);
}


/* Find the beginning of a data line.  Returns the ASCII code for the first
 * character of the data line; returns EOF at the end of the file. */
int find_line_2(fp)
     FILE *fp;
{
  int letter;            /* The letter just read from the file. */

  while ((letter = getc(fp)) != EOF)
    {
      if (isspace(letter) != NO);
      else if ((letter == ';') || (letter == '%') || (letter == '#'))
	while ((letter != '\n') && (letter != EOF)) letter = getc(fp);
      else if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr, "An \"ungetc\" error while reading the ");
	  fprintf(stderr, "matrix information.\n");
	  exit(1);
	}
      else return(letter);
    }

  /* Double check of EOF. */
  if (feof(fp) == NO)
    {
      fprintf(stderr, "A  \"getc\" error occurred while reading input.\n");
      exit(1);
    }
  else return(EOF);
}
