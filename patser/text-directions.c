/* Copyright 1990, 1994, 1995, 1996, 2000, 2001 Gerald Z. Hertz
 * May be copied for noncommercial purposes.
 *
 * Author:
 *   Gerald Z. Hertz
 *   hertz@colorado.edu
 */


#ifndef OPTIONS
#include <stdio.h>
#else
#include "options.h"
#endif

/* The text of the directions. */
void text_directions(FILE *fp)
{
  fprintf(fp, "Copyright 1990, 1994, 1995, 1996, 2000, 2001 Gerald Z. ");
  fprintf(fp, "Hertz\nMay be copied for noncommercial purposes.\n\n");

  fprintf(fp, "Author:\n  Gerald Z. Hertz\n  hertz@colorado.edu\n\n");

  fprintf(fp, "\n");

  fprintf(fp, "\n");

  fprintf(fp, "PATSER (version 3d)\n\n");

  fprintf(fp, "This program scores the L-mers (subsequences of length ");
  fprintf(fp, "L) of the\nindicated sequences against the indicated al");
  fprintf(fp, "ignment or weight matrix.\nThe elements of an alignment");
  fprintf(fp, " matrix are simply the number of times\nthat the indica");
  fprintf(fp, "ted letter is observed at the indicated position of a\n");
  fprintf(fp, "sequence alignment.  Such elements must be processed be");
  fprintf(fp, "fore the matrix\ncan be used to score an L-mer (e.g., H");
  fprintf(fp, "ertz and Stormo, 1999,\nBioinformatics, 15:563-577).  A");
  fprintf(fp, " weight matrix is a matrix whose\nelements are in a for");
  fprintf(fp, "m considered appropriate for scoring an L-mer.\n\n");

  fprintf(fp, "Each element of an alignment matrix is converted to an ");
  fprintf(fp, "element of a\nweight matrix by first adding pseudo-coun");
  fprintf(fp, "ts in proportion to the a\npriori probability of the co");
  fprintf(fp, "rresponding letter (see option \"-b\" in\nsection 1 bel");
  fprintf(fp, "ow) and dividing by the total number of sequences plus\n");
  fprintf(fp, "the total number of pseudo-counts.  The resulting frequ");
  fprintf(fp, "ency is\nnormalized by the a priori probability for the");
  fprintf(fp, " corresponding letter.\nThe final quotient is converted");
  fprintf(fp, " to an element of a weight matrix by\ntaking its natura");
  fprintf(fp, "l logarithm.  The use of pseudo-counts here differs\nfr");
  fprintf(fp, "om previous versions of this program by being proportio");
  fprintf(fp, "nal to the a\npriori probability.\n\n");

  fprintf(fp, "Version 3 of this program differs from previous version");
  fprintf(fp, "s by also\nnumerically estimating the p-value of the sc");
  fprintf(fp, "ores.  The p-value\ncalculated here is the probability ");
  fprintf(fp, "of observing a particular score or\nhigher at a particu");
  fprintf(fp, "lar sequence position and does NOT account for the\ntot");
  fprintf(fp, "al amount of sequence being scored.  P-values are estim");
  fprintf(fp, "ated by the\nmethod described in Staden, 1989, CABIOS, ");
  fprintf(fp, "p. 89--96.  The relative\nvalue for each element of the");
  fprintf(fp, " weight matrix is approximated by\nintegers in a range ");
  fprintf(fp, "determined by the \"-R\" and \"-M\" options (section 6\n");
  fprintf(fp, "below).  The p-value is calculated for each possible in");
  fprintf(fp, "teger score and\nthe values are stored.  The actual sco");
  fprintf(fp, "res for the sequences are\ndetermined from the true wei");
  fprintf(fp, "ght matrix.  The true scores are converted\nto their co");
  fprintf(fp, "rresponding integer values and their p-values are looke");
  fprintf(fp, "d up.\n\n");

  fprintf(fp, "Matrices can be either horizontal or vertical.  In a ho");
  fprintf(fp, "rizontal\nmatrix, the columns correspond to the positio");
  fprintf(fp, "ns within the pattern,\nand the rows correspond to the ");
  fprintf(fp, "letters.  Each row begins with the\ncorresponding lette");
  fprintf(fp, "r (or integer, if the \"-i\" option is used).  In a\nve");
  fprintf(fp, "rtical matrix, the rows correspond to the positions wit");
  fprintf(fp, "hin the\npattern, and the columns correspond to the let");
  fprintf(fp, "ters.  The first row\ncontains the letters (or integers");
  fprintf(fp, ", if the \"-i\" option is used)\ncorresponding to each ");
  fprintf(fp, "column.  In both types of matrices, spaces,\ntabs, and ");
  fprintf(fp, "vertical bars (|) are ignored.  The output of the \"con");
  fprintf(fp, "sensus\"\nand \"wconsensus\" programs consists of horiz");
  fprintf(fp, "ontal alignment matrices.\n\n");

  fprintf(fp, "The input files can contain comments according to the f");
  fprintf(fp, "ollowing\nconvention.  The portion of a line following ");
  fprintf(fp, "a ';', '%%', or '#' is\nconsidered a comment and is ign");
  fprintf(fp, "ored.  Comments can begin anywhere in a\nline and alway");
  fprintf(fp, "s end at the end of the line.  The output of this\nprog");
  fprintf(fp, "ram is sent to the standard output.\n\n");

  fprintf(fp, "The following options can be determined on the command ");
  fprintf(fp, "line.\n\n");

  fprintf(fp, "  0) -h: print these directions.\n\n");

  fprintf(fp, "  1) Matrix options.\n     -m filename: (default name i");
  fprintf(fp, "s \"matrix\") file containing the matrix.\n     -w: the");
  fprintf(fp, " matrix is a weight matrix (default: alignment matrix)\n");
  fprintf(fp, "     -b number: a non-negative number indicating the to");
  fprintf(fp, "tal number of\n         pseudo-counts added to each ali");
  fprintf(fp, "gnment position (default: 1).\n         Before converti");
  fprintf(fp, "ng an alignment matrix to a weight matrix, the\n       ");
  fprintf(fp, "  total pseudo-counts multiplied by the a priori probab");
  fprintf(fp, "ility\n         (see section 3 below) of the correspond");
  fprintf(fp, "ing letter is added\n         to each matrix element.\n");
  fprintf(fp, "     -v: the matrix is a vertical matrix (default: hori");
  fprintf(fp, "zontal matrix).\n     -p: print the weight matrix deriv");
  fprintf(fp, "ed from the alignment matrix.\n\n");

  fprintf(fp, "  2) -f filename: this file (default: read from the sta");
  fprintf(fp, "ndard input) contains\n        the names of the sequenc");
  fprintf(fp, "es.  The corresponding sequence may follow\n        its");
  fprintf(fp, " name if the sequence is enclosed between backslashes (");
  fprintf(fp, "\\).\n        Otherwise, the sequence is assumed to be ");
  fprintf(fp, "in a separate file having\n        the indicated name.\n");
  fprintf(fp, "\n");

  fprintf(fp, "        In the sequences, whitespace, slashes (/), peri");
  fprintf(fp, "ods, dashes (unless\n        part of an integer when th");
  fprintf(fp, "e \"-i\" option is used), and comments\n        beginni");
  fprintf(fp, "ng with ';', '%%', or '#' are ignored.  When using lett");
  fprintf(fp, "er\n        characters (i.e., with the \"-a\" or \"-A\"");
  fprintf(fp, " alphabet option), integers\n        are also ignored s");
  fprintf(fp, "o that the sequence file can contain positional\n      ");
  fprintf(fp, "  information.  When using integer characters (i.e., wi");
  fprintf(fp, "th the \"-i\"\n        alphabet option) the integers mu");
  fprintf(fp, "st be separated by whitespace.\n\n");

  fprintf(fp, "        A \"-c\" preceding the name of a sequence file ");
  fprintf(fp, "indicates that the\n        corresponding sequence is c");
  fprintf(fp, "ircular.\n\n");

  fprintf(fp, "  3) Alphabet options---the three options in this secti");
  fprintf(fp, "on are mutually\n     exclusive (default: \"-a alphabet");
  fprintf(fp, "\").  The a priori probabilities mentioned\n     below ");
  fprintf(fp, "are used when converting an alignment matrix to a weigh");
  fprintf(fp, "t matrix.\n     -a filename: file containing the alphab");
  fprintf(fp, "et and normalization information.\n\n");

  fprintf(fp, "        Each line contains a letter (a symbol in the al");
  fprintf(fp, "phabet) followed by an\n        optional normalization ");
  fprintf(fp, "number (default: 1.0).  The normalization is\n        b");
  fprintf(fp, "ased on the relative a priori probabilities of the lett");
  fprintf(fp, "ers.  For\n        nucleic acids, this might be be the ");
  fprintf(fp, "genomic frequency of the bases\n        or the frequenc");
  fprintf(fp, "ies observed in the data used to generate the alignment");
  fprintf(fp, ".\n        In nucleic acid alphabets, a letter and its ");
  fprintf(fp, "complement appear on the\n        same line, separated ");
  fprintf(fp, "by a colon (a letter can be its own complement,\n      ");
  fprintf(fp, "  e.g. when using a dimer alphabet). Complementary lett");
  fprintf(fp, "ers may use the\n        same normalization number.  On");
  fprintf(fp, "ly the standard 26 letters are\n        permissible; ho");
  fprintf(fp, "wever, when the \"-CS\" option is used, the alphabet is");
  fprintf(fp, "\n        case sensitive so that a total of 52 differen");
  fprintf(fp, "t characters are possible.\n\n");

  fprintf(fp, "        POSSIBLE LINE FORMATS WITHOUT COMPLEMENTARY LET");
  fprintf(fp, "TERS:\n        letter\n        letter normalization\n\n");

  fprintf(fp, "        POSSIBLE LINE FORMATS WITH COMPLEMENTARY LETTER");
  fprintf(fp, "S:\n        letter:complement\n        letter:complemen");
  fprintf(fp, "t normalization\n        letter:complement normalizatio");
  fprintf(fp, "n:complement's_normalization\n\n");

  fprintf(fp, "     -i filename: same as the \"-a\" option, except tha");
  fprintf(fp, "t the symbols of\n        the alphabet are represented ");
  fprintf(fp, "by integers rather than by letters.\n        Any intege");
  fprintf(fp, "r permitted by the machine is a permissible symbol.\n\n");

  fprintf(fp, "     -A alphabet_and_normalization_information: same as");
  fprintf(fp, " \"-a\" option, except\n        information appears on ");
  fprintf(fp, "the command line (e.g., -A a:t 3 c:g 2).\n\n");

  fprintf(fp, "  4) Alphabet modifiers indicating whether ascii alphab");
  fprintf(fp, "ets are case\n     sensitive---the two options in this ");
  fprintf(fp, "section are mutually exclusive\n     with each other an");
  fprintf(fp, "d with the \"-i\" option (default: ascii alphabets are\n");
  fprintf(fp, "     case insensitive).\n     -CS: ascii alphabets are ");
  fprintf(fp, "case sensitive.\n     -CM: ascii alphabets are case ins");
  fprintf(fp, "ensitive, but mark the location\n          of lowercase");
  fprintf(fp, " letters by printing a line containing their locations.");
  fprintf(fp, "\n          This option is useful when lowercase letter");
  fprintf(fp, "s indicate a functional\n          landmark such as a t");
  fprintf(fp, "ranscriptional start in a DNA sequence.\n\n");

  fprintf(fp, "  5) Options for adjusting or restricting which scores ");
  fprintf(fp, "are printed.\n     The \"-ls\", \"-li\", and \"-lp\" op");
  fprintf(fp, "tions are mutually exclusive.\n     -c: also score the ");
  fprintf(fp, "complementary sequences.  The complements are\n        ");
  fprintf(fp, "determined by the program and are not explicitly stated");
  fprintf(fp, " in the\n        sequence files.\n     -ls number: lowe");
  fprintf(fp, "r threshold for printing scores, inclusive\n           ");
  fprintf(fp, "      (formerly the -l option).\n     -li: assume that ");
  fprintf(fp, "the maximum ln(p-value) for printing scores equals\n   ");
  fprintf(fp, "       the negative of the sample-size adjusted informa");
  fprintf(fp, "tion content;\n          indirectly determines the lowe");
  fprintf(fp, "r threshold for printing scores.\n     -lp number: the ");
  fprintf(fp, "maximum ln(p-value) for printing scores; indirectly\n  ");
  fprintf(fp, "               determines the lower threshold for print");
  fprintf(fp, "ing scores.\n     -u number: upper threshold for printi");
  fprintf(fp, "ng scores, exclusive.\n\n");

  fprintf(fp, "     -t: just print the top score for each sequence.\n ");
  fprintf(fp, "    -t number: print the indicated number of top scores");
  fprintf(fp, " for each sequence.\n     -ds: if the \"-t number\" opt");
  fprintf(fp, "ion is used, print the top scores for each\n          s");
  fprintf(fp, "equence in the order of decreasing score (default: orde");
  fprintf(fp, "r the\n          scores according to their position wit");
  fprintf(fp, "hin the sequence).\n     -e number: the small differenc");
  fprintf(fp, "e for considering 2 scores equal\n                (defa");
  fprintf(fp, "ult: 0.000001)\n\n");

  fprintf(fp, "  6) Options indicating how unrecognized symbols are tr");
  fprintf(fp, "eated (default: -d1).\n     Symbols are letters when op");
  fprintf(fp, "tion \"-a\" or \"-A\" is used;\n     symbols are intege");
  fprintf(fp, "rs when option \"-i\" is used.\n     The following thre");
  fprintf(fp, "e options are mutually exclusive.\n     -d0: treat unre");
  fprintf(fp, "cognized symbols as errors and exit the program.\n     ");
  fprintf(fp, "-d1: treat unrecognized symbols as discontinuities, but");
  fprintf(fp, " print a warning.\n          Treating a symbol as a dis");
  fprintf(fp, "continuity means that any L-mer\n          containing t");
  fprintf(fp, "he unrecognized symbol will be ignored.\n     -d2: trea");
  fprintf(fp, "t unrecognized symbols as discontinuities, and print NO");
  fprintf(fp, " warning.\n\n");

  fprintf(fp, "  7) Options for adjusting the estimation of p-value.\n");
  fprintf(fp, "     If the -R option is set to zero, the p-value is no");
  fprintf(fp, "t estimated.\n     -R number: the range for approximati");
  fprintf(fp, "ng a column of the weight matrix with\n                ");
  fprintf(fp, "integers (default: 10000).  This number is the differen");
  fprintf(fp, "ce\n                between the largest and smallest in");
  fprintf(fp, "tegers used to estimate\n                the scores.  H");
  fprintf(fp, "igher values increase precision, but will take\n       ");
  fprintf(fp, "         longer to calculate the table of p-values.\n  ");
  fprintf(fp, "   -M number: the minimum score for approximating p-val");
  fprintf(fp, "ues (default: 0).\n                Higher values will i");
  fprintf(fp, "ncrease precision,\n                but may miss intere");
  fprintf(fp, "sting scores.\n");
}
