/* Global alignment by dynamic programming:
 * output the optimal alignment score, and one point in the optimal alignment:
 * the best alignment for the midpoint (residue M/2) of sequence 1.
 *
 * To build the "global" program:
 *     gcc -g -Wall -o global global.c fasta.c
 * To run:
 *     ./global HBA_HUMAN HBB_HUMAN
 *     
 * SRE, Mon Sep  9 13:06:02 2002    
 */    

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fasta.h"

int
main(int argc, char **argv)
{
  char      *file1, *file2;     /* name of FASTA files 1 and 2 from command line */
  char      *s;                 /* tmp storage for a seq before we shift it +1 into 1..L coord system */
  char      *seq1, *seq2;       /* sequences: 1..M and 1..N coord system */
  char      *name1, *name2;     /* names of the two sequences (unused) */
  int        M,N;		/* lengths of seq 1 and seq2 */
  FASTAFILE *ffp;               /* ptr to an open FASTA file that we're reading */
  int      **mx;                /* dynamic programming matrix */
  int        i,j;               /* residue coords in seq1 and seq2 */
  int        mid;		/* midpoint in seq1, that we're finding optimal alignment for */
  int        usedj;             /* optimal alignment for residue seq1[mid] */
  int        gapcost;           /* linear gap cost to apply to an insertion */
  int        mismat;		/* substitution penalty */
  int        mat;		/* match score */
  int        sc;		/* tmp variable used to find best path */

  /* Set up our simple scoring system.
   */
  gapcost = -6;
  mat     = 5;
  mismat  = -2;

  /* Parse the command line
   */
  if (argc != 3) 
    {
      fprintf(stderr, "Wrong number of command line arguments.\nUsage: nw <file1> <file2>\n");
      exit(1);
    }
  file1 = argv[1];
  file2 = argv[2];

  /* Get the first sequence out of both files.
   * There's a little tidbit in here that shifts the sequence from 0..L-1 coordinates
   * to 1..L, which will make it convenient to set up our initialization conditions in
   * row 0 and column 0. Some texts would simply do "s--;", decrementing the pointer by 1,
   * which will usually work as long as you never dereference s[0] - but it's a trick
   * that's technically illegal C.                                                               
   */
  if ((ffp = OpenFASTA(file1)) == NULL) { fprintf(stderr, "Failed to open file %s\n", file1); exit(1); }
  if (! ReadFASTA(ffp, &s, &name1, &M)) { fprintf(stderr, "No data in file %s\n", file1);     exit(1); }
  seq1 = malloc(sizeof(char) * (M+2));
  strcpy(seq1+1, s);		/* we can now index the residues in seq1 as 1..M, not 0..M-1 */
  free(s);
  CloseFASTA(ffp);

  if ((ffp = OpenFASTA(file2)) == NULL) { fprintf(stderr, "Failed to open file %s\n", file2); exit(1); }
  if (! ReadFASTA(ffp, &s, &name2, &N)) { fprintf(stderr, "No data in file %s\n", file2);     exit(1); }
  seq2 = malloc(sizeof(char) * (N+2));
  strcpy(seq2+1, s);		/* we can now index the residues in seq1 as 1..M, not 0..M-1 */
  free(s);
  CloseFASTA(ffp);

  /* Allocate the dynamic programming matrix, 0..M rows by 0..N columns.
   */
  mx = malloc(sizeof(int *) * (M+1));
  for (i = 0; i <= M; i++)
    mx[i] = malloc(sizeof(int) * (N+1));
  
  /* Pick the midpoint in seq1: the problem asks us to determine
   * what residue seq1[mid] is aligned to in seq2.
   */
  mid = M/2;			

  /* Initialize the dynamic programming matrix. 
   */
  mx[0][0] = 0;
  for (i = 1; i <= M; i++)  mx[i][0] = i * gapcost;
  for (j = 1; j <= N; j++)  mx[0][j] = j * gapcost;

  /* The dynamic programming, global alignment (Needleman/Wunsch) recursion.
   */  
  for (i = 1; i <= M; i++)
    for (j = 1; j <= N; j++)
      {
				              /* score if i,j are aligned */
	if (seq1[i] == seq2[j]) mx[i][j] = mx[i-1][j-1] + mat;
	else                    mx[i][j] = mx[i-1][j-1] + mismat;

	sc = mx[i-1][j] + gapcost;            /* score if i unaligned  */
	if (sc > mx[i][j]) mx[i][j] = sc;

	sc = mx[i][j-1] + gapcost;            /* score if j unaligned */
	if (sc > mx[i][j]) mx[i][j] = sc;
      }

  /* The result is now in mx[M][N].
   */

  /* Trace back the optimal alignment. 
   * This is a simplified traceback that doesn't recover the whole alignment.
   * All we want is to know what residue j was optimally aligned to the middle row mid,
   * (or, mid is unaligned, e.g. inserted).
   * We only have to trace back to the middle row, mid, and figure
   * out what we did there. We record this info in an integer, usedj. When usedj==-1,
   * this means i is unaligned.
   * 
   * Although this traceback is simplified, it's trivial to extend it to
   * a complete alignment recovery.
   */
  i = M; 
  j = N;
  while (i >= mid) 
    {
      sc = mx[i-1][j-1];
      if (seq1[i] == seq2[j]) sc += mat; else sc += mismat;
      if (sc == mx[i][j]) { usedj = j; i--; j--; continue; } /* residue i aligns to j */

      sc = mx[i-1][j] + gapcost;      
      if (sc == mx[i][j]) { usedj = -1; i--; continue; }     /* residue i is inserted */

      sc = mx[i][j-1] + gapcost;      
      if (sc == mx[i][j]) { j--; continue;}                  /* residue j is inserted */
    } 

  /* Output the results.
   */
  printf("Overall alignment score for sequences 1 and 2: %d\n", mx[M][N]);
  if (usedj != -1) 
    printf("Midpoint residue %d [%c] in seq1 optimally aligns to residue %d [%c] in seq2\n", 
	   mid, seq1[mid], usedj, seq2[usedj]);
  else
    printf("Midpoint residue %d [%c] in seq1 is an insertion in the optimal alignment.\n",
	   mid, seq1[mid]);

  /* Free memory and exit.
   */
  free(seq1);
  free(name1);
  free(seq2);
  free(name2);
  for (i = 0; i <= M; i++) free(mx[i]);
  free(mx);
  exit(0);
}


  
