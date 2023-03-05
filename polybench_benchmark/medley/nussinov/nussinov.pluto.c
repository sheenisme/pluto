#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* nussinov.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "nussinov.h"

/* RNA bases represented as chars, range is [0,3] */
typedef char base;

#define match(b1, b2) (((b1) + (b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)

/* Array initialization. */
static void init_array(int n,
                       base POLYBENCH_1D(seq, N, n),
                       DATA_TYPE POLYBENCH_2D(table, N, N, n, n))
{
   int i, j;

   // base is AGCT/0..3
   for (i = 0; i < n; i++)
   {
      seq[i] = (base)((i + 1) % 4);
   }

   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         table[i][j] = 0;
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int n,
                        DATA_TYPE POLYBENCH_2D(table, N, N, n, n))

{
   int i, j;
   int t = 0;

   POLYBENCH_DUMP_START;
   POLYBENCH_DUMP_BEGIN("table");
   for (i = 0; i < n; i++)
   {
      for (j = i; j < n; j++)
      {
         if (t % 20 == 0)
            fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, table[i][j]);
         t++;
      }
   }
   POLYBENCH_DUMP_END("table");
   POLYBENCH_DUMP_FINISH;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/*
  Original version by Dave Wonnacott at Haverford College <davew@cs.haverford.edu>,
  with help from Allison Lake, Ting Zhou, and Tian Jin,
  based on algorithm by Nussinov, described in Allison Lake's senior thesis.
*/
static void kernel_nussinov(int n, base POLYBENCH_1D(seq, N, n),
                            DATA_TYPE POLYBENCH_2D(table, N, N, n, n))
{
   int i, j, k;

/* Copyright (C) 1991-2020 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 register int lbv, ubv;
/* Start of CLooG code */
if (_PB_N >= 2) {
  for (t1=1;t1<=_PB_N-1;t1++) {
    for (t8=0;t8<=floord(t1-1,32);t8++) {
      lbv=32*t8;
      ubv=min(t1-1,32*t8+31);
#pragma ivdep
#pragma vector always
      for (t9=lbv;t9<=ubv;t9++) {
        table[t9][t1] = max_score(table[t9][t1], table[t9][t1 - 1]);;
      }
    }
    for (t3=0;t3<=floord(t1-1,32);t3++) {
      for (t4=t3;t4<=floord(t1-1,32);t4++) {
        if ((t1 >= 34) && (t3 == 0) && (t4 == 0)) {
          table[0][t1] = max_score(table[0][t1], table[0 + 1][t1]);;
          table[0][t1] = max_score(table[0][t1], table[0 + 1][t1 - 1] + match(seq[0], seq[t1]));;
          for (t7=1;t7<=t1-1;t7++) {
            table[t7][t1] = max_score(table[t7][t1], table[t7 + 1][t1]);;
          }
          for (t6=1;t6<=31;t6++) {
            table[0][t1] = max_score(table[0][t1], table[0][t6] + table[t6 + 1][t1]);;
            table[t6][t1] = max_score(table[t6][t1], table[t6 + 1][t1 - 1] + match(seq[t6], seq[t1]));;
          }
          for (t6=32;t6<=t1-2;t6++) {
            table[t6][t1] = max_score(table[t6][t1], table[t6 + 1][t1 - 1] + match(seq[t6], seq[t1]));;
          }
        }
        if ((t1 >= 3) && (t1 <= 32) && (t3 == 0) && (t4 == 0)) {
          table[0][t1] = max_score(table[0][t1], table[0 + 1][t1]);;
          table[0][t1] = max_score(table[0][t1], table[0 + 1][t1 - 1] + match(seq[0], seq[t1]));;
          for (t7=1;t7<=t1-1;t7++) {
            table[t7][t1] = max_score(table[t7][t1], table[t7 + 1][t1]);;
          }
          for (t6=1;t6<=t1-2;t6++) {
            table[0][t1] = max_score(table[0][t1], table[0][t6] + table[t6 + 1][t1]);;
            table[t6][t1] = max_score(table[t6][t1], table[t6 + 1][t1 - 1] + match(seq[t6], seq[t1]));;
          }
          table[0][t1] = max_score(table[0][t1], table[0][(t1-1)] + table[(t1-1) + 1][t1]);;
        }
        if ((t1 == 33) && (t3 == 0) && (t4 == 0)) {
          table[0][33] = max_score(table[0][33], table[0 + 1][33]);;
          table[0][33] = max_score(table[0][33], table[0 + 1][33 - 1] + match(seq[0], seq[33]));;
          for (t7=1;t7<=32;t7++) {
            table[t7][33] = max_score(table[t7][33], table[t7 + 1][33]);;
          }
          for (t6=1;t6<=31;t6++) {
            table[0][33] = max_score(table[0][33], table[0][t6] + table[t6 + 1][33]);;
            table[t6][33] = max_score(table[t6][33], table[t6 + 1][33 - 1] + match(seq[t6], seq[33]));;
          }
        }
        if ((t1 == 2) && (t3 == 0) && (t4 == 0)) {
          table[0][2] = max_score(table[0][2], table[0 + 1][2]);;
          table[0][2] = max_score(table[0][2], table[0 + 1][2 - 1] + match(seq[0], seq[2]));;
          table[1][2] = max_score(table[1][2], table[1 + 1][2]);;
          table[0][2] = max_score(table[0][2], table[0][1] + table[1 + 1][2]);;
        }
        if ((t1 == 1) && (t3 == 0) && (t4 == 0)) {
          table[0][1] = max_score(table[0][1], table[0 + 1][1]);;
        }
        if ((t3 == 0) && (t4 >= 1)) {
          for (t6=32*t4;t6<=min(t1-1,32*t4+31);t6++) {
            table[0][t1] = max_score(table[0][t1], table[0][t6] + table[t6 + 1][t1]);;
          }
        }
        for (t5=max(1,32*t3);t5<=min(min(t1-2,32*t3+31),32*t4+30);t5++) {
          for (t6=max(32*t4,t5+1);t6<=min(t1-1,32*t4+31);t6++) {
            table[t5][t1] = max_score(table[t5][t1], table[t5][t6] + table[t6 + 1][t1]);;
          }
        }
        if ((t3 == 0) && (t4 == 0)) {
          table[(t1-1)][t1] = max_score(table[(t1-1)][t1], table[(t1-1) + 1][t1 - 1]);;
        }
      }
    }
  }
}
/* End of CLooG code */
}

int main(int argc, char **argv)
{
   /* Retrieve problem size. */
   int n = N;

   /* Variable declaration/allocation. */
   POLYBENCH_1D_ARRAY_DECL(seq, base, N, n);
   POLYBENCH_2D_ARRAY_DECL(table, DATA_TYPE, N, N, n, n);

   /* Initialize array(s). */
   init_array(n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

   /* Start timer. */
   polybench_start_instruments;

   /* Run kernel. */
   kernel_nussinov(n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

   /* Stop and print timer. */
   polybench_stop_instruments;
   polybench_print_instruments;

   /* Prevent dead-code elimination. All live-out data must be printed
      by the function call in argument. */
   polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(table)));

   /* Be clean. */
   POLYBENCH_FREE_ARRAY(seq);
   POLYBENCH_FREE_ARRAY(table);

   return 0;
}
