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
/* adi.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "adi.h"

/* Array initialization. */
static void init_array(int n,
                       DATA_TYPE POLYBENCH_2D(u, N, N, n, n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      u[i][j] = (DATA_TYPE)(i + n - j) / n;
    }
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int n,
                        DATA_TYPE POLYBENCH_2D(u, N, N, n, n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      if ((i * n + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
    }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static void kernel_adi(int tsteps, int n,
                       DATA_TYPE POLYBENCH_2D(u, N, N, n, n),
                       DATA_TYPE POLYBENCH_2D(v, N, N, n, n),
                       DATA_TYPE POLYBENCH_2D(p, N, N, n, n),
                       DATA_TYPE POLYBENCH_2D(q, N, N, n, n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;

  DATA_TYPE two = SCALAR_VAL(2.0);
  DATA_TYPE one = SCALAR_VAL(1.0);
  DATA_TYPE zero = SCALAR_VAL(0.0);

  DX = one / (DATA_TYPE)_PB_N;
  DY = one / (DATA_TYPE)_PB_N;
  DT = one / (DATA_TYPE)_PB_TSTEPS;
  B1 = two;
  B2 = one;
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 / two;
  b = one + mul1;
  c = a;
  d = -mul2 / two;
  e = one + mul2;
  f = d;

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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23;
 register int lbv, ubv;
/* Start of CLooG code */
if ((_PB_N >= 3) && (_PB_TSTEPS >= 1)) {
  for (t1=1;t1<=_PB_TSTEPS;t1++) {
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        p[t5][0] = zero;;
        v[_PB_N - 1][t5] = one;;
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      for (t13=0;t13<=floord(_PB_N-2,32);t13++) {
        for (t14=max(1,32*t4);t14<=min(_PB_N-2,32*t4+31);t14++) {
          for (t23=max(1,32*t13);t23<=min(_PB_N-2,32*t13+31);t23++) {
            p[t14][t23] = -c / (a * p[t14][t23 - 1] + b);;
          }
        }
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        v[0][t5] = one;;
        q[t5][0] = v[0][t5];;
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      for (t12=0;t12<=floord(_PB_N-2,32);t12++) {
        for (t13=max(1,32*t4);t13<=min(_PB_N-2,32*t4+31);t13++) {
          for (t21=max(1,32*t12);t21<=min(_PB_N-2,32*t12+31);t21++) {
            q[t13][t21] = (-d * u[t21][t13 - 1] + (one + two * d) * u[t21][t13] - f * u[t21][t13 + 1] - a * q[t13][t21 - 1]) / (a * p[t13][t21 - 1] + b);;
          }
        }
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      lbv=max(1,32*t4);
      ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t5=lbv;t5<=ubv;t5++) {
        u[t5][0] = one;;
        p[t5][0] = zero;;
        u[t5][_PB_N - 1] = one;;
        q[t5][0] = u[t5][0];;
      }
    }
    for (t4=1;t4<=2*_PB_N-4;t4++) {
      if (t4 >= 2) {
        for (t9=max(0,ceild(t4-_PB_N-29,32));t9<=min(floord(_PB_N-2,32),floord(t4-1,32));t9++) {
          lbv=max(max(1,32*t9),t4-_PB_N+2);
          ubv=min(min(_PB_N-2,t4-1),32*t9+31);
#pragma ivdep
#pragma vector always
          for (t10=lbv;t10<=ubv;t10++) {
            q[(t4-t10)][t10] = (-a * v[(t4-t10) - 1][t10] + (one + two * a) * v[(t4-t10)][t10] - c * v[(t4-t10) + 1][t10] - d * q[(t4-t10)][t10 - 1]) / (d * p[(t4-t10)][t10 - 1] + e);;
          }
        }
      }
      if (t4 <= _PB_N-2) {
        for (t11=1;t11<=_PB_N-2;t11++) {
          v[t11][t4] = p[t4][t11] * v[t11 + 1][t4] + q[t4][t11];;
        }
      }
      if (t4 <= _PB_N-2) {
        for (t10=1;t10<=_PB_N-2;t10++) {
          p[t4][t10] = -f / (d * p[t4][t10 - 1] + e);;
        }
      }
    }
    for (t4=0;t4<=floord(_PB_N-2,32);t4++) {
      for (t8=0;t8<=floord(_PB_N-2,32);t8++) {
        for (t9=max(1,32*t4);t9<=min(_PB_N-2,32*t4+31);t9++) {
          for (t13=max(1,32*t8);t13<=min(_PB_N-2,32*t8+31);t13++) {
            u[t9][t13] = p[t9][t13] * u[t9][t13 + 1] + q[t9][t13];;
          }
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array(n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi(tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}
