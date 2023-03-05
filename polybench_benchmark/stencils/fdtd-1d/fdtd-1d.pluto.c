#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/* fdtd-1d.c: this file is part of PolyBench/C, which is extend by sheen song*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-1d.h"

/* Array initialization. */
static void init_array(int tsteps, int n, DATA_TYPE POLYBENCH_1D(H, N, n),
                       DATA_TYPE POLYBENCH_1D(E, N, n))
{
    int i;

    for (i = 0; i < n; i++)
    {
        H[i] = ((DATA_TYPE)i) / n;
        E[i] = ((DATA_TYPE)i) / n;
    }
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int tsteps, int n, DATA_TYPE POLYBENCH_1D(H, N, n))
{
    int i;

    POLYBENCH_DUMP_START;
    POLYBENCH_DUMP_BEGIN("H");
    for (i = 0; i < n; i++)
    {
        if (i % 20 == 0)
            fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, H[i]);
    }
    POLYBENCH_DUMP_END("H");
    POLYBENCH_DUMP_FINISH;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static void kernel_fdtd_1d(int tsteps, int n, DATA_TYPE POLYBENCH_1D(H, N, n),
                           DATA_TYPE POLYBENCH_1D(E, N, n))
{
    int t, i;
    DATA_TYPE coeff1 = SCALAR_VAL(0.5);
    DATA_TYPE coeff2 = SCALAR_VAL(0.7);

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
  int t1, t2, t3, t4, t5;
 register int lbv, ubv;
/* Start of CLooG code */
if ((_PB_N >= 1) && (_PB_TSTEPS >= 1)) {
  for (t1=0;t1<=floord(_PB_TSTEPS-1,32);t1++) {
    for (t2=t1;t2<=min(floord(_PB_TSTEPS+_PB_N-1,32),floord(32*t1+_PB_N+31,32));t2++) {
      if ((_PB_N >= 2) && (t1 <= floord(32*t2-_PB_N,32))) {
        H[(_PB_N-1)] = H[(_PB_N-1)] - coeff2 * (E[(_PB_N-1) + 1] - E[(_PB_N-1)]);;
      }
      if (_PB_N == 1) {
        for (t3=max(32*t1,32*t2-1);t3<=min(min(_PB_TSTEPS-1,32*t1+31),32*t2+30);t3++) {
          H[0] = H[0] - coeff2 * (E[0 + 1] - E[0]);;
        }
      }
      if (_PB_N >= 2) {
        for (t3=max(32*t1,32*t2-_PB_N+1);t3<=min(min(_PB_TSTEPS-1,32*t1+31),32*t2-_PB_N+31);t3++) {
          for (t4=max(32*t2,t3+1);t4<=t3+_PB_N-1;t4++) {
            E[(-t3+t4)] = E[(-t3+t4)] - coeff1 * (H[(-t3+t4)] - H[(-t3+t4) - 1]);;
            H[(-t3+t4-1)] = H[(-t3+t4-1)] - coeff2 * (E[(-t3+t4-1) + 1] - E[(-t3+t4-1)]);;
          }
          H[(_PB_N-1)] = H[(_PB_N-1)] - coeff2 * (E[(_PB_N-1) + 1] - E[(_PB_N-1)]);;
        }
      }
      for (t3=max(32*t1,32*t2-_PB_N+32);t3<=min(min(_PB_TSTEPS-1,32*t1+31),32*t2+30);t3++) {
        for (t4=max(32*t2,t3+1);t4<=32*t2+31;t4++) {
          E[(-t3+t4)] = E[(-t3+t4)] - coeff1 * (H[(-t3+t4)] - H[(-t3+t4) - 1]);;
          H[(-t3+t4-1)] = H[(-t3+t4-1)] - coeff2 * (E[(-t3+t4-1) + 1] - E[(-t3+t4-1)]);;
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
    POLYBENCH_1D_ARRAY_DECL(H, DATA_TYPE, N, n);
    POLYBENCH_1D_ARRAY_DECL(E, DATA_TYPE, N, n);

    /* Initialize array(s). */
    init_array(tsteps, n, POLYBENCH_ARRAY(H), POLYBENCH_ARRAY(E));

    /* Start timer. */
    polybench_start_instruments;

    /* Run kernel. */
    kernel_fdtd_1d(tsteps, n, POLYBENCH_ARRAY(H), POLYBENCH_ARRAY(E));

    /* Stop and print timer. */
    polybench_stop_instruments;
    polybench_print_instruments;

    /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
    polybench_prevent_dce(print_array(tsteps, n, POLYBENCH_ARRAY(H)));

    /* Be clean. */
    POLYBENCH_FREE_ARRAY(H);
    POLYBENCH_FREE_ARRAY(E);

    return 0;
}
