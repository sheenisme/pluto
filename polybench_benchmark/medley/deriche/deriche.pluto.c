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
/* deriche.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "deriche.h"

/* Array initialization. */
static void init_array(int w, int h, DATA_TYPE *alpha,
                       DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
                       DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h))
{
    int i, j;

    *alpha = 0.25; // parameter of the filter

    // input should be between 0 and 1 (grayscale image pixel)
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++)
            imgIn[i][j] = (DATA_TYPE)((313 * i + 991 * j) % 65536) / 65535.0f;
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int w, int h,
                        DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h))

{
    int i, j;

    POLYBENCH_DUMP_START;
    POLYBENCH_DUMP_BEGIN("imgOut");
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++)
        {
            if ((i * h + j) % 20 == 0)
                fprintf(POLYBENCH_DUMP_TARGET, "\n");
            fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i][j]);
        }
    POLYBENCH_DUMP_END("imgOut");
    POLYBENCH_DUMP_FINISH;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Original code provided by Gael Deest */
static void kernel_deriche(int w, int h, DATA_TYPE alpha,
                           DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
                           DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
                           DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
                           DATA_TYPE POLYBENCH_2D(y2, W, H, w, h))
{
    int i, j;
    DATA_TYPE xm1, tm1, ym1, ym2;
    DATA_TYPE xp1, xp2;
    DATA_TYPE tp1, tp2;
    DATA_TYPE yp1, yp2;

    DATA_TYPE k;
    DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
    DATA_TYPE b1, b2, c1, c2;

    DATA_TYPE two = SCALAR_VAL(2.0);
    DATA_TYPE one = SCALAR_VAL(1.0);
    DATA_TYPE zero = SCALAR_VAL(0.0);
    DATA_TYPE neg_two = SCALAR_VAL(-2.0);

    k = (one - EXP_FUN(-alpha)) * (one - EXP_FUN(-alpha)) / (one + two * alpha * EXP_FUN(-alpha) - EXP_FUN(two * alpha));
    a1 = a5 = k;
    a2 = a6 = k * EXP_FUN(-alpha) * (alpha - one);
    a3 = a7 = k * EXP_FUN(-alpha) * (alpha + one);
    a4 = a8 = -k * EXP_FUN(neg_two * alpha);
    b1 = POW_FUN(two, -alpha);
    b2 = -EXP_FUN(neg_two * alpha);
    c1 = c2 = 1;

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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=0;t2<=_PB_W-1;t2++) {
  yp1 = zero;;
  yp2 = zero;;
  xp1 = zero;;
  xp2 = zero;;
  if (_PB_H >= 1) {
    for (t17=0;t17<=_PB_H-1;t17++) {
      for (t18=0;t18<=-t17+_PB_H-1;t18++) {
        for (t19=0;t19<=-t17-t18+_PB_H-1;t19++) {
          for (t20=0;t20<=-t17-t18-t19+_PB_H-1;t20++) {
            for (t21=0;t21<=floord(-t17-t18-t19-t20+_PB_H-1,32);t21++) {
              if ((t17 == 0) && (t18 == 0) && (t19 == 0) && (t20 == 0) && (t21 == 0)) {
                y2[t2][0] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;;
                xp2 = xp1;;
                xp1 = imgIn[t2][0];;
                yp2 = yp1;;
                yp1 = y2[t2][0];;
              }
              if ((t17 == 0) && (t18 == 0) && (t19 == 0) && (t20 == 0)) {
                for (t22=max(1,32*t21);t22<=min(_PB_H-1,32*t21+31);t22++) {
                  y2[t2][t22] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;;
                }
              }
              if ((t17 == 0) && (t18 == 0) && (t19 == 0) && (t20 >= 1) && (t21 == 0)) {
                xp2 = xp1;;
              }
              if ((t17 == 0) && (t18 == 0) && (t19 >= 1) && (t20 == 0) && (t21 == 0)) {
                xp1 = imgIn[t2][t19];;
              }
              if ((t17 == 0) && (t18 >= 1) && (t19 == 0) && (t20 == 0) && (t21 == 0)) {
                yp2 = yp1;;
              }
              if ((t17 >= 1) && (t18 == 0) && (t19 == 0) && (t20 == 0) && (t21 == 0)) {
                yp1 = y2[t2][t17];;
              }
            }
          }
        }
      }
    }
  }
  ym1 = zero;;
  ym2 = zero;;
  xm1 = zero;;
  if (_PB_H >= 1) {
    for (t22=0;t22<=_PB_H-1;t22++) {
      for (t23=0;t23<=-t22+_PB_H-1;t23++) {
        for (t24=0;t24<=-t22-t23+_PB_H-1;t24++) {
          for (t25=0;t25<=floord(-t22-t23-t24+_PB_H-1,32);t25++) {
            if ((t22 == 0) && (t23 == 0) && (t24 == 0) && (t25 == 0)) {
              y1[t2][0] = a1 * imgIn[t2][0] + a2 * xm1 + b1 * ym1 + b2 * ym2;;
              xm1 = imgIn[t2][0];;
              ym2 = ym1;;
              ym1 = y1[t2][0];;
            }
            if ((t22 == 0) && (t23 == 0) && (t24 == 0)) {
              for (t26=max(1,32*t25);t26<=min(_PB_H-1,32*t25+31);t26++) {
                y1[t2][t26] = a1 * imgIn[t2][t26] + a2 * xm1 + b1 * ym1 + b2 * ym2;;
              }
            }
            if ((t22 == 0) && (t23 == 0) && (t24 >= 1) && (t25 == 0)) {
              xm1 = imgIn[t2][t24];;
            }
            if ((t22 == 0) && (t23 >= 1) && (t24 == 0) && (t25 == 0)) {
              ym2 = ym1;;
            }
            if ((t22 >= 1) && (t23 == 0) && (t24 == 0) && (t25 == 0)) {
              ym1 = y1[t2][t22];;
            }
          }
        }
      }
    }
  }
  if (_PB_H >= 1) {
    for (t16=0;t16<=floord(_PB_H-1,32);t16++) {
      for (t17=32*t16;t17<=min(_PB_H-1,32*t16+31);t17++) {
        imgOut[t2][t17] = c1 * (y1[t2][t17] + y2[t2][t17]);;
      }
    }
  }
}
for (t2=0;t2<=_PB_H-1;t2++) {
  tp1 = zero;;
  tp2 = zero;;
  yp1 = zero;;
  yp2 = zero;;
  if (_PB_W >= 1) {
    for (t7=0;t7<=_PB_W-1;t7++) {
      for (t8=0;t8<=-t7+_PB_W-1;t8++) {
        for (t9=0;t9<=-t7-t8+_PB_W-1;t9++) {
          for (t10=0;t10<=-t7-t8-t9+_PB_W-1;t10++) {
            for (t11=0;t11<=floord(-t7-t8-t9-t10+_PB_W-1,32);t11++) {
              if ((t10 == 0) && (t11 == 0) && (t7 == 0) && (t8 == 0) && (t9 == 0)) {
                y2[0][t2] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;;
                tp2 = tp1;;
                tp1 = imgOut[0][t2];;
                yp2 = yp1;;
                yp1 = y2[0][t2];;
              }
              if ((t10 == 0) && (t7 == 0) && (t8 == 0) && (t9 == 0)) {
                for (t12=max(1,32*t11);t12<=min(_PB_W-1,32*t11+31);t12++) {
                  y2[t12][t2] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;;
                }
              }
              if ((t10 >= 1) && (t11 == 0) && (t7 == 0) && (t8 == 0) && (t9 == 0)) {
                tp2 = tp1;;
              }
              if ((t10 == 0) && (t11 == 0) && (t7 == 0) && (t8 == 0) && (t9 >= 1)) {
                tp1 = imgOut[t9][t2];;
              }
              if ((t10 == 0) && (t11 == 0) && (t7 == 0) && (t8 >= 1) && (t9 == 0)) {
                yp2 = yp1;;
              }
              if ((t10 == 0) && (t11 == 0) && (t7 >= 1) && (t8 == 0) && (t9 == 0)) {
                yp1 = y2[t7][t2];;
              }
            }
          }
        }
      }
    }
  }
  tm1 = zero;;
  ym1 = zero;;
  ym2 = zero;;
  if (_PB_W >= 1) {
    for (t12=0;t12<=_PB_W-1;t12++) {
      for (t13=0;t13<=-t12+_PB_W-1;t13++) {
        for (t14=0;t14<=-t12-t13+_PB_W-1;t14++) {
          for (t15=0;t15<=floord(-t12-t13-t14+_PB_W-1,32);t15++) {
            if ((t12 == 0) && (t13 == 0) && (t14 == 0) && (t15 == 0)) {
              y1[0][t2] = a5 * imgOut[0][t2] + a6 * tm1 + b1 * ym1 + b2 * ym2;;
              tm1 = imgOut[0][t2];;
              ym2 = ym1;;
              ym1 = y1[0][t2];;
            }
            if ((t12 == 0) && (t13 == 0) && (t14 == 0)) {
              for (t16=max(1,32*t15);t16<=min(_PB_W-1,32*t15+31);t16++) {
                y1[t16][t2] = a5 * imgOut[t16][t2] + a6 * tm1 + b1 * ym1 + b2 * ym2;;
              }
            }
            if ((t12 == 0) && (t13 == 0) && (t14 >= 1) && (t15 == 0)) {
              tm1 = imgOut[t14][t2];;
            }
            if ((t12 == 0) && (t13 >= 1) && (t14 == 0) && (t15 == 0)) {
              ym2 = ym1;;
            }
            if ((t12 >= 1) && (t13 == 0) && (t14 == 0) && (t15 == 0)) {
              ym1 = y1[t12][t2];;
            }
          }
        }
      }
    }
  }
  for (t6=0;t6<=floord(_PB_W-1,32);t6++) {
    for (t7=32*t6;t7<=min(_PB_W-1,32*t6+31);t7++) {
      imgOut[t7][t2] = c2 * (y1[t7][t2] + y2[t7][t2]);;
    }
  }
}
/* End of CLooG code */
}

int main(int argc, char **argv)
{
    /* Retrieve problem size. */
    int w = W;
    int h = H;

    /* Variable declaration/allocation. */
    DATA_TYPE alpha;
    POLYBENCH_2D_ARRAY_DECL(imgIn, DATA_TYPE, W, H, w, h);
    POLYBENCH_2D_ARRAY_DECL(imgOut, DATA_TYPE, W, H, w, h);
    POLYBENCH_2D_ARRAY_DECL(y1, DATA_TYPE, W, H, w, h);
    POLYBENCH_2D_ARRAY_DECL(y2, DATA_TYPE, W, H, w, h);

    /* Initialize array(s). */
    init_array(w, h, &alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut));

    /* Start timer. */
    polybench_start_instruments;

    /* Run kernel. */
    kernel_deriche(w, h, alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut), POLYBENCH_ARRAY(y1), POLYBENCH_ARRAY(y2));

    /* Stop and print timer. */
    polybench_stop_instruments;
    polybench_print_instruments;

    /* Prevent dead-code elimination. All live-out data must be printed
       by the function call in argument. */
    polybench_prevent_dce(print_array(w, h, POLYBENCH_ARRAY(imgOut)));

    /* Be clean. */
    POLYBENCH_FREE_ARRAY(imgIn);
    POLYBENCH_FREE_ARRAY(imgOut);
    POLYBENCH_FREE_ARRAY(y1);
    POLYBENCH_FREE_ARRAY(y2);

    return 0;
}
