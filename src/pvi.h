/* *****************************************************************************
   Macro library to implement various methods to solve the
   initial value problem.

   Last modified: April 08, 2024
   *****************************************************************************
   E-mail: ismlxd@gmail.com
   Site: https://ismdamiao.github.io/
   *****************************************************************************
   Copyright (C) 2024 I.F.F. dos Santos

   Permission is hereby granted, free of charge, to any person obtaining a copy 
   of this software and associated documentation files (the “Software”), to 
   deal in the Software without restriction, including without limitation the 
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
   sell copies of the Software, and to permit persons to whom the Software is 
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in 
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
   IN THE SOFTWARE.
***************************************************************************** */
#ifndef PVI_H
#define PVI_H 1

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__STDC__)
# error "Utilize um compilador compativel com o ISO C."
#endif /* __STDC__ */

#include <stddef.h>
#include <stdlib.h>

#define PVI_CORPUS double
#define PVI_FAC_ALIQUID()
#define PVI_ALLOCARE() (PVI_CORPUS*)malloc(pvi_dimensio*sizeof(PVI_CORPUS))

/* For the value od pvi_h I recommend to use numbers that can be written as
   as um of term of the form k / 2 ** n where k and n are
   non-negative integers, for example:
   * 0.875
   * 0.75
   * 0.5
   * 0.25
   * 0.125
   * 0.0625
   * 0.015625
   * 0.00390625
   I also recommend pvi_finalis to be divisible by pvi_h. */

static size_t pvi_dimensio = (size_t)1;
static double pvi_h = 0.25, pvi_finalis = 1.0;

/* ------------------------------------
   Metodos de Runge-Kutta
----------------------------------- */

#define PVI_INTEGRATOR_EULER(t, X, X_punctum) \
{\
   size_t pvi_index;\
   PVI_CORPUS *pvi_inclinatio = NULL;\
\
   pvi_inclinatio = PVI_ALLOCARE();\
\
   while(t < pvi_finalis){\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += pvi_inclinatio[pvi_index] * pvi_h;\
      t += pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio);\
}


#define PVI_INTEGRATOR_RK2(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_h2;\
   PVI_CORPUS \
      *pvi_inclinatio[2] = { NULL, NULL },\
      *pvi_Xaux = NULL;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
\
   pvi_h2 = 0.5 * pvi_h;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
         pvi_Xaux[pvi_index] = \
            (X)[pvi_index] + pvi_inclinatio[0][pvi_index] * pvi_h;\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += pvi_h2 * \
            (pvi_inclinatio[0][pvi_index] + pvi_inclinatio[1][pvi_index]);\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_RK4(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_h2, pvi_h6;\
   double pvi__t_aux;\
   PVI_CORPUS \
      *pvi_inclinatio[4] = { NULL, NULL, NULL, NULL },\
      *pvi_Xaux[2] = {NULL, NULL};\
\
   (void)pvi__t_aux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_Xaux[0] = PVI_ALLOCARE();\
   pvi_Xaux[1] = PVI_ALLOCARE();\
\
   pvi_h2 = 0.5 * pvi_h;\
   pvi_h6 = pvi_h / 6.0;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
         pvi_Xaux[0][pvi_index] = \
            (X)[pvi_index] + pvi_inclinatio[0][pvi_index] * pvi_h2;\
      }\
      pvi__t_aux = t + pvi_h2;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[1][pvi_index] = \
            X_punctum(pvi_index, pvi__t_aux, pvi_Xaux[0]);\
         pvi_Xaux[1][pvi_index] = \
            (X)[pvi_index] + pvi_inclinatio[1][pvi_index] * pvi_h2;\
      }\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[2][pvi_index] = \
            X_punctum(pvi_index, pvi__t_aux, pvi_Xaux[1]);\
         pvi_Xaux[0][pvi_index] = \
            (X)[pvi_index] + pvi_inclinatio[2][pvi_index] * pvi_h;\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux[0]);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += pvi_h6 * (\
            pvi_inclinatio[0][pvi_index] + pvi_inclinatio[3][pvi_index] +\
            2.0 * (pvi_inclinatio[1][pvi_index] + pvi_inclinatio[2][pvi_index])\
          );\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_Xaux[0]);\
   free(pvi_Xaux[1]);\
}



/* ------------------------------------
   Metodos de multipasso
----------------------------------- */

/* https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods */
#define PVI_INTEGRATOR_AB2(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[2], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[2] = { NULL, NULL };\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
\
   pvi_hh[0] = - 0.5 * pvi_h;\
   pvi_hh[1] = 1.5 * pvi_h;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
\
   pvi_finalis1 = pvi_finalis;\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK2(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
}

/* https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods */
#define PVI_INTEGRATOR_AB3(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[3], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[3] = { NULL, NULL, NULL };\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
\
   pvi_hh[0] = pvi_h * 5.0 / 12.0;\
   pvi_hh[1] = pvi_h * (-4.0 / 3.0);\
   pvi_hh[2] = pvi_h * 23.0 / 12.0;\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
}

#define PVI_INTEGRATOR_AB4(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[4], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[4] = { NULL, NULL, NULL, NULL };\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
\
   pvi_hh[0] = pvi_h * (-9.0 / 24.0);\
   pvi_hh[1] = pvi_h * (37.0 / 24.0);\
   pvi_hh[2] = pvi_h * (-59.0 / 24.0);\
   pvi_hh[3] = pvi_h * (55.0 / 24.0);\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
}

#define PVI_INTEGRATOR_AB5(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[5], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[5] = { NULL, NULL, NULL, NULL, NULL };\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_inclinatio[4] = PVI_ALLOCARE();\
\
   pvi_hh[0] = pvi_h * 251.0 / 720.0;\
   pvi_hh[1] = pvi_h * (-1274.0 / 720.0);\
   pvi_hh[2] = pvi_h * 2616.0 / 720.0;\
   pvi_hh[3] = pvi_h * (-2774.0 / 720.0);\
   pvi_hh[4] = pvi_h * 1901.0 / 720.0;\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[4][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[4][pvi_index] * pvi_hh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
         pvi_inclinatio[3][pvi_index] = pvi_inclinatio[4][pvi_index];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_inclinatio[4]);\
}

/* Dr. F.A.B.F. de Moura */
#define PVI_INTEGRATOR_AB10(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[10], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[10] = { \
      NULL, NULL, NULL, NULL, NULL, \
      NULL, NULL, NULL, NULL, NULL };\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_inclinatio[4] = PVI_ALLOCARE();\
   pvi_inclinatio[5] = PVI_ALLOCARE();\
   pvi_inclinatio[6] = PVI_ALLOCARE();\
   pvi_inclinatio[7] = PVI_ALLOCARE();\
   pvi_inclinatio[8] = PVI_ALLOCARE();\
   pvi_inclinatio[9] = PVI_ALLOCARE();\
\
   pvi_hh[0] = pvi_h * (-2082753.0 / 7257600.0);\
   pvi_hh[1] = pvi_h * (20884811.0 / 7257600.0);\
   pvi_hh[2] = pvi_h * (-94307320.0 / 7257600.0);\
   pvi_hh[3] = pvi_h * (252618224.0 / 7257600.0);\
   pvi_hh[4] = pvi_h * (-444772162.0 / 7257600.0);\
   pvi_hh[5] = pvi_h * (538363838.0 / 7257600.0);\
   pvi_hh[6] = pvi_h * (-454661776.0 / 7257600.0);\
   pvi_hh[7] = pvi_h * (265932680.0 / 7257600.0);\
   pvi_hh[8] = pvi_h * (-104995189.0 / 7257600.0);\
   pvi_hh[9] = pvi_h * (30277247.0 / 7257600.0);\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[4][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[5][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[6][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[7][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[8][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[9][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[9][pvi_index] * pvi_hh[9] + \
            pvi_inclinatio[8][pvi_index] * pvi_hh[8] + \
            pvi_inclinatio[7][pvi_index] * pvi_hh[7] + \
            pvi_inclinatio[6][pvi_index] * pvi_hh[6] + \
            pvi_inclinatio[5][pvi_index] * pvi_hh[5] + \
            pvi_inclinatio[4][pvi_index] * pvi_hh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
         pvi_inclinatio[3][pvi_index] = pvi_inclinatio[4][pvi_index];\
         pvi_inclinatio[4][pvi_index] = pvi_inclinatio[5][pvi_index];\
         pvi_inclinatio[5][pvi_index] = pvi_inclinatio[6][pvi_index];\
         pvi_inclinatio[6][pvi_index] = pvi_inclinatio[7][pvi_index];\
         pvi_inclinatio[7][pvi_index] = pvi_inclinatio[8][pvi_index];\
         pvi_inclinatio[8][pvi_index] = pvi_inclinatio[9][pvi_index];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_inclinatio[4]);\
   free(pvi_inclinatio[5]);\
   free(pvi_inclinatio[6]);\
   free(pvi_inclinatio[7]);\
   free(pvi_inclinatio[8]);\
   free(pvi_inclinatio[9]);\
}

/* ------------------------------------
   Metodos previsor-corretor
----------------------------------- */

#define PVI_INTEGRATOR_ABM1(t, X, X_punctum) \
{\
   size_t pvi_index;\
   PVI_CORPUS *pvi_inclinatio = NULL, *pvi_Xaux = NULL;\
\
   pvi_inclinatio = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_Xaux[pvi_index] = \
         (X)[pvi_index] +  pvi_inclinatio[pvi_index] * pvi_h;\
      t += pvi_h;\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += pvi_inclinatio[pvi_index] * pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_ABM2(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[2], pvi_hhh[3], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[2] = { NULL, NULL }, *pvi_Xaux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
   /* Coeficiente de Adans-Bashforth */\
   pvi_hh[0] = pvi_h * (-0.5);\
   pvi_hh[1] = pvi_h * 1.5;\
   /* Coeficiente de Adans-Moulton */\
   pvi_hhh[0] = pvi_h * 0.5;\
   pvi_hhh[1] = pvi_h * 0.5;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
\
   pvi_finalis1 = pvi_finalis;\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK2(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_Xaux[pvi_index] = (X)[pvi_index] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[1][pvi_index] * pvi_hhh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hhh[0];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_ABM3(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[3], pvi_hhh[3], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[3] = { NULL, NULL, NULL }, *pvi_Xaux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
   /* Coeficiente de Adans-Bashforth */\
   pvi_hh[0] = pvi_h * 5.0 / 12.0;\
   pvi_hh[1] = pvi_h * (-4.0 / 3.0);\
   pvi_hh[2] = pvi_h * 23.0 / 12.0;\
   /* Coeficiente de Adans-Moulton */\
   pvi_hhh[0] = pvi_h * (-1.0 / 12.0);\
   pvi_hhh[1] = pvi_h * (8.0 / 12.0);\
   pvi_hhh[2] = pvi_h * (5.0 / 12.0);\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_Xaux[pvi_index] = (X)[pvi_index] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
      }\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[2][pvi_index] * pvi_hhh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hhh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hhh[0];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_ABM4(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[4], pvi_hhh[4], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[4] = { NULL, NULL, NULL, NULL }, *pvi_Xaux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
   /* Coeficiente de Adans-Bashforth */\
   pvi_hh[0] = pvi_h * (-9.0 / 24.0);\
   pvi_hh[1] = pvi_h * (37.0 / 24.0);\
   pvi_hh[2] = pvi_h * (-59.0 / 24.0);\
   pvi_hh[3] = pvi_h * (55.0 / 24.0);\
   /* Coeficiente de Adans-Moulton */\
   pvi_hhh[0] = pvi_h * (1.0 / 24.0);\
   pvi_hhh[1] = pvi_h * (-5.0 / 24.0);\
   pvi_hhh[2] = pvi_h * (19.0 / 24.0);\
   pvi_hhh[3] = pvi_h * (9.0 / 24.0);\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_Xaux[pvi_index] = (X)[pvi_index] + \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
      }\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[3][pvi_index] * pvi_hhh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hhh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hhh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hhh[0];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_ABM5(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[5], pvi_hhh[5], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[5] = { NULL, NULL, NULL, NULL, NULL }, *pvi_Xaux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_inclinatio[4] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
   /* Coeficiente de Adans-Bashforth */\
   pvi_hh[0] = pvi_h * 251.0 / 720.0;\
   pvi_hh[1] = pvi_h * (-1274.0 / 720.0);\
   pvi_hh[2] = pvi_h * 2616.0 / 720.0;\
   pvi_hh[3] = pvi_h * (-2774.0 / 720.0);\
   pvi_hh[4] = pvi_h * 1901.0 / 720.0;\
   /* Coeficiente de Adans-Moulton */\
   pvi_hhh[0] = pvi_h * (-19.0 / 720.0);\
   pvi_hhh[1] = pvi_h * (106.0 / 720.0);\
   pvi_hhh[2] = pvi_h * (-264.0 / 720.0);\
   pvi_hhh[3] = pvi_h * (646.0 / 720.0);\
   pvi_hhh[4] = pvi_h * 251.0 / 720.0;\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[4][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_Xaux[pvi_index] = (X)[pvi_index] + \
            pvi_inclinatio[4][pvi_index] * pvi_hh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
         pvi_inclinatio[3][pvi_index] = pvi_inclinatio[4][pvi_index];\
      }\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[4][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[4][pvi_index] * pvi_hhh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hhh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hhh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hhh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hhh[0];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_inclinatio[4]);\
   free(pvi_Xaux);\
}

#define PVI_INTEGRATOR_ABM10(t, X, X_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[10], pvi_hhh[10], pvi_finalis1;\
   PVI_CORPUS \
      *pvi_inclinatio[10] = { \
      NULL, NULL, NULL, NULL, NULL, \
      NULL, NULL, NULL, NULL, NULL }, *pvi_Xaux;\
\
   pvi_inclinatio[0] = PVI_ALLOCARE();\
   pvi_inclinatio[1] = PVI_ALLOCARE();\
   pvi_inclinatio[2] = PVI_ALLOCARE();\
   pvi_inclinatio[3] = PVI_ALLOCARE();\
   pvi_inclinatio[4] = PVI_ALLOCARE();\
   pvi_inclinatio[5] = PVI_ALLOCARE();\
   pvi_inclinatio[6] = PVI_ALLOCARE();\
   pvi_inclinatio[7] = PVI_ALLOCARE();\
   pvi_inclinatio[8] = PVI_ALLOCARE();\
   pvi_inclinatio[9] = PVI_ALLOCARE();\
   pvi_Xaux = PVI_ALLOCARE();\
   /* Coeficiente de Adans-Bashforth */\
   pvi_hh[0] = pvi_h * (-2082753.0 / 7257600.0);\
   pvi_hh[1] = pvi_h * (20884811.0 / 7257600.0);\
   pvi_hh[2] = pvi_h * (-94307320.0 / 7257600.0);\
   pvi_hh[3] = pvi_h * (252618224.0 / 7257600.0);\
   pvi_hh[4] = pvi_h * (-444772162.0 / 7257600.0);\
   pvi_hh[5] = pvi_h * (538363838.0 / 7257600.0);\
   pvi_hh[6] = pvi_h * (-454661776.0 / 7257600.0);\
   pvi_hh[7] = pvi_h * (265932680.0 / 7257600.0);\
   pvi_hh[8] = pvi_h * (-104995189.0 / 7257600.0);\
   pvi_hh[9] = pvi_h * (30277247.0 / 7257600.0);\
   /* Coeficiente de Adans-Moulton */\
   pvi_hhh[0] = pvi_h * (57281.0 / 7257600.0);\
   pvi_hhh[1] = pvi_h * (-583435.0 / 7257600.0);\
   pvi_hhh[2] = pvi_h * (2687864.0 / 7257600.0);\
   pvi_hhh[3] = pvi_h * (-7394032.0 / 7257600.0);\
   pvi_hhh[4] = pvi_h * (13510082.0 / 7257600.0);\
   pvi_hhh[5] = pvi_h * (-17283646.0 / 7257600.0);\
   pvi_hhh[6] = pvi_h * (16002320.0 / 7257600.0);\
   pvi_hhh[7] = pvi_h * (-11271304.0 / 7257600.0);\
   pvi_hhh[8] = pvi_h * (9449717.0 / 7257600.0);\
   pvi_hhh[9] = pvi_h * (2082753.0 / 7257600.0);\
\
   pvi_finalis1 = pvi_finalis;\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[0][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[1][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[2][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[3][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[4][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[5][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[6][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[7][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
\
   for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
      pvi_inclinatio[8][pvi_index] = X_punctum(pvi_index, t, X);\
   pvi_finalis = t + pvi_h;\
   PVI_INTEGRATOR_RK4(t, X, X_punctum);\
   pvi_finalis = pvi_finalis1;\
\
   while(t < pvi_finalis){\
      /* ----- Preditor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[9][pvi_index] = X_punctum(pvi_index, t, X);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_Xaux[pvi_index] = (X)[pvi_index] + \
            pvi_inclinatio[9][pvi_index] * pvi_hh[9] + \
            pvi_inclinatio[8][pvi_index] * pvi_hh[8] + \
            pvi_inclinatio[7][pvi_index] * pvi_hh[7] + \
            pvi_inclinatio[6][pvi_index] * pvi_hh[6] + \
            pvi_inclinatio[5][pvi_index] * pvi_hh[5] + \
            pvi_inclinatio[4][pvi_index] * pvi_hh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hh[0];\
      }\
      t += pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         pvi_inclinatio[0][pvi_index] = pvi_inclinatio[1][pvi_index];\
         pvi_inclinatio[1][pvi_index] = pvi_inclinatio[2][pvi_index];\
         pvi_inclinatio[2][pvi_index] = pvi_inclinatio[3][pvi_index];\
         pvi_inclinatio[3][pvi_index] = pvi_inclinatio[4][pvi_index];\
         pvi_inclinatio[4][pvi_index] = pvi_inclinatio[5][pvi_index];\
         pvi_inclinatio[5][pvi_index] = pvi_inclinatio[6][pvi_index];\
         pvi_inclinatio[6][pvi_index] = pvi_inclinatio[7][pvi_index];\
         pvi_inclinatio[7][pvi_index] = pvi_inclinatio[8][pvi_index];\
         pvi_inclinatio[8][pvi_index] = pvi_inclinatio[9][pvi_index];\
      }\
      /* ----- Corretor ----- */ \
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         pvi_inclinatio[9][pvi_index] = X_punctum(pvi_index, t, pvi_Xaux);\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index){\
         (X)[pvi_index] += \
            pvi_inclinatio[9][pvi_index] * pvi_hhh[9] + \
            pvi_inclinatio[8][pvi_index] * pvi_hhh[8] + \
            pvi_inclinatio[7][pvi_index] * pvi_hhh[7] + \
            pvi_inclinatio[6][pvi_index] * pvi_hhh[6] + \
            pvi_inclinatio[5][pvi_index] * pvi_hhh[5] + \
            pvi_inclinatio[4][pvi_index] * pvi_hhh[4] + \
            pvi_inclinatio[3][pvi_index] * pvi_hhh[3] + \
            pvi_inclinatio[2][pvi_index] * pvi_hhh[2] + \
            pvi_inclinatio[1][pvi_index] * pvi_hhh[1] + \
            pvi_inclinatio[0][pvi_index] * pvi_hhh[0];\
      }\
      PVI_FAC_ALIQUID();\
   }\
   free(pvi_inclinatio[0]);\
   free(pvi_inclinatio[1]);\
   free(pvi_inclinatio[2]);\
   free(pvi_inclinatio[3]);\
   free(pvi_inclinatio[4]);\
   free(pvi_inclinatio[5]);\
   free(pvi_inclinatio[6]);\
   free(pvi_inclinatio[7]);\
   free(pvi_inclinatio[8]);\
   free(pvi_inclinatio[9]);\
   free(pvi_Xaux);\
}

/* ------------------------------------
   Metodos simpleticos
----------------------------------- */

#define PVI_INTEGRATOR_EULER_S(t, X, Y, X_punctum, Y_punctum) \
{\
   size_t pvi_index;\
   while(t < pvi_finalis){\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_h;\
      t += pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
}

#define PVI_INTEGRATOR_VERLET(t, X, Y, X_punctum, Y_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh; \
   pvi_hh = pvi_h * 0.5;\
\
   while(t < pvi_finalis){\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_h;\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh;\
      t += pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
}

#define PVI_INTEGRATOR_RUTH3(t, X, Y, X_punctum, Y_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[6]; \
   pvi_hh[0] = pvi_h * (7.0 / 24.0);\
   pvi_hh[1] = pvi_h * (2.0 / 3.0);\
   pvi_hh[2] = pvi_h * 0.75;\
   pvi_hh[3] = pvi_h * (-2.0 / 3.0);\
   pvi_hh[4] = pvi_h * (-1.0 / 24.0);\
   pvi_hh[5] = pvi_h;\
\
   while(t < pvi_finalis){\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[0];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[1];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[2];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[3];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[4];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[5];\
      t += pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
}

#define PVI_RAIZ_CUBICA_2 1.25992104989

#define PVI_INTEGRATOR_RUTH4(t, X, Y, X_punctum, Y_punctum) \
{\
   size_t pvi_index;\
   double pvi_hh[4]; \
   pvi_hh[0] = pvi_h * (0.5 / (2.0 - PVI_RAIZ_CUBICA_2));\
   pvi_hh[1] = pvi_h * (1.0 / (2.0 - PVI_RAIZ_CUBICA_2));\
   pvi_hh[2] = pvi_h * ((1.0 - PVI_RAIZ_CUBICA_2) * 0.5 / (2.0 - PVI_RAIZ_CUBICA_2));\
   pvi_hh[3] = pvi_h * (-PVI_RAIZ_CUBICA_2 / (2.0 - PVI_RAIZ_CUBICA_2));\
\
   while(t < pvi_finalis){\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[0];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[1];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[2];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[3];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[2];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (Y)[pvi_index] += Y_punctum(pvi_index, X) * pvi_hh[1];\
      for(pvi_index = (size_t)0; pvi_index < pvi_dimensio; ++pvi_index)\
         (X)[pvi_index] += X_punctum(pvi_index, Y) * pvi_hh[0];\
      t += pvi_h;\
      PVI_FAC_ALIQUID();\
   }\
}

#ifdef __cplusplus
}
#endif

#endif /* PVI_H */
