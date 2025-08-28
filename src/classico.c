/* *****************************************************************************
   Copyright (c) 2025 I.F.F. dos Santos <ismaellxd@gmail.com>

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
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pvi.h"
/* ---
   Programa escrito durante o ciclo 2024-2025 do PIBIC da UFAL
   para resolver numericamente as equações de movimento de uma rede 1D
   com acoplamento harmonico entre primeiros vizinhos.
   Ultima alteracao: 25 de agosto de 2025
--- */

/* Buffer para a memoria alocada nesta unidade de translacao,
   ele eh inicializado por `preparar_sistema` e aponta para
   uma regiao de memoria liberada por `main`. */
static double * buffer;

#define SIZE_C(x) ((size_t)(x))

typedef double *double_p;

static size_t contar_linhas(FILE *arquivo);

/* As grandezas fisicas de relevante interesse sao declaradas
   como variaveis globais, a fim de que sejam acessiveis a partir
   de qualquer subrotina nessa unidade de translacao. */
static size_t N; /* numero de corpos oscilando */
static double t; /* variavel independente */
static double_p massa, kappa; /* parametros */
static double_p Q, P; /* variaveis dependentes */
static double E; /* energia do sistema */

/* componentes do campo vetorial hamiltoniano */
static double dot_Q(size_t n, double *P);
static double dot_P(size_t n, double *Q);
static double hamiltoniano(void);

static int preparar_sistema(char *nome_arquivo);
static int escrever(double t);

#undef PVI_FAC_ALIQUID
#define PVI_FAC_ALIQUID() {\
   ++i; \
   if(i < i_max) continue;\
   i = 0; \
   if(escrever(t) != 0) break;\
}

int main(int argc, char **argv){
   unsigned i, i_max;
   double *buffer;
   FILE *arquivo;
   int status;

   if(argc < 2){
      fputs("ERRO: Argumentos inv" "\xc3\xa1" "lidos.\n", stderr);
      fprintf(stderr, "%s [arquivo] <tempo final>\n", argv[0]);
      return EXIT_FAILURE;
   }

   status = preparar_sistema(argv[1]);
   if(status != EXIT_SUCCESS) return status;

   /* Resolver numericamente o PVI */
   pvi_dimensio = N;

   t = 0.0;
   pvi_h = (argc > 3 ? atof(argv[3]) : 0.5);
   pvi_finalis = (argc > 2 ? atof(argv[2]) : 10.0);

   i = 0U;
   i_max = (unsigned)(1.0 / (2.0 * pvi_h)); // escreve duas vezes por segundo
   PVI_INTEGRATOR_RUTH4(t, Q, P, dot_Q, dot_P);

   free(buffer);
   return EXIT_SUCCESS;
}

static int preparar_sistema(char *nome_arquivo){
   FILE *arquivo;

   arquivo = fopen(nome_arquivo, "rb");
   if(arquivo == NULL){
      fputs(
         "ERRO: "
         "N" "\xC3\xA3" "o foi poss" "\xC3\xAD" "vel "
         "abrir o arquivo para leitura.\n",
         stderr
      );
      return EXIT_FAILURE;
   }

   N = contar_linhas(arquivo);

   buffer = malloc((4 * N + 3) * sizeof(*buffer));
   if(buffer == NULL){
      fputs(
         "ERRO: "
         "N" "\xC3\xA3" "o h" "\xC3\xA1" " suficiente mem" "\xC3\xB3" "ria.\n",
         stderr
      );
      fclose(arquivo);
      return EXIT_FAILURE;
   }
   massa = buffer;
   kappa = buffer + N + 1;
   Q = buffer + 2*N + 2;
   P = buffer + 3*N + 3;

   /* por hipotese o arquivo deve conter N linhas com quatro colunas cada:
      * a primeira coluna deve conter o valor da massa do corpo;
      * a segunda deve conter o valor da constante de acoplamento harmonico
        entre o corpo e o seguinte vizinho;
      * a terceira o valor inicial do deslocamento do corpo;
      * a quarta o valor inicial do momento linear do corpo. */
   kappa[-1] = 0.0;
   Q[-1] = Q[N] = 0.0;
   for(size_t n = SIZE_C(0); n < N; ++n){
      fscanf(
         arquivo, "%lf %lf %lf %lf",
         massa + n, kappa + n, Q + n, P + n
      );
   }
   fclose(arquivo);

   /* calculo da energia inicial */
   E = hamiltoniano();

   return EXIT_SUCCESS;
}

static int escrever(double t){
   if(fabs(E - hamiltoniano()) > 1.0e-8){
      fprintf(stderr, "A energia n" "\xC3\xA3" "o foi conservada.\n");
      return 1;
   }
   for(size_t n = SIZE_C(0); n < N; ++n){
      fprintf(stdout, "%g %u %g %g\n", t, (unsigned)n, Q[n], P[n]);\
   }
   fprintf(stdout, "\n");
   return 0;
}

static double dot_Q(size_t n, double *P){
   return P[n] / massa[n];
}
static double dot_P(size_t n, double *Q){
   return kappa[n] * (Q[n+1] - Q[n]) - kappa[n-1] * (Q[n] - Q[n-1]);
}

static double square(double x){
   return x*x;
}
static double hamiltoniano(void){
   double H = 0.0;
   for(size_t n = SIZE_C(0); n < N; ++n){
      // energia cinetica
      H += 0.5 * square(P[n]) / massa[n];
      // energia potencial
      H += 0.25 * kappa[n] * square(Q[n+1] - Q[n]);
      H += 0.25 * kappa[n-1] * square(Q[n] - Q[n-1]);
   }
   return H;
}

static size_t contar_linhas(FILE *arquivo){
   int byte;
   long int offset;
   size_t linhas = SIZE_C(1);

   offset = ftell(arquivo);
   fseek(arquivo, 0L, SEEK_SET);

   loop: {
      byte = fgetc(arquivo);
      if(byte != EOF){
         if(byte == (int)'\n') ++linhas;
         goto loop;
      }
   }

   fseek(arquivo, offset, SEEK_SET);

   return linhas;
}
