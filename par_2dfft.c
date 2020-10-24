#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./includes/complex.h"

/*
Description: operates the Fast Fourier Transform
Input: discrete time signal *sig
       discrete frequency signal *f
       length N
Result: f* constains the Discrete Fourier Transform of *sig
*/
void fft(const complex *sig, complex *f, int s, int N) {
    if (N == 1) {
        *f = *sig;
    } else {
        complex ep = euler_formula(2 * -PI / (double)N), ei;
        complex *pi = &ei, *pp = &ep;
        int i, hn = N >> 1;
        fft(sig, f, s << 1, hn);
        fft(sig + s, f + hn, s << 1, hn);
        pi->a = 1;
        pi->b = 0;
        for (i = 0; i < hn; i++) {
            complex even = f[i], *pe = f + i, *po = pe + hn;
            mul_complex_self((*po), (*pi));
            pe->a += po->a;
            pe->b += po->b;
            po->a = even.a - po->a;
            po->b = even.b - po->b;
            mul_complex_self((*pi), (*pp));
        }
    }
}

/*
Output: returns current clock time
*/
double now() {
    const double ONE_BILLION = 1000000000.0;
    struct timespec current_time;

    clock_gettime(CLOCK_REALTIME, &current_time);

    return current_time.tv_sec + (current_time.tv_nsec / ONE_BILLION);
}

void transpose(complex *vec_matrix, int N) {
    complex(*matrix)[N] = vec_matrix;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            swap(&matrix[i][j], &matrix[j][i]);
        }
    }
}

void par_2dfft(complex *sig, complex *f, int N) {
    int n_sqrt = sqrt(N);

    complex(*par_f)[n_sqrt] = f;
    complex(*par_sig)[n_sqrt] = sig;

    printf("\nTranspondo sinal de entrada...\n");
    transpose(sig, n_sqrt);

    printf("\nCalculando FFT das colunas...\n");
    for (int i = 0; i < n_sqrt; i++) {
        fft(par_sig[i], par_f[i], 1, n_sqrt);
    }

    printf("\n## F ##\n");
    for (int i = 0; i < N; i++) print_complex(f[i]);

    printf("\nTranspondo matriz resultante...\n");
    transpose(f, n_sqrt);

    printf("\nCalculando FFT das linhas...\n");
    for (int i = 0; i < n_sqrt; i++) {
        fft(par_f[i], par_sig[i], 1, n_sqrt);
    }

    printf("\n## PAR ##\n");
    for (int i = 0; i < N; i++) print_complex(sig[i]);
}

int main() {
    int n, i, k;
    complex *sig, *f, *f_seq;

    scanf("%d", &k);

    n = 1 << k;
    sig = (complex *)calloc(sizeof(complex), (size_t)n);
    f = (complex *)calloc(sizeof(complex), (size_t)n);

    for (i = 0; i < n; i++) {
        sig[i].a = rand() % 10;
        sig[i].b = 0;
    }

    printf("## Sinal de entrada ##\n");
    for (i = 0; i < n; i++) print_complex(sig[i]);

    par_2dfft(sig, f, n);

    return 0;
}