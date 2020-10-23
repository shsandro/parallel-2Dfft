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
void fft(const complex *sig, complex *f, int s, int N)
{
    int i, hn = N >> 1;
    complex ep = euler_formula(-PI / (double)hn), ei;
    complex *pi = &ei, *pp = &ep;
    if (!hn)
        *f = *sig;
    else
    {
        fft(sig, f, s << 1, hn);
        fft(sig + s, f + hn, s << 1, hn);
        pi->a = 1;
        pi->b = 0;
        for (i = 0; i < hn; i++)
        {
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
double now()
{
    const double ONE_BILLION = 1000000000.0;
    struct timespec current_time;

    clock_gettime(CLOCK_REALTIME, &current_time);

    return current_time.tv_sec + (current_time.tv_nsec / ONE_BILLION);
}

void transpose(complex *vec_matrix, int N)
{
    complex(*matrix)[N] = vec_matrix;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            swap(&matrix[i][j], &matrix[j][i]);
        }
    }
}

void par_fft(complex *sig, complex *f, int N)
{
    int n_sqrt = sqrt(N);

    complex(*par_f)[n_sqrt] = f;
    complex(*par_sig)[n_sqrt] = sig;

    transpose(sig, n_sqrt);
}

int main()
{
    int n, i, k;
    complex *sig, *f;

    scanf("%d", &k);

    n = 1 << k;
    sig = (complex *)malloc(sizeof(complex) * (size_t)n);
    f = (complex *)malloc(sizeof(complex) * (size_t)n);

    for (i = 0; i < n; i++)
    {
        sig[i].a = rand() % 10;
        sig[i].b = 0;
    }

    printf("## Antes ##\n");
    for (i = 0; i < n; i++)
        print_complex(sig[i]);
    printf("#####################\n");
    par_fft(sig, f, n);
    printf("## Depois ##\n");
    for (i = 0; i < n; i++)
        print_complex(sig[i]);

    return 0;
}