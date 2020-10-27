#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "./includes/complex.h"

/*
Description: operates the Discrete Fourier Transform
Input: discrete time signal *sig
       discrete frequency signal *f
       length N
Result: f* constains the Discrete Fourier Transform of *sig
*/
void dft(const complex *sig, complex *f, int N) {
    complex e = euler_formula(2 * -PI / (double)N);
    complex wi = {.a = 1, .b = 0};
    complex wj = {.a = 1, .b = 0};

    for (int i = 0; i < N; i++) {
        f[i] = create_complex(0, 0);

        for (int j = 0; j < N; j++) {
            f[i] = add_complex(f[i], mul_complex(sig[j], wj));
            mul_complex_self(wj, wi);
        }
        mul_complex_self(wi, e);
    }
}

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

void transpose(complex *vec_matrix, int N) {
    complex(*matrix)[N] = vec_matrix;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            swap(&matrix[i][j], &matrix[j][i]);
        }
    }
}

void _2dfft(complex *sig, complex *f, int n) {
    complex(*par_f)[n] = f;
    complex(*par_sig)[n] = sig;

    transpose(sig, n);

    for (int i = 0; i < n; i++) {
        fft(par_sig[i], par_f[i], 1, n);
    }

    transpose(f, n);

    for (int i = 0; i < n; i++) {
        fft(par_f[i], par_sig[i], 1, n);
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

void help_menu(char *prog_name) {
    printf("Usage: %s [flags]\n", prog_name);
    printf("    -h               prints this usage guide\n");
    printf(
        "    -n <number> generate random a matrix of size number x number\n");
}

int main(int argc, char **argv) {
    int n, k, ch;
    complex *sig, *f, *f_seq;
    double start_time, end_time;

    if (argc == 1) {
        help_menu(argv[0]);
        exit(EXIT_SUCCESS);
    }

    while ((ch = getopt(argc, argv, "n:h")) != -1) {
        switch (ch) {
            case 'n':
                n = atoi(optarg);
                break;

            case 'h':
            default:
                help_menu(argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }

    int total_size = n * n;
    sig = (complex *)calloc(sizeof(complex), (size_t)total_size);
    f = (complex *)calloc(sizeof(complex), (size_t)total_size);

    for (int i = 0; i < total_size; i++) {
        sig[i].a = rand() % 10;
        sig[i].b = 0;
    }
    start_time = now();
    _2dfft(sig, f, n);
    end_time = now();

    // for (int i = 0; i < total_size; i++) print_complex(sig[i]);

    printf("Time elapsed: %lf\n", end_time - start_time);

    return 0;
}
