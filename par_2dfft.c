#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "./includes/complex.h"

#define MASTER 0      /* taskid of first task */
#define FROM_MASTER 1 /* setting a message type */
#define FROM_WORKER 2 /* setting a message type */

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

void par_2dfft(complex *sig, complex *f, int n, int rows) {
    complex(*par_f)[n] = f;
    complex(*par_sig)[n] = sig;

    for (int i = 0; i < rows; i++) {
        fft(par_sig[i], par_f[i], 1, n);
    }
}

void help_menu(char *prog_name) {
    printf("Usage: %s [flags]\n", prog_name);
    printf("    -h               prints this usage guide\n");
    printf("    -o <file_name> name for the optional output file\n");
    printf(
        "    -n <number> generate random a matrix of size number x number\n");
}

int main(int argc, char **argv) {
    int numtasks = 1,  /* number of tasks in partition */
        taskid,        /* a task identifier */
        numworkers,    /* number of worker tasks */
        mtype,         /* message type */
        n,             /* size of square matrix */
        rows_per_task, /* number of rows per each task */
        extra_rows,    /* extra rows  */
        offset = 0, ch, rows, total_size;

    double start_time_1, end_time_1, start_time_2, end_time_2, start_time,
        start_time_3, end_time_3, end_time;

    complex *sig, *f, *proccess_sig, *proccess_f;
    FILE *output = NULL;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if (numtasks < 2) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    while ((ch = getopt(argc, argv, "o:n:h")) != -1) {
        switch (ch) {
            case 'n':
                n = atoi(optarg);

                break;

            case 'o':
                output = fopen(optarg, "w");

                break;

            case 'h':
            default:
                help_menu(argv[0]);
                MPI_Finalize();
                exit(1);
                break;
        }
    }

    int offsets[numtasks];
    int all_rows[numtasks];

    rows_per_task = n / numtasks;
    extra_rows = n % numtasks;
    offset = 0;

    for (int dest = 0; dest < numtasks; dest++) {
        rows = (dest < extra_rows) ? rows_per_task + 1 : rows_per_task;
        offsets[dest] = offset;
        all_rows[dest] = rows * n * 2;
        offset += all_rows[dest];
    }

    if (taskid == MASTER) {
        total_size = n * n;

        sig = (complex *)calloc(sizeof(complex), (size_t)total_size);
        f = (complex *)calloc(sizeof(complex), (size_t)total_size);

        for (int i = 0; i < total_size; i++) {
            sig[i].a = rand() % 10;
            sig[i].b = 0;
        }

        transpose(sig, n);  // First transpose
    }

    proccess_sig =
        (complex *)calloc(sizeof(complex), (size_t)all_rows[taskid] / 2);
    proccess_f =
        (complex *)calloc(sizeof(complex), (size_t)all_rows[taskid] / 2);

    if (taskid == MASTER) start_time = now();

    MPI_Scatterv(sig, all_rows, offsets, MPI_DOUBLE, proccess_sig,
                 all_rows[taskid], MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    start_time_1 = now();
    par_2dfft(proccess_sig, proccess_f, n, (all_rows[taskid] / 2) / n);
    end_time_1 = now();

    MPI_Gatherv(proccess_f, all_rows[taskid], MPI_DOUBLE, f, all_rows, offsets,
                MPI_DOUBLE, MASTER,
                MPI_COMM_WORLD);  // proccess_f to f

    if (taskid == MASTER) {
        start_time_3 = now();
        transpose(f, n);
        end_time_3 = now();
    }

    MPI_Scatterv(f, all_rows, offsets, MPI_DOUBLE, proccess_f, all_rows[taskid],
                 MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    start_time_2 = now();
    par_2dfft(proccess_f, proccess_sig, n, (all_rows[taskid] / 2) / n);
    end_time_2 = now();

    MPI_Gatherv(proccess_sig, all_rows[taskid], MPI_DOUBLE, sig, all_rows,
                offsets, MPI_DOUBLE, MASTER,
                MPI_COMM_WORLD);  // process_sig to sig

    if (taskid == MASTER) end_time = now();

    double delta = (end_time_1 - start_time_1) + (end_time_2 - start_time_2);
    double reduced = 0;

    MPI_Reduce(&delta, &reduced, 1, MPI_DOUBLE, MPI_SUM, MASTER,
               MPI_COMM_WORLD);

    if (taskid == MASTER) {
        printf("Time elapsed   : %lf\n",
               (reduced / numtasks) + (end_time_3 - start_time_3));
        printf("Time with comms: %lf\n", end_time - start_time);

        if (output != NULL) {
            complex(*par_sig)[n] = sig;

            for (int i = 0; i < n; i++) {
                fprintf(output, "[ ");
                for (int j = 0; j < n; j++) {
                    print_complex(par_sig[i][j], output);
                }
                fprintf(output, "]\n");
            }
        }
    }

    MPI_Finalize();
    return 0;
}
