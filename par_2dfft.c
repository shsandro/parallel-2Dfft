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
        offset, ch, rows;

    double start_time_1, end_time_1, start_time_2, end_time_2, start_time_3,
        end_time_3, start_time, end_time;

    double delta;
    double reduced;
    double delta_comms;
    double reduced_comms;

    complex *sig, *f, *f_seq;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if (numtasks < 2) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    /* ******************** master task ************************ */
    if (taskid == MASTER) {
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
                    MPI_Finalize();
                    exit(1);
                    break;
            }
        }

        int total_size = n * n;
        int offsets[numtasks];
        int all_rows[numtasks];

        sig = (complex *)calloc(sizeof(complex), (size_t)total_size);
        f = (complex *)calloc(sizeof(complex), (size_t)total_size);

        for (int i = 0; i < total_size; i++) {
            sig[i].a = rand() % 10;
            sig[i].b = 0;
        }

        rows_per_task = n / numtasks;
        extra_rows = n % numtasks;
        offset = rows_per_task * n;

        transpose(sig, n);

        start_time = now();
        start_time_1 = now();
        par_2dfft(sig, f, n, rows_per_task);
        end_time_1 = now();

        for (int dest = 1; dest < numtasks; dest++) {
            rows = (dest <= extra_rows) ? rows_per_task + 1 : rows_per_task;
            MPI_Send(&n, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(sig + offset, rows * n * 2, MPI_DOUBLE, dest, FROM_MASTER,
                     MPI_COMM_WORLD);
            offsets[dest] = offset;
            all_rows[dest] = rows;
            offset += (rows * n);
        }

        for (int source = 1; source < numtasks; source++) {
            MPI_Recv(f + offsets[source], all_rows[source] * n * 2, MPI_DOUBLE,
                     source, FROM_WORKER, MPI_COMM_WORLD, &status);
        }

        start_time_3 = now();
        transpose(f, n);
        end_time_3 = now();

        start_time_2 = now();
        par_2dfft(f, sig, n, rows_per_task);
        end_time_2 = now();

        for (int dest = 1; dest < numtasks; dest++) {
            MPI_Send(f + offsets[dest], all_rows[dest] * n * 2, MPI_DOUBLE,
                     dest, FROM_MASTER, MPI_COMM_WORLD);
        }

        for (int source = 1; source < numtasks; source++) {
            MPI_Recv(sig + offsets[source], all_rows[source] * n * 2,
                     MPI_DOUBLE, source, FROM_WORKER, MPI_COMM_WORLD, &status);
        }
        end_time = now();

        delta = (end_time_1 - start_time_1) + (end_time_2 - start_time_2);
        reduced = 0;

        MPI_Reduce(&delta, &reduced, 1, MPI_DOUBLE, MPI_SUM, MASTER,
                   MPI_COMM_WORLD);

        delta_comms = end_time - start_time;
        reduced_comms = 0;

        MPI_Reduce(&delta_comms, &reduced_comms, 1, MPI_DOUBLE, MPI_SUM, MASTER,
                   MPI_COMM_WORLD);

        // for (int i = 0; i < total_size; i++) print_complex(sig[i]);
        // printf("Time elapsed: %lf\n",
        //        (end_time_1 - start_time_1) + (end_time_2 - start_time_2));
        // printf("Time elapsed with comms: %lf\n", (end_time - start_time));
    } else {
        /* ******************** worker task ************************ */
        MPI_Recv(&n, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);

        MPI_Recv(&offset, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD,
                 &status);

        MPI_Recv(&rows, 1, MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD,
                 &status);

        sig = (complex *)calloc(sizeof(complex), (size_t)n * rows);
        f = (complex *)calloc(sizeof(complex), (size_t)n * rows);

        MPI_Recv(sig, rows * n * 2, MPI_DOUBLE, MASTER, FROM_MASTER,
                 MPI_COMM_WORLD, &status);

        start_time = now();
        start_time_1 = now();
        par_2dfft(sig, f, n, rows);
        end_time_1 = now();

        MPI_Send(f, rows * n * 2, MPI_DOUBLE, MASTER, FROM_WORKER,
                 MPI_COMM_WORLD);

        MPI_Recv(f, rows * n * 2, MPI_DOUBLE, MASTER, FROM_MASTER,
                 MPI_COMM_WORLD, &status);

        start_time_2 = now();
        par_2dfft(f, sig, n, rows);
        end_time_2 = now();

        MPI_Send(sig, rows * n * 2, MPI_DOUBLE, MASTER, FROM_WORKER,
                 MPI_COMM_WORLD);

        end_time = now();

        delta = (end_time_1 - start_time_1) + (end_time_2 - start_time_2);
        reduced = 0;

        MPI_Reduce(&delta, &reduced, 1, MPI_DOUBLE, MPI_SUM, MASTER,
                   MPI_COMM_WORLD);

        delta_comms = end_time - start_time;
        reduced_comms = 0;

        MPI_Reduce(&delta_comms, &reduced_comms, 1, MPI_DOUBLE, MPI_SUM, MASTER,
                   MPI_COMM_WORLD);
    }

    if (taskid == MASTER) {
        printf("Time elapsed   : %lf\n",
               (reduced / numtasks) + end_time_3 - start_time_3);
        printf("Time with comms: %lf\n", reduced_comms / numtasks);
    }

    MPI_Finalize();
    return 0;
}
