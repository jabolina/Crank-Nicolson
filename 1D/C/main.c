#include <stdio.h>
#include <stdlib.h>
#include <math.h>


long double analytic_solution(long double x, long double t) {
    return sin(x) * exp(-t);
}

void initial_value(long double* arr, long double dx, int N) {
    int i;

    for (i = 0; i<N ; i++)
        arr[i] = analytic_solution((i * dx), 0.0);
}

long double** alloc_matrix(int N, long double s) {
    int i, j;
    long double** m = (long double**)malloc(N * sizeof(long double*));

    for(i=0; i<N; i++) {
        m[i] = (long double*)malloc(N * sizeof(long double));
        for(j=0;j<N;j++) {
            if(j == i) {
                m[i][j] = (1 + 2*s);
            }
            else if((j == i-1) || (j == i+1)) {
                m[i][j] = -s;
            }
            else {
                m[i][j] = 0;
            }
        }
    }

    return m;
}

long double** free_matrix(long double** A, int N) {
    int i;

    for(i=0; i<N; i++)
        free(A[i]);

    free(A);
    return A;
}

long double** rearrange(long double** m, int N) {

    int i;

    for(i=0; i<N; i++) {
        m[0][i] = 0;
        m[N-1][i] = 0;
    }
    m[0][0] = 1;
    m[N-1][N-1] = 1;

    return m;
}

void fill_values(long double* an, int N, long double dx, long double ta) {
    int i;

    for(i=0; i<N; i++)
        an[i] = analytic_solution((i * dx), ta);
}

void copy(const long double* v, long double* T, int N) {
    int i;
    for(i=0; i<N; i++)
        T[i] = v[i];
}

void right_hand_updt(long double* nu, int N, long double s) {
    int i;
    long double* x = (long double*)calloc((size_t) N, sizeof(long double));

    x[0] = 0.0;
    nu[0] = 0.0;

    for(i=1; i<(N - 1); i++)
        x[i] = (s * nu[i-1]) + ((1.0 - 2.0*s) * nu[i]) + (s * nu[i+1]);

    x[N - 1] = 0.0;
    nu[N - 1] = 0.0;

    copy(x, nu, N);
    free(x);
}

void gauss_seidel(long double** A, long double* b, int N) {
    int i, j, iter = 5 * N;
    long double aux;
    long double* x = (long double*)malloc(N * sizeof(long double));

    for (i=0; i<N; i++) x[i] = 1.0;

    while (iter > 0) {
        for(i=0; i<N; i++) {
            aux = 0;
            for(j=0; j<N; j++) {
                if(j != i)
                    aux += A[i][j] * x[j];
            }
            x[i] = (b[i] - aux) / A[i][i];
        }
        iter--;
    }

    copy(x, b, N);
    free(x);
}

long double infinity_norm(const long double* nu, const long double* an, int N) {
    int i;
    long double* l = (long double*)calloc(N, sizeof(long double));
    long double greatest;

    for(i=0; i<N; i++) {
        l[i] = sqrt((nu[i] - an[i]) * (nu[i] - an[i]));
    }

    greatest = l[0];
    for(i=1; i<N; i++)
        (l[i] > greatest) ? greatest = l[i]: greatest;

    return greatest;
}

void save_file(const long double* nu, const long double* an, int N, long double dx) {
    int i;

    FILE* fa = fopen("g_a.dat", "a");
    FILE* fn = fopen("g_n.dat", "a");

    for(i=0; i<N; i++) {
        fprintf(fa, "%Lf\t\t%Lf\n", (i * dx), an[i]);
        fprintf(fn, "%Lf\t\t%Lf\n", (i * dx), nu[i]);
    }

    fprintf(fa, "\n");
    fprintf(fn, "\n");

    fclose(fa);
    fclose(fn);
}

void collect_data(unsigned int *N, long double* dt, int* ts) {

    printf("Cells: ");
    scanf("%u", N);

    printf("Time steps: ");
    scanf("%d", ts);

    printf("Delta t: ");
    scanf("%Lf", dt);

    system("clear");
}


int main() {

    remove("infinity_norm.txt");

    /*
     * N -> Nodes in the mesh
     */
    unsigned int N;
    int o;

    for(o=2; o<=8; o++) {
        N = (unsigned int) pow(2, o);

        /*
         * L -> Length, [0, 2pi]
         * dx -> Delta x
         * dt -> Delta t
         * s -> auxiliary variable
         */
        long double L, dx, dt, s;

        // ts -> Time steps
        int ts;

        long double **A;

        system("clear");
        remove("g_a.dat");
        remove("g_n.dat");

        //collect_data(&N, &dt, &ts);

        ts = 100;
        dt = 0.001;
        L = 2 * M_PI;
        dx = L / (N - 1);
        s = dt / (2.0 * (dx * dx));

        system("clear");

        A = alloc_matrix(N, s);
        A = rearrange(A, N);

        ////////////////////////////// S O L V I N G ///////////////////////////////

        int i;

        /*
         * nu -> Numeric solution
         * x  -> Auxiliary
         * an -> Analytic solution
         */
        long double *nu, *x, *an, *error;

        nu = (long double *) calloc(N, sizeof(long double));
        an = (long double *) calloc(N, sizeof(long double));
        x = (long double *) calloc(N, sizeof(long double));
        error = (long double *) calloc((size_t) ts, sizeof(long double));

        initial_value(nu, dx, N);
        initial_value(an, dx, N);

        for (i = 1; i < ts; i++) {
            save_file(nu, an, N, dx);
            error[i - 1] = infinity_norm(nu, an, N);

            fill_values(an, N, dx, i * dt);
            right_hand_updt(nu, N, s);
            gauss_seidel(A, nu, N);
        }

        ////////////////////////////////////////////////////////////////////////////

        long double greatest = error[0];
        for (i = 1; i < ts; i++)
            (error[i] > greatest) ? greatest = error[i] : greatest;

        FILE *fe = fopen("infinity_norm.txt", "a");
        fprintf(fe, "%Lf\n", greatest);
        fclose(fe);

        ////////////////////////////////////////////////////////////////////////////

        free_matrix(A, N);
        free(nu);
        free(an);
        free(x);
        free(error);

    }

    return N;
}

