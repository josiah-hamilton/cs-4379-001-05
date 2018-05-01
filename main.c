#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#ifndef NODES
    #define NODES 8
#endif
#define RANDMAX 20.0
#define FMAX 20.0
#define FMIN 0.5

void bound(int*,int*,int,int,int);
void distgen(double *, int, time_t);

int main(int argc, char **argv) {

    int n = NODES;
    int i,j,hii, hij, loi, loj;
    int rank,root=0,size;
    double *x,diff, tmp; 
    double *F, *Ftmp;

    x = (double *) calloc(n, sizeof(int));
    F = (double *) calloc(n, sizeof(double));
    Ftmp = (double *) calloc(n, sizeof(double));

    distgen(x,n,time(NULL));
    for (i=0; i<n; i++) { 
        F[i] = 0.0; 
        Ftmp[i] = 0.0; 
    }
#ifdef DEBUG
    fprintf(stderr,"Generated particles\n");
#endif

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    bound(&loi,&loj,n,rank,size);
    bound(&hii,&hij,n,rank+1,size);
    if (rank != size-1) {
        for (i=loi; i<hii; i++) {
            for (j=loj; j<i; j++) {
                diff=x[i]-x[j];
                tmp = 1.0/(diff*diff*diff);
                Ftmp[i] += tmp;
                Ftmp[j] -= tmp;
            }
        }
        for (j=loj; j<hij; j++) {
            diff=x[i]-x[j];
            tmp = 1.0/(diff*diff*diff);
            Ftmp[i] += tmp;
            Ftmp[j] -= tmp;
        }
    } else { // cover any final indeces that may have slipped through bounds calc
        for (i=loi; i<n; i++) {
            for (j=loj; j<i; j++) {
                diff=x[i]-x[j];
                tmp = 1.0/(diff*diff*diff);
                Ftmp[i] += tmp;
                Ftmp[j] -= tmp;
            }
        }
    }
    
#ifdef DEBUG
    fprintf(stderr,"%d: Calculated forces\n",rank);
#endif


    MPI_Reduce(Ftmp,F,n,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
#ifdef DEBUG
    fprintf(stderr,"%d: Sent to root\n",rank);
#endif

    MPI_Finalize();
    return 0;
}

void bound (int *i, int *j, int n, int rank, int size) {
    int k,computation_bound;
    double bound, numexec=0; // sigma bound, number of executions

    // SUM(k=1,n, k)
    for (k=0; k<n; k++) {
        numexec+=k;
    }

    // example, for 8 particles, rank 1 will perform calculations [7,14).
    computation_bound = rank * (numexec / size);

    // formula for determining the bounds of a sigma series, adapted to partitioning
    bound = 1/2 * ( sqrt(computation_bound*8+1) - 1 );

    // break the bounds into i and j iterator components
    *i = (int)ceil(bound);
    *j = (int)((bound - floor(bound)) * (ceil(bound)));
}

void distgen(double *x, int n, time_t seed) {
    srand(seed);
    int i;
    for (i=0; i<n; i++) {
        x[n] = (double) (FMIN + (rand() / RANDMAX) * (FMAX - FMIN));
    }
}
