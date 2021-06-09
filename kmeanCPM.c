#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define N 600000
#define G 200

long V[N];
long R[G];
long R2[G];
int A[G];
/* rank del proces    */
int    el_meu_rank;
/* numero de processos        */
int    world_size;

void kmean(int fN, int fK, long fV[], long fR[], int fA[]) {
    int i, j, min, iter = 0;
    long dif, t, dif_global;
    long fS[G];
    //int fD[N];
    // int fD_local[N];
    int *fD = (int *)malloc(sizeof(int) * N);
    int *fD_local = (int *)malloc(sizeof(int) * N);

    do {
        //for (i = 0; i < fN; i++) {
        for(i=(el_meu_rank)*(fN/world_size); i<(fN/world_size)*(el_meu_rank+1); i++) {
            min = 0;
            dif = abs(fV[i] - fR[0]);
            for (j = 1; j < fK; j++)
                if (abs(fV[i] - fR[j]) < dif) {
                    min = j;
                    dif = abs(fV[i] - fR[j]);
                }
            fD_local[i] = min;
        }
        MPI_Allgather(
            &fD_local[(el_meu_rank)*(fN/world_size)],
            fN/world_size,
            MPI_INT,
            fD,
            fN/world_size,
            MPI_INT,
            MPI_COMM_WORLD
        );


        for (i = 0; i < fK; i++)
            fS[i] = fA[i] = 0;

        for (i = 0; i < fN; i++) {
            fS[fD[i]] += fV[i];
            fA[fD[i]]++;
        }

        dif = 0;
        for (i = 0; i < fK; i++) {
            t = fR[i];
            if (fA[i]) fR[i] = fS[i] / fA[i];
            dif += abs(t - fR[i]);
        }
        iter++;
    } while (dif);

    if (el_meu_rank == 0) printf("iter %d\n", iter);
}

void qs(int ii, int fi, long fV[], int fA[]) {
    int i, f, j;
    long pi, pa, vtmp, vta, vfi, vfa;

    pi = fV[ii];
    pa = fA[ii];
    i = ii + 1;
    f = fi;
    vtmp = fV[i];
    vta = fA[i];

    while (i <= f) {
        if (vtmp < pi) {
            fV[i - 1] = vtmp;
            fA[i - 1] = vta;
            i++;
            vtmp = fV[i];
            vta = fA[i];
        } else {
            vfi = fV[f];
            vfa = fA[f];
            fV[f] = vtmp;
            fA[f] = vta;
            f--;
            vtmp = vfi;
            vta = vfa;
        }
    }
    fV[i - 1] = pi;
    fA[i - 1] = pa;

    if (ii < f) qs(ii, f, fV, fA);
    if (i < fi) qs(i, fi, fV, fA);
}

int main() {
    int i;
    /* Inicialitzar MPI */
    MPI_Init(NULL, NULL);
    /* Obtenir el rank del proces  */
    MPI_Comm_rank(MPI_COMM_WORLD, &el_meu_rank);
    /* Obtenir el numero total de processos */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    for (i=0;i<N;i++) V[i] = (rand()%rand())/N;

    // primers candidats
    for (i = 0; i < G; i++) R[i] = V[i];

    // calcular els G mes representatius
    kmean(N, G, V, R, A);

    qs(0, G - 1, R, A);

    if (el_meu_rank == 0)
    {
        for (i = 0; i < G; i++)
            printf("R[%d] : %ld te %d agrupats\n", i, R[i], A[i]);
    }

    
    MPI_Finalize();
    return (0);
}