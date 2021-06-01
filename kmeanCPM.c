#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#define N 600000
#define G 200

long V[N];
long R[G];
int A[G];

void kmean(int fN, int fK, long fV[], long fR[], int fA[]) {
    int i, j, min, iter = 0;
    long dif, t;
    long fS[G];
    int *fD = malloc(sizeof(int) * N);
    int *temp_fD = malloc(sizeof(int) * N);
    int world_size, el_meu_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &el_meu_rank);


    do {
        MPI_Bcast(fR,G,MPI_LONG,0,MPI_COMM_WORLD);
        MPI_Bcast(fV,N,MPI_LONG,0,MPI_COMM_WORLD);
        // MPI_Send(fR, G, MPI_LONG, 1, 0, MPI_COMM_WORLD);
        // MPI_Send(fV, N, MPI_LONG, 1, 0, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD); 
        for (int p = 1; p < world_size; p++)
        {
            MPI_Recv(temp_fD, N, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //memcpy(fD + (el_meu_rank-1)*(N/(world_size-1)) * sizeof(int), 
            //    temp_fD + (el_meu_rank-1)*(N/(world_size-1)) * sizeof(int),
            //    sizeof(int)*(N/(world_size-1)));
        }
        fD = temp_fD;
        
        //printf("0 - %i\n", fD[iter]);
        
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
        MPI_Send(&dif, 1, MPI_LONG, 1, 0, MPI_COMM_WORLD);
    } while (dif);

    printf("iter %d\n", iter);
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
    /* rank del proces    */
    int    el_meu_rank;
    /* numero de processos        */
    int    world_size;
    /* rank de l'emissor */
    int    font;
    /* rank del receptor */
    int    desti;
    /* etiqueta dels missatges */
    int    etiq = 0;
    /* espai per al missatge      */
    char missatge[100];    
    /* estat de la recepcio       */
    MPI_Status  estat;

    /* Inicialitzar MPI */
   MPI_Init(NULL, NULL);
    /* Obtenir el rank del proces  */
   MPI_Comm_rank(MPI_COMM_WORLD, &el_meu_rank);
    /* Obtenir el numero total de processos */
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    if (el_meu_rank == 0) {
        int i;

        for (i = 0; i < N; i++) V[i] = (rand() % rand()) / N;

        // primers candidats
        for (i = 0; i < G; i++) R[i] = V[i];

        sprintf(missatge, "Salutacions desde el proces %d ",el_meu_rank);
        //MPI_Bcast(missatge,101,MPI_CHAR,0,MPI_COMM_WORLD);
        // calcular els G mes representatius

        kmean(N, G, V, R, A);

        qs(0, G - 1, R, A);

        for (i = 0; i < G; i++)
            printf("R[%d] : %ld te %d agrupats\n", i, R[i], A[i]);
    }
    else {
        int i=0, min, j, iter=0;
        long dif,  temp ;
        long *fR = malloc(sizeof(long) * G) ;
        long *fV = malloc(sizeof(long) * N);
        int *fD = malloc(sizeof(int) * N);
        do
        {
            MPI_Bcast(fR,G,MPI_LONG,0,MPI_COMM_WORLD);
            MPI_Bcast(fV,N,MPI_LONG,0,MPI_COMM_WORLD);
            // MPI_Recv(fR,G,MPI_LONG,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // MPI_Recv(fV,N,MPI_LONG,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Barrier(MPI_COMM_WORLD);
            
            for (i = ((el_meu_rank-1)*(N/(world_size-1))); i < (N/(world_size-1)); i++) {
                min = 0;
                temp = abs(fV[i] - fR[0]);
                for (j = 1; j < G; j++)
                    if (abs(fV[i] - fR[j]) < temp) {
                        min = j;
                        temp = abs(fV[i] - fR[j]);
                    }
                fD[i] = min;
            }
            //printf("1 - %i\n", world_size);
            MPI_Send(fD, N, MPI_INT, 0, 0, MPI_COMM_WORLD);

            //printf("1 - %d\n", fD[iter]);
            iter++;
            // Para seguir bucles
            MPI_Recv(&dif,1,MPI_LONG,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } while (dif);


        //kmean(N, G, V, R, A);
        //MPI_Bcast(missatge,101,MPI_CHAR,0,MPI_COMM_WORLD);
        //printf("Soy: %i %s\n",el_meu_rank, missatge);
    }


    MPI_Finalize();
    return (0);
}