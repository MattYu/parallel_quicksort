#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

void srandom(unsigned seed);

/*
Reference Paper:
http://www.winlab.rutgers.edu/~pkataria/pdc.pdf
*/


int partition(int* A, int start, int end, int pivot);
void localQuickSort(int *a, int start, int end);


int main(int argc, char*argv[])
{

    if (argc != 5)
    {
        printf("Enter value size argument - Number of Dimension - Printing mode (0 or 1) - Gather (0 or 1)", argv[0]);
        return 0;
    }

    int t_size = atoi(argv[1]);
    int temp2 = t_size;
    int *a = malloc(t_size * sizeof(int));
    
    MPI_Comm childcomm;
    int id;
    int D = atoi(argv[2]);
    int size = t_size/pow(2,D);
    int rank;
    int pivot;
    int printMode = atoi(argv[3]);
    int gather = atoi(argv[4]);


    // Pre-determine common pivots used at each dimension
    int *pivotProcessor = malloc((D+1) * sizeof(int));
    for (int i=1; i<=D;i++){
        pivotProcessor[i] = rand()%(int)pow(2,D);
    }
    clock_t start_time, end_time;
    start_time = clock();
    MPI_Init(&argc, &argv);
    // Processor 0 generate random values bounded within 100 and broadcast it to the rest of the processors 
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0){
        printf("\nInitializing random array of size %d", t_size);
        for (int i = 0; i< t_size; i++){
        a[i] = rand()%100;
        }
    }
    MPI_Bcast(a, t_size, MPI_INT, 0, MPI_COMM_WORLD);
    if (id == 0){
        printf("\nRandom array broadcasted");
    }

    //
    if (id == 0 && printMode == 1){
        printf("\nInitial unsorted array\n");
        for (int i = 0; i <t_size; i++){
            printf("\t%d", a[i]);
        }
    }

    int *container = malloc(size*sizeof(int));
    int j =0;
    for (int i= size*id; i < size*id + size; i++){
        container[j++] = a[i];
    }

    // Beginning of hypercube traversal
    for (int i = D; i > 0;i--){
        //printf("\n%d Processor %d contains %d elements", i, id, size);
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0){
            printf("\nDimension %d", i);
        }
        // One processor selects a pivot. If the processor contains value, it will pick one of its random values. Otherwise, the pivot will be a random int bounded within 100
        if (id == pivotProcessor[i]){
            if (size> 0){
                pivot = container[rand()%size];
            }
            else {
                pivot = rand()%100;
            }
        }

        MPI_Bcast(&pivot, 1, MPI_INT, pivotProcessor[i], MPI_COMM_WORLD);

        /*if (id == 0){
            printf("\nPivot value for dimension %d = %d", i, pivot);
        }*/
        // Each processor is paired with another processor who's ith ID bit is complementary to its own
        int pairedId = id^(1 << (i-1));
        int smallSize = 0;
        int largeSize = 0;
        int *smallArray;
        int *largeArray;
        //printf("\n!!%d %d is paired with %d", i, id, pairedId);
        // Each processor locally sort around the pivot
        for (int i = 0; i < size; i++){
            if (container[i] <= pivot) smallSize++;
            else largeSize++;
        }

        smallArray = malloc(smallSize*sizeof(int));
        largeArray = malloc(largeSize*sizeof(int));
        int j = 0;
        int k = 0;
        for (int i=0; i < size; i++){
            if (*(container +i) <= pivot) smallArray[j++] = *(container + i);
            else largeArray[k++] = *(container + i);
        }
        //Communication between pair processors happens below and each processor does a compare-split exchange with its paired processor
        //Note to avoid a deadlock, one processor send data first than receives data, the other processor receives data first and then send data
        free(container);
        if (!(id&(1 << (i - 1)))){
            MPI_Send(&largeSize, 1, MPI_INT, pairedId, 0, MPI_COMM_WORLD);
            //MPI_Bcast(&largeSize, 1, MPI_INT, tempId, childcomm);
            if (largeSize > 0){
                MPI_Send(largeArray, largeSize, MPI_INT, pairedId, 1, MPI_COMM_WORLD);
                //MPI_Bcast(largeArray, largeSize, MPI_INT, tempId, childcomm);
            }
            int oldSmall = smallSize;
            MPI_Recv(&smallSize, 1, MPI_INT, pairedId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Bcast(&smallSize, 1, MPI_INT, pairedTempId, childcomm);
            container = malloc((oldSmall + smallSize)*sizeof(int));
            for (int i = 0; i < oldSmall; i++){
                container[i] = smallArray[i];
            } 
            if (smallSize > 0){
                free(smallArray);
                smallArray = malloc(smallSize* sizeof(int));
                MPI_Recv(smallArray, smallSize, MPI_INT, pairedId, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Bcast(smallArray, smallSize, MPI_INT, pairedTempId, childcomm);
                
                for (int i = 0; i < smallSize; i++){
                    container[oldSmall++] = *(smallArray + i);
                } 
            }
            size = oldSmall;
            free(largeArray);
            free(smallArray);

            printf("\nFor dimension %d, %d is the complement of %d\n", i, pairedId, id);
            //printf("\nTesting- %d : %d : %d received smallsize=%d - large=%d - size=%d - pivot=%d\n", i, id, pairedId, smallSize, largeSize, size, pivot);
        }
        else{
            int oldLarge = largeSize;
            MPI_Recv(&largeSize, 1, MPI_INT, pairedId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Bcast(&largeSize, 1, MPI_INT, pairedTempId, childcomm);
            container = malloc((oldLarge + largeSize)*sizeof(int));
            for (int i = 0; i < oldLarge; i++){
                container[i] = largeArray[i];
            } 
            
            if (largeSize > 0){
                free(largeArray);
                largeArray = malloc(largeSize * sizeof(int));
                MPI_Recv(largeArray, largeSize, MPI_INT, pairedId, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Bcast(largeArray, largeSize, MPI_INT, pairedTempId, childcomm);
                for (int i = 0; i < largeSize; i++){
                    container[oldLarge++] = *(largeArray + i);
                } 
            }
            size = oldLarge;
            //MPI_Bcast(&smallSize, 1, MPI_INT, tempId, childcomm);
            MPI_Send(&smallSize, 1, MPI_INT, pairedId, 0, MPI_COMM_WORLD);
            if (smallSize > 0){
                MPI_Send(smallArray, smallSize, MPI_INT, pairedId, 1, MPI_COMM_WORLD);
                //MPI_Bcast(smallArray, smallSize, MPI_INT, tempId, childcomm);
            }
            free(smallArray);
            free(largeArray);
            printf("\nFor dimension %d, %d is the complement of %d\n", i, pairedId, id);
            //printf("\nTesting- %d : %d : %d received smallsize=%d - large=%d - size=%d - pivot=%d\n", i, id, pairedId, smallSize, largeSize, size, pivot);
        }
        printf("\n**%d Processor %d exiting dimension with %d elements at time %d", i, id, size, clock()-start_time);
        if (id == 0){
            printf("\nExiting Dimension %d", i);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Every processor does a local sort at the end
    localQuickSort(container, 0, size-1);

    if (id == 0){
            //printf("\nInitial Array results:\n", id, pivot);
            for (int i = 0; i <t_size; i++){
                //printf("\t%d", a[i]);
            }
        }

    // Since the data within each processor may be unbalanced, we collect the sorted data to processor 0 using gatherv. 
    // The implementation of gatherv is adapted from https://stackoverflow.com/questions/31890523/how-to-use-mpi-gatherv-for-collecting-strings-of-diiferent-length-from-different
    const int root = 0;
    
    int nProcessor = pow(2,D);
    int *displs = NULL;
    int totlen = 0;
    //MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    int *recvcounts = NULL;
    int *ans = NULL;
    
    if (id == 0){
        recvcounts = malloc(nProcessor * sizeof(int));
        recvcounts[0] = size;
    }
    if (gather == 1){
        MPI_Gather(&size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        if (id == root) {
            /*for (int i = 0; i < nProcessor; i++){
                printf("\nCount for processor %d : %d", i, recvcounts[i]);
            }*/

            displs = malloc( nProcessor * sizeof(int) );

            displs[0] = 0;

            for (int i=0; i<nProcessor; i++) { 
            displs[i] = totlen;
            totlen += recvcounts[i]; 
            }

            /*for (int i = 0; i < nProcessor; i++){
                printf("\nDisplacement for processor %d : %d", i, displs[i]);
            }*/

            ans = malloc(totlen*sizeof(int));
        }
        //if (size != 0){
        MPI_Gatherv(container, size, MPI_INT,
                    ans, recvcounts, displs, MPI_INT,
                    root, MPI_COMM_WORLD);
    }
    

	if (id == 0 && printMode == 1 && gather == 1){
        printf("\nFinal Array results %d:\n", totlen);
        for (int i = 0; i <t_size; i++){
                printf("\t%d", *(ans+i));
        }
    }
    

    if (id == 0){
        end_time = clock();
        printf("Time taken is %f seconds \n", ((double)(end_time - start_time))/CLOCKS_PER_SEC);
    }
    
    MPI_Finalize();
    return 0;
}



int partition(int *A, int start, int end, int pivot){
    // Partionning used by local quick sort
    int i = start;
    int j = end;

    //User temp buffer to perform a rotational swap
    int temp = A[start];

    while (i < j)
    {
        //Element of Programming Interview "Dutch Flag" and "Quick Select" implementation
        while (i < j && A[j] >= pivot) j--;
        A[i] = A[j];
        while (i < j && A[i] <= pivot) i++;
        A[j] = A[i];
    }

    A[i] = temp;
    //if (temp <= pivot) i++;

    return i;
}

void localQuickSort(int *a, int start, int end){
    if (start >= end) return;

    int pivot = a[start];

    int part = partition(a, start, end, pivot);

    localQuickSort(a, start, part-1);
    localQuickSort(a, part+1, end);
}
