#include <stdio.h>
#include <limits.h>
#include <omp.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#ifndef V
#define V 5000
#endif

int num;

int minKeyParallel(int key[], int visited[])
{
    int min = INT_MAX, index, i;
#pragma omp parallel
    {
        num = omp_get_num_threads();
        int index_local = index;
        int min_local = min;
#pragma omp for nowait
        for (i = 0; i < V; i++)
        {
            if (visited[i] == 0 && key[i] < min_local)
            {
                min_local = key[i];
                index_local = i;
            }
        }
#pragma omp critical
        {
            if (min_local < min)
            {
                min = min_local;
                index = index_local;
            }
        }
    }
    return index;
}

int minKeySequential(int key[], int visited[])
{
    int min = INT_MAX, index, i;
    for (i = 0; i < V; i++)
    {
        if (visited[i] == 0 && key[i] < min)
        {
            min = key[i];
            index = i;
        }
    }
    return index;
}

void printMST(int from[], int n, int** graph)
{
    int i;
    printf("Edge   Weight\n");
    for (i = 1; i < V; i++)
        printf("%d - %d    %d \n", from[i], i, graph[i][from[i]]);
}

void primsMST_Parallel(int** graph)
{
    int from[V];
    int key[V], num_threads;
    int visited[V];
    int i, count;
    for (i = 0; i < V; i++)
        key[i] = INT_MAX, visited[i] = 0;

    key[0] = 0;
    from[0] = -1;

    for (count = 0; count < V - 1; count++)
    {
        int u = minKeyParallel(key, visited);
        visited[u] = 1;

        int v;
#pragma omp parallel for schedule(static)
        for (v = 0; v < V; v++)
        {
            if (graph[u][v] && visited[v] == 0 && graph[u][v] < key[v])
                from[v] = u, key[v] = graph[u][v];
        }
    }
}

void primsMST_Sequential(int** graph)
{
    int from[V];
    int key[V], num_threads;
    int visited[V];
    int i, count;
    for (i = 0; i < V; i++)
        key[i] = INT_MAX, visited[i] = 0;

    key[0] = 0;
    from[0] = -1;

    for (count = 0; count < V - 1; count++)
    {
        int u = minKeySequential(key, visited);
        visited[u] = 1;

        int v;
        for (v = 0; v < V; v++)
        {
            if (graph[u][v] && visited[v] == 0 && graph[u][v] < key[v])
                from[v] = u, key[v] = graph[u][v];
        }
    }
}

int main()
{
    int** graph = (int**)malloc(V * sizeof(int*));
    for (int x = 0; x < V; x++)
        graph[x] = (int*)malloc(V * sizeof(int));
    int i, j;

    srand(time(NULL));
    for (i = 0; i < V; i++)
        for (j = 0; j < V; j++)
            graph[i][j] = rand() % 10;

    for (i = 0; i < V; i++)
    {
        graph[i][i] = 0;
    }

    for (i = 0; i < V; i++)
        for (j = 0; j < V; j++)
            graph[j][i] = graph[i][j];

    double start, end;

    start = omp_get_wtime();
    primsMST_Sequential(graph);
    end = omp_get_wtime();
    printf("Time for seq = %f\n", end - start);

    start = omp_get_wtime();
    primsMST_Parallel(graph);
    end = omp_get_wtime();
    printf("Time for par = %f\nThreads = %d\n", end - start, num);

    return 0;
}