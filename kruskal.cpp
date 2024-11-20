#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>
#include <omp.h>

class DSU {
    int num;
    std::vector<int> parent;
    std::vector<int> size;

    public:
    DSU(const int& n) {
        num = n + 5;
        parent.resize(num, 0);
        size.resize(num, 0);
    }

    void make_set(int v) {
        parent[v] = v;
        size[v] = 1;
    }

    int find_set(int x) {
        if (x == parent[x]) return x;
        return parent[x] = find_set(parent[x]);
    }

    int set_size(int x) {
        return size[find_set(x)];
    }

    void union_set(int a, int b) {
        a = find_set(a);
        b = find_set(b);

        if (a != b) {
            if (size[a] < size[b]) std::swap(a, b);
            parent[b] = a;
            size[a] += size[b];
        }
    }
};

#ifndef N
#define N 100
#endif

std::vector<int> temp[N * N];
int graph[N][N] = {0};

void merge_sequential(std::vector<int> a[], int n, std::vector<int> b[], int m, std::vector<int> c[]) {
    for (int i = 0; i < n; i++) {
        temp[
            i + (std::upper_bound(b, b + m, a[i]) - b)
        ] = a[i];
    }
    for (int i = 0; i < m; i++) {
        temp[
            i + (std::upper_bound(a, a + n, b[i]) - a)
        ] = b[i];
    }
    for (int i = 0; i < n + m; i++)
        c[i] = temp[i];

}

void merge_parallel(std::vector<int> a[], int n, std::vector<int> b[], int m, std::vector<int> c[]) {
    if (n + m < 10) {
        merge_sequential(a, n, b, m, c);
        return;
    }
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < n; i++) {
            temp[
                i + (std::upper_bound(b, b + m, a[i]) - b)
            ] = a[i];
        }
#pragma omp for
        for (int i = 0; i < m; i++) {
            temp[
                i + (std::upper_bound(a, a + n, b[i]) - a)
            ] = b[i];
        }
#pragma omp for
        for (int i = 0; i < n + m; i++)
            c[i] = temp[i];
    }
}

void merge_sort_sequential(std::vector<int> a[], int n) {
    if (n == 1)
        return;
    int half = n / 2;
    merge_sort_sequential(a, half);
    merge_sort_sequential(a + half, n - half);
    merge_sequential(a, half, a + half, n - half, a);
}

void merge_sort_parallel(std::vector<int> a[], int n) {
    if (n == 1)
        return;
    int half = n / 2;
    merge_sort_parallel(a, half);
    merge_sort_parallel(a + half, n - half);
    merge_parallel(a, half, a + half, n - half, a);
}

void kruskal_sequential(std::vector<int> arr[], int n) {
    merge_sort_sequential(arr, n);
    DSU ds(n);
    for (int i = 0; i < n; i++)
        ds.make_set(i);

    int u, v;

    for (int i = 0; i < n; i++) {
        u = arr[i][1];
        v = arr[i][2];

        if (ds.find_set(u) != ds.find_set(v)) {
            ds.union_set(u, v);
        }
    }
}

void kruskal_parallel(std::vector<int> arr[], int n) {
    merge_sort_parallel(arr, n);
    DSU ds(n);
    for (int i = 0; i < n; i++)
        ds.make_set(i);

    int u, v;

    for (int i = 0; i < n; i++) {
        u = arr[i][1];
        v = arr[i][2];

        if (ds.find_set(u) != ds.find_set(v)) {
            ds.union_set(u, v);
        }
    }
}

int main() {

    int ed = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            graph[i][j] = rand() % 10;
            ed += (graph[i][j] > 0);
        }

    std::vector<int> arr1[ed];
    std::vector<int> arr2[ed];

    int cnt = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            if (graph[i][j]) {
                arr1[cnt++] = {graph[i][j], i, j};
            }
        }

    for (int i = 0; i < ed; i++)
        arr2[i] = arr1[i];

    double start, end;

    start = omp_get_wtime();
    kruskal_sequential(arr1, ed);
    end = omp_get_wtime();
    std::cout << "Sequential Time: " << end - start << std::endl;

    start = omp_get_wtime();
    kruskal_parallel(arr2, ed);
    end = omp_get_wtime();
    std::cout << "Parallel Time: " << end - start << std::endl;

}