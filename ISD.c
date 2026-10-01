#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

#define K 500
#define N 1000
#define T 20

// fonction pour calculé le poids
int weight_calc(int e[K]) {
    int w = 0;
    #pragma omp parallel for reduction(+:w)
    for (int i = 0; i < K; i++)
        if (e[i] == 1) w++;
    return w;
}
// fire une pérmutation aléatoire (pour l'appliquer sur H)
void permutation(int p[N]) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) 
        p[i] = i;

    for (int i = N - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = p[i]; p[i] = p[j]; p[j] = temp;
    }
}
// utilisé dans gauss pour construire le pivot
void swap_rows(int h_copy[K][N], int s_copy[K], int r1, int r2) {
    #pragma omp parallel for schedule(static)
    for (int col = 0; col < N; col++) {
        int temp = h_copy[r1][col];
        h_copy[r1][col] = h_copy[r2][col];
        h_copy[r2][col] = temp;
    }
    int temp_s = s_copy[r1];
    s_copy[r1] = s_copy[r2];
    s_copy[r2] = temp_s;
}
//  pour mettre la matrice sous forme systimatique
bool gauss(int h_copy[K][N], int s_copy[K]) {
    for (int i = 0; i < K; i++) {
        int index = K + i;
        if (h_copy[i][index] == 0) {
            bool found = false;
            for (int j = i + 1; j < K; j++) {
                if (h_copy[j][index] == 1) {
                    swap_rows(h_copy, s_copy, i, j);
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        for (int j = 0; j < K; j++) {
            if (i != j && h_copy[j][index] == 1) {
                for (int col = 0; col < N; col++)
                    h_copy[j][col] ^= h_copy[i][col];
                s_copy[j] ^= s_copy[i];
            }
        }
    }
    return true;
}
// mettre tous les fonctions necéssaire pour ISD
void ISD(int H[K][N], int s[K], int e_found[N]) {
    while (1) {
        int h_prime[K][N];
        int s_copy[K];

        #pragma omp parallel for
        for (int i = 0; i < K; i++) {
            s_copy[i] = s[i];
            for (int j = 0; j < N; j++)
                h_prime[i][j] = H[i][j];
        }

        int p[N];
        permutation(p);

        int h_perm[K][N];
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < K; i++)
            for (int j = 0; j < N; j++)
                h_perm[i][j] = h_prime[i][p[j]];

        if (!gauss(h_perm, s_copy)) continue;

        if (weight_calc(s_copy) <= T) {
            printf("Solution found\n");
            printf("Weight of the found error vector: %d\n", weight_calc(s_copy));

            int e_perm[N];
            #pragma omp parallel for
            for (int i = 0; i < N; i++) e_perm[i] = 0;
            for (int i = 0; i < K; i++)
                e_perm[K + i] = s_copy[i];
            #pragma omp parallel for
            for (int i = 0; i < N; i++) e_found[i] = 0;
            #pragma omp parallel for
            for (int j = 0; j < N; j++)
                e_found[p[j]] = e_perm[j];

            return;
        }
    }
}

int main() {
    srand(time(NULL));

    int H[K][N];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            H[i][j] = rand() % 2;

    int e[N];
    for (int i = 0; i < N; i++) e[i] = 0;

    int count = 0;
    while (count < T) {
        int pos = rand() % N;
        if (e[pos] == 0) { e[pos] = 1; count++; }
    }

    printf("The error vector: [");
    for (int i = 0; i < N; i++) printf(" %d,", e[i]);
    printf("]\n");

    int s[K];
    #pragma omp parallel for
    for (int i = 0; i < K; i++) {
        s[i] = 0;
        for (int j = 0; j < N; j++)
            s[i] ^= (H[i][j] & e[j]);
    }

    printf("Syndrome calculated:\n");
    for (int i = 0; i < K; i++) printf("%d ", s[i]);
    printf("\n");

    printf("Lancement de l'algorithme ISD...\n\n");

    double start = omp_get_wtime();

    int e_found[N];
    ISD(H, s, e_found);

    printf("Error founded [");
    for (int j = 0; j < N; j++) printf("%d ,", e_found[j]);
    printf("]\n");

    bool match = true;
    for (int j = 0; j < N; j++) {
        if (e_found[j] != e[j]) { match = false; break; }
    }
    printf("Match with original error: %s\n", match ? "YES" : "NO");

    double end = omp_get_wtime();
    printf("Time taken: %f seconds\n", end - start);

    return 0;
}