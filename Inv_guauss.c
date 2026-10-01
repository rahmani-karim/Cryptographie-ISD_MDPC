#include<stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define N 9

void affichage_matrice(int M[N][N]){
    for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                printf("%d ", M[i][j]);
                }
            printf("\n");
            }
    }
void generation_matrice(int M[N][N]){
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            M[i][j] = rand() % 2;
        }

bool inverse_gauss(int M[N][N], int inv_M[N][N]) {
    
    // Construire la matrice augmentée  
    int M2[N][2*N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            M2[i][j] = M[i][j];
            M2[i][N + j] = (i == j) ? 1 : 0;
        }
    }

    // Pivot de Gauss-Jordan
    for (int i = 0; i < N; i++) {

        // Chercher le pivot
        if (M2[i][i] == 0) {
            bool found = false;
            for (int j = i + 1; j < N; j++) {
                if (M2[j][i] == 1) {
                    // Échanger les lignes i et j
                    for (int col = 0; col < 2*N; col++) {
                        int temp = M2[i][col];
                        M2[i][col] = M2[j][col];
                        M2[j][col] = temp;
                    }
                    found = true;
                    break;
                }
            }
            if (!found) 
                return false; // Matrice non inversible
        }

        // Éliminer sur toutes les autres lignes
        for (int j = 0; j < N; j++) {
            if (j != i && M2[j][i] == 1) {
                for (int col = 0; col < 2*N; col++) {
                    M2[j][col] ^= M2[i][col];
                }
            }
        }
    }

    // Extraire la partie droite = M⁻¹
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inv_M[i][j] = M2[i][N + j];

    return true;
}

int main() {
    srand(time(NULL));

    int M[N][N];
    int inv_M[N][N];
    generation_matrice(M);
    printf("On va inverser la matrice suivante \n");
    affichage_matrice( M);
    if (!inverse_gauss(M, inv_M)) {
        printf("Matrice non inversible\n");
        return 1;
    }
    if (inverse_gauss(M,inv_M)){
        printf(" Voici l'inverse de la matrice M \n");
         affichage_matrice( inv_M);
    }

    return 0;
}