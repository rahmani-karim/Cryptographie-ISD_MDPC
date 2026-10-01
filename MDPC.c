#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define W 35
#define N 4813
#define T 26
#define E 78
// Une fonction qui génere des vecteurs aléatoire avec un poids voulu  
void generation_vecteurs(int *v, int weight) {
    int compt = 0;
    while(compt< weight){
        int pos=rand()%N;
        if (v[pos]==0){
            v[pos]=1;
            compt++;
        }
    }
}
 //Fonction de calculs de poids 
int poids(int v[N]) {
    int w = 0;
    for (int i = 0; i < N; i++) {
        if (v[i] == 1) w++;
    }
    return w;
}

void multiplication_circulaire(int *res,int A[N],int B[N]) {
    for(int i=0;i<N;++i){
        if (A[i]==1){
            for(int j=0; j<N;++j){
                if(B[j]==1){
                    int index=(i+j)%N;
                    res[index]^=1;
                }
            }
        }
    }
}

int get_degre(int poly[], int taille) {
    for (int i = taille - 1; i >= 0; i--) {
        if (poly[i] == 1) return i;
    }
    return -1; 
}

bool inverse_poly(int y[N], int y_inv[N]) {
    int r0[N + 1] = {0};
    int r1[N + 1] = {0};
    // u0 et u1 restent strictement de taille N !
    int u0[N] = {0}; 
    int u1[N] = {0};

    r0[0] = 1; r0[N] = 1; // r0 = X^N + 1
    for (int i = 0; i < N; i++) r1[i] = y[i];
    u1[0] = 1;

    while (1) {
        int d1 = get_degre(r1, N + 1);
        if (d1 == -1) break;

        int d0 = get_degre(r0, N + 1);

        if (d0 < d1) {
            for (int i = 0; i <= N; i++) {
                int temp = r0[i]; r0[i] = r1[i]; r1[i] = temp;
            }
            for (int i = 0; i < N; i++) {
                int temp = u0[i]; u0[i] = u1[i]; u1[i] = temp;
            }
            int temp_d = d0; d0 = d1; d1 = temp_d;
        }

        int shift = d0 - d1;
        
        // Division sur le reste (classique)
        for (int i = 0; i <= d1; i++) {
            if (r1[i] == 1) r0[i + shift] ^= 1;
        }

        // Division sur le coefficient de Bezout (Modulo X^N - 1)
        for (int i = 0; i < N; i++) {
            if (u1[i] == 1) {
                // LA MAGIE EST ICI : Le modulo empêche la perte du bit X^N
                int index = (i + shift) % N; 
                u0[index] ^= 1;
            }
        }
    }

    if (get_degre(r0, N + 1) != 0 || r0[0] != 1) return false;

    for (int i = 0; i < N; i++) y_inv[i] = u0[i];
    return true;
}
// Fonction pour la géneration des clés 
void generation_des_cles(int x[N], int y[N], int H[N]) {
    generation_vecteurs(x, W);
    
    bool inversible=false;
    int y_inv[N]={0};
    //On doit vérifier que y est inversible car on a besoin de y^(-1) 
    while (inversible == false) {
        // On remet y à zéro avant chaque tentative (très important !)
        for (int i = 0; i < N; i++) y[i] = 0;
        
        generation_vecteurs(y,W);
        inversible = inverse_poly(y, y_inv); 
    }
    //On calcule notre H= x*y^(-1)
    multiplication_circulaire(H, x, y_inv);
}


 


void chiffrer(int c0[N], int c1[N], int m[N], int g[N]) {
    int e0[N] = {0};
    int e1[N] = {0};
    // On distribue les erreurs de façons aléatoire sur nos deux vecteurs e_0 et e_1
    int p_e0 =  1+ rand()%(E-1);
    int p_e1 = E-p_e0;
    generation_vecteurs(e0, p_e0);
    generation_vecteurs(e1,p_e1);

    for (int i = 0; i < N; i++) {
        c0[i] = m[i] ^ e0[i];
    }

    int mg[N] = {0};
    multiplication_circulaire(mg, m, g);
    for (int i = 0; i < N; i++) {
        c1[i] = mg[i] ^ e1[i];
    }
}
// calcul du syndrome
void calculer_syndrome(int s[N], int c0[N], int c1[N], int x[N], int y[N]) {
    int s0[N] = {0};
    int s1[N] = {0};
    multiplication_circulaire(s0, c0, x);
    multiplication_circulaire(s1, c1, y);
    for (int i = 0; i < N; i++){ 
        s[i] = s0[i] ^ s1[i];
        }
}
//Déchiffrement en utilisant le bitflipping
void dechiffrer(int m_retrouve[N], int c0[N], int c1[N], int x[N], int y[N]) {
    int s[N] = {0};
    calculer_syndrome(s, c0, c1, x, y);

    int iter = 0;
    int max_iter = 100;

    // Algorithme BitFlip
    while (poids(s) > 0 && iter < max_iter) {
        int scores_c0[N] = {0};
        int scores_c1[N] = {0};

        // Comptage des bits qui impliquent des erreurs 
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (s[j] == 1) {
                    int index = (j - i + N) % N;
                    if (x[index] == 1) scores_c0[i]++;
                    if (y[index] == 1) scores_c1[i]++;
                }
            }
        }



        // Le Flip avec le seuil T
        for (int i = 0; i < N; i++) {
            if (scores_c0[i] >= T) c0[i] ^= 1;
            if (scores_c1[i] >= T) c1[i] ^= 1;
        }

        calculer_syndrome(s, c0, c1, x, y);
        iter++;
    }

    if (poids(s) == 0) {
        printf("BitFlipiing réussie en %d iterations \n", iter);
    } else {
        printf("BitFlipping échoué \n");
    }

    for (int i = 0; i < N; i++) 
        m_retrouve[i] = c0[i];
}

int main() {

    int x[N] = {0}, y[N] = {0}, g[N] = {0};
    int m[N] = {0};
    int c0[N] = {0}, c1[N] = {0};
    int msg_dechiff[N] = {0};
    srand(time(NULL));

     for (int i=0; i<1000; i++) {
         int pos = rand() % N;
         if (m[pos] == 0)
             m[pos] = 1;
     }

    printf(" Génération des clés:\n");
    clock_t start = clock();
    generation_des_cles(x, y, g);
    clock_t end = clock();
    double temps_keygen = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Temps pris pour la Génération des clés est %f: \n", temps_keygen);
    printf("\n");
    printf("Chiffrement:\n");
    start = clock();
    chiffrer(c0, c1, m, g);
    end = clock();
    double temps_encrypt = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Temps pris pour le chiffrement est %f: \n", temps_encrypt);
    printf("\n");
    printf("Déchiffrement:\n");
    start = clock();
    dechiffrer(msg_dechiff, c0, c1, x, y);
    end = clock();
    double temps_decrypt = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Temps pris pour le déchiffrement est %f:\n ", temps_decrypt);

    bool succes = true;
    for(int i = 0; i < N; i++) {
        if (m[i] != msg_dechiff[i]) {
            succes = false;
            break;
        }
    }
    printf("\n");
    if (succes) printf("Le message est déchiffré avec succées \n");
    else printf("Echec de déchiffrement \n");
    
        
    return 0;
}