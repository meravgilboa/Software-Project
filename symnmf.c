#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

double **alloc_matrix(int rows, int cols);
void free_matrix(double **matrix, int rows);
double euclidean_distance(double *p1, double *p2, int dimension);
double **pow_half_D(double **D, int n);
double **multiply_mat(double **A, double **B, int rows_A, int cols_A, int rows_b, int cols_B);
double frobenius_distance(double **A, double **B, int rows, int cols);
double **transpose(double **Mat, int rows, int cols);
double **calc_A(double **points, int n, int dimension);
double **calc_D(double **A, int n);
double **calc_W(double **D, double **A, int n);
double **update_H(double **H, double **W, int n, int k);
double **converge_H(double **H, double **W, int n, int k);
int *derive_solution(double **H, int n, int k);
void print_result(double** result, int n);
/*Allocate a matrix and returns it*/
double **alloc_matrix(int rows, int cols) {
    double **matrix = malloc(rows*sizeof(double*));
    int i;
    for (i=0; i<rows; i++){
        matrix[i] = malloc(cols*sizeof(double));
    }
    return matrix;
}
/*Free the matrix*/
void free_matrix(double **matrix, int rows) {
    int i;
    for (i=0; i<rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}
/*Calculate the euclidean distance between two points*/
double euclidean_distance(double* p1, double* p2, int dimension) {
    double sum;
    int i;
    sum = 0;
    for (i=0; i<dimension; i++){
        sum += ((p1[i]-p2[i]) * (p1[i]-p2[i]));
    }
    return sum;
}
/*Calculate the D^(-0.5) matrix*/
double **pow_half_D(double** D, int n) { 
    double** E;
    int i, j;
    E = alloc_matrix(n, n);
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            E[i][j] = 0;
        }
        E[i][i] = 1.0 / sqrt(D[i][i]);
    }
    return E;
}
/*Multiply between 2 matrices, and return a new matrix*/
double **multiply_mat(double **A, double **B, int rows_A, int cols_A, int rows_b, int cols_B) {
    double **C;
    int i, j, k;
    C = alloc_matrix(rows_A, cols_B);
    if (cols_A != rows_b){
        printf("An Error Has Occured");
        exit(1);
    }
    for (i=0; i<rows_A; i++){
        for (j=0; j<cols_B; j++){
            C[i][j] = 0;
            for (k=0; k<cols_A; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}
/*Calculate the frobenius distance between two matrices*/
double frobenius_distance(double **A, double **B, int rows, int cols) {
    double sum = 0;
    int i, j;
    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++){
            sum += pow(A[i][j]-B[i][j], 2);
        }
    }
    return sum;
}
/*Transpose a matrix*/
double **transpose(double **Mat, int rows, int cols) {
    double** T;
    int i, j;
    T = alloc_matrix(cols, rows);
    for (i=0; i<cols; i++)
    {
        for (j=0; j<rows; j++)
        {
            T[i][j] = Mat[j][i];
        }
    }
    return T;
}
/*Calculate the A matrix*/
double **calc_A(double **points, int n, int dimension) {
    double** A;
    int i,j;
    A = malloc(n*sizeof(double*));
    for (i=0; i<n; i++){
        A[i] = malloc(n*sizeof(double));
        for (j=0; j<n; j++){
            if (i==j){
                A[i][j] = 0;
            }
            else {
                A[i][j] = exp(-0.5 * euclidean_distance(points[i], points[j], dimension));
            }
        }
    }
    return A;
}
/*Calculate the D matrix*/
double **calc_D(double **A, int n) {
    double** D;
    int i,j;
    D = malloc(n*sizeof(double*));
    for (i=0; i<n; i++){
        D[i] = calloc(n, sizeof(double));
        for (j=0; j<n; j++){
            D[i][i] += A[i][j];
        }
    }
    return D;
}
/*Calculate the W matrix*/
double **calc_W(double **D, double **A, int n) {
    double **C, **E, **W;
    E = pow_half_D(D, n);
    C = multiply_mat(E, A, n, n, n, n);
    W = multiply_mat(C, E, n, n, n, n);
    free_matrix(C, n);
    free_matrix(E, n);
    return W;
}
/*Update the H matrix*/
double **update_H(double **H, double **W, int n, int k) {
    float beta = 0.5;
    double **WH, **Ht, **HHt, **HHtH, **H1;
    int i, j;
    WH = multiply_mat(W, H, n, n, n, k);
    Ht = transpose(H, n, k);
    HHt = multiply_mat(H, Ht, n, k, k, n);
    HHtH = multiply_mat(HHt, H, n, n, n, k);
    H1 = alloc_matrix(n, k);
    for (i=0; i<n; i++){
        for (j=0; j<k; j++){
            H1[i][j] = H[i][j] * ((1-beta) + beta * (WH[i][j] / HHtH[i][j]));
        }
    }
    free_matrix(WH, n);
    free_matrix(Ht, k);
    free_matrix(HHt, n);
    free_matrix(HHtH, n);
    return H1;
}
/*Converge and calculate the H matrix*/
double **converge_H(double **H, double **W, int n, int k) {
    int i, j, l;
    double epsilon, max_iter;
    double** H1;
    epsilon = 0.0001;
    max_iter = 300;
    for (i=0; i<max_iter; i++){
        H1 = update_H(H, W, n, k);
        /*if converge:*/
        if (frobenius_distance(H1, H, n, k) < epsilon){
            return H1;
        }
        /*else, update H and free H1 to call the function again*/
        for (j=0; j<n; j++){
            for (l=0; l<k; l++){
                H[j][l] = H1[j][l];
            }
        }
        free_matrix(H1, n);
    }
    return H;
}
/*Derive the clustering solution from the H matrix*/
int *derive_solution(double **H, int n, int k) {
    int* solution;
    int i, j;
    solution = (int *)malloc(n*sizeof(int));
    for (i=0; i<n; i++){
        solution[i] = 0;
        for (j=0; j<k; j++){
            if (H[i][j] > H[i][solution[i]]){
                solution[i] = j;
            }
        }
    }
    return solution;
}
/*Print the result*/
void print_result(double** result, int n) {
    int i, j;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){ 
            printf("%.4f", result[i][j]);
            if (j<n-1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]){
    double **points;
    double **A, **D, **W;
    FILE *f;
    char *filename, *goal;
    char ch;
    int n, dimension, i, j;
    n = 0;
    dimension = 0;
    /*check CMD arguments*/
    if (argc != 3){
        printf("An Error Has Occured");
        exit(1);
    }
    goal = argv[1];
    filename = argv[2];
    /*Extract the points from the file:*/
    f = fopen(filename, "r");
    if (f == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    while ((ch = fgetc(f)) != EOF){ /*Count the number of points and the dimension*/
        if (ch == '\n'){
            n++;
        }
        else if ((n==0) && (ch==',')){
            dimension++;
        }
    }
    dimension++;
    rewind(f);
    points = malloc(n*sizeof(double*));
    for (i=0; i<n; i++){
        points[i] = malloc(dimension*sizeof(double));
        for (j=0; j<dimension; j++){
            if (fscanf(f, "%lf,", &points[i][j]) != 1)
            {
                printf("An Error Has Occured");
                exit(1);
            }
        }
    }
    fclose(f);
    /*Check what goal equals to using strcmpr and call the right function*/
    A = calc_A(points, n, dimension);
    if (strcmp(goal, "sym") == 0){
        print_result(A, n);
        free_matrix(A, n);
    }
    else {
        D = calc_D(A, n);
        if (strcmp(goal, "ddg") == 0){
            print_result(D, n);
            free_matrix(D, n);
        }
        else if (strcmp(goal, "norm") == 0){
            W = calc_W(D, A, n);
            print_result(W, n);
            free_matrix(D, n);
            free_matrix(W, n);
        }
        else{
            printf("An Error Has Occured");
            exit(1);
        }
    }
    /*free all memory*/
    free_matrix(points, n);
    return 0;
}
