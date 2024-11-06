void free_matrix(double **matrix, int rows);
double **calc_A(double **points, int n, int dimension);
double **calc_D(double **A, int n);
double **calc_W(double **D, double **A, int n);
double **update_H(double **H, double **W, int n, int k);
double **converge_H(double **H, double **W, int n, int k);
int *derive_solution(double **H, int n, int k);