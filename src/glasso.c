#include <stdlib.h>
#include <stdio.h>
#include "glasso.h"
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <float.h>

#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif

int     izero       =  0;
int     ione        =  1;
double  dzero       =  0.0;
double  done        =  1.0;
double  dminusone   = -1.0;

// for debugging //

bool isnan_mat(const double* mat, int n){
    bool tmp = false;
    for (int i = 0; i < n; i++){
        if(isnan(mat[i])) { tmp = true; break; }
    }
    return tmp;
}

void print_bmatrix(const bool* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%d, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_dmatrix(const double* matrix,int m,int n){
    int i,j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%5.4f, ", matrix[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}
// end for debugging //

void dmat_elemprod(const int n, const double *x, const double *y, double *z)
{
    if (y != z) {
        F77_CALL(dcopy)(&n, y, &ione, z, &ione);
    }
    F77_CALL(dtbmv)("U", "N", "N", &n, &izero, x, &ione, z, &ione);
}


void dmat_vcopy(const int n, const double *src, double *dst)
{
    F77_CALL(dcopy)(&n, src, &ione, dst, &ione);
}

void dmat_vset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w)

{
#if 1
    if (w != x && w != y)
    {
        dmat_vset(n, 0, w);
        F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
        F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else if (w == x && w == y)
    {
        double tmp;
        tmp = alpha+beta;
        F77_CALL(dscal)(&n, &tmp, w, &ione);
    }
    else if (w == x /*&& w != y */)
    {
        if (alpha != 1.0) F77_CALL(dscal)(&n, &alpha, w, &ione);
        if (beta  != 0.0) F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else /* if (w == y && w != x ) */
    {
        if (beta  != 1.0) F77_CALL(dscal)(&n, &beta , w, &ione);
        if (alpha != 0.0) F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
    }
#else
    int i;
    
    if (beta == 0.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i];
        }
    }
    else if (beta == 1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] + y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] + y[i];
        }
    }
    else if (beta == -1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] - y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] - y[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            w[i] = alpha*x[i] + beta*y[i];
    }
#endif
}

double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}

double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}

void dmat_C_ATB(int m, int n1, int n2,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,&done,B,&n2,A,&n1,&dzero,C,&n2);
}

// not tested //
/* A is a m1 by n matrix; B is a m2 by n matrix. */
void dmat_C_ABT(int n, int m1, int m2, double* beta,
                const double* A, const double* B, double* C)
{
    F77_CALL(dgemm)("T","N",&m2,&m1,&n,&done,B,&n,A,&n,beta,C,&m2);
}

void dmat_B_AAT(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","T",&m,&n,&done,A,&n,&dzero,B,&m);
}

// not tested //
/* A n by m, D n by n */
void dmat_B_ATDA(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    dmat_vset(m*n,.0,Tmp);
    for (i = 0; i < n; i++) {
        // dmat_elemprod(n, A+i*n, D, Tmp+i*n);
        dmat_waxpby(m,D[i], A+i*m,.0,Tmp+i*m,Tmp+i*m);
    }
    // if (isnan_mat(A,n*m)) { printf("A is nan! \n"); }
    // if (isnan_mat(Tmp,n*m)) { printf("Tmp is nan! \n"); }
    dmat_C_ATB(n, m, m, A, Tmp, B);
    // if (isnan_mat(B,m*m)) { printf("B is nan! \n"); }
    free(Tmp);
}


// not tested //
/* A is a m by n matrix, D is a n by n matrix, B is a m by m matrix. */
void dmat_B_ADAT(int m, int n, double *A, double* D, double *B)
{
    int i;
    double* Tmp = (double*) malloc(sizeof(double)*m*n);
    // scale cols of A by D and stores the scaled matrix in Tmp //
    for (i = 0; i < m; i++) {
        dmat_elemprod(n, A+i*n, D, Tmp+i*n);
        // dmat_waxpby(m,D[i], A+i*m,.0, Tmp,Tmp+i*m);
    }
    dmat_C_ABT(n, m, m, &dzero, A, Tmp, B);
    free(Tmp);
}


void proj_L1_gen(double *x, const double *c, const double *z, const int* nn)
{
    int n = nn[0];
    double f_0 = -z[0]; 
    int converge = 1; 
    for (int i = 0; i < n; ++i)
    {
        f_0 = f_0 + c[i] * fabs(x[i]); 
    }
    if (f_0 > .0)
    {
        converge = 0; 
        double lambda_min = .0, lambda_max = .0, tmp = .0, obj = .0, new_lam = .0; 
        
        for (int i = 0; i < n; ++i)
        {
            lambda_max = fmax(lambda_max , fabs(x[i])); 
        }
        
        for (int i = 0; i < 1e4; ++i)
        {
            new_lam = .5 * (lambda_max + lambda_min);
            // printf("current solution is %f at iteration %d. \n", new_lam, i);
            obj = -z[0];
            for (int j = 0; j < n; ++j)
            {
                obj += c[j] * fmax(fabs(x[j]) - c[j]*new_lam, .0); 
            }
            if (obj > 0)
            {
                lambda_min = new_lam; 
            } else 
            {
                lambda_max = new_lam;
            }
            if (fabs(lambda_max - lambda_min)/(lambda_max+1e-10) < 1e-10)
            {
                // printf("converges in %d steps. \n", i);
                converge = 1; 
                break;
            }
        }
        for (int i = 0; i < n; ++i)
        {
            x[i] = sign(x[i]) * fmax(fabs(x[i]) - c[i]*lambda_max, .0);
        }
    }
}

void soft_threshold_gen(double *x, bool* thred, const int n, const double lam)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        if (thred[i]) x[i] = fmax(0, x[i] - lam) - fmax(0, -x[i] - lam);
    }
}

void soft_threshold(double *x, const int n, const double lam)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i] - lam) - fmax(0, -x[i] - lam);
    }
}

void soft_threshold_vec(double *x, const int n, const double* lam_mat)
{
    int i;
    // #pragma omp parallel for
    for (i = 0; i < n; i++) {
        x[i] = fmax(0, x[i] - lam_mat[i]) - fmax(0, -x[i] - lam_mat[i]);
    }
}

/*    X = eigvec' * eigval * eigvec    */
void eigen_decomp(int n, double* X, double *eigvec, double *eigval) 
{
    
    double *WORK;
    double abstol, WORKopt, vl, vu;
    int *IWORK;
    int numeig, sizeWORK, sizeIWORK, IWORKopt, il, iu,info;
    vl = 0.0;
    vu = 0.0;
    il = 0;
    iu = 0;
    /*  The support of the eigenvectors. We will not use this but the routine needs it  */
    int ISUPPZ[2*n];
    abstol = -1.0; // default tolerance
    
    /*  Query the Lapack routine for optimal sizes for workspace arrays  */
    sizeWORK = -1;
    sizeIWORK = -1;
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, &WORKopt, &sizeWORK, &IWORKopt, &sizeIWORK,&info);
    sizeWORK = (int)WORKopt;
    sizeIWORK = IWORKopt;
    WORK = (double*)malloc (sizeWORK*sizeof(double));
    IWORK = (int*)malloc (sizeIWORK*sizeof(int));
    /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, WORK, &sizeWORK, IWORK, &sizeIWORK,&info);
    /*  Cleanup  */
    free((void*)(WORK)); free((void*)IWORK);
}


void create_tmp_vars(tmpvars** tmp_,int p)
{
    int p_square = p*p;
    tmpvars* tmp    = (tmpvars*) malloc(sizeof(tmpvars));
    tmp->Delta      = (double *) malloc(sizeof(double) * p_square);
    tmp->Delta_old  = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp1       = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp2       = (double *) malloc(sizeof(double) * p_square);
    tmp->Lambda     = (double *) malloc(sizeof(double) * p);
    tmp->lam_mat    = (double *) malloc(sizeof(double) * p_square);
    *tmp_ = tmp;
}

void free_tmp_vars_constrained(tmpvars* tmp){
    if (tmp->Delta != NULL) free(tmp->Delta);
    if (tmp->Delta_old != NULL) free(tmp->Delta_old);
    if (tmp->tmp1 != NULL) free(tmp->tmp1);
    if (tmp->tmp2 != NULL) free(tmp->tmp2);
    if (tmp->Lambda != NULL) free(tmp->Lambda);
    if (tmp->lam_mat  != NULL) free(tmp->lam_mat);
    if (tmp->constrained_idx_matrix) free(tmp->constrained_idx_matrix);
    if (tmp != NULL) free(tmp);
}

void create_tmp_vars_constrained(tmpvars** tmp_,int p)
{
    int p_square = p*p;
    tmpvars* tmp    = (tmpvars*) malloc(sizeof(tmpvars));
    tmp->Delta      = (double *) malloc(sizeof(double) * p_square);
    tmp->Delta_old  = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp1       = (double *) malloc(sizeof(double) * p_square);
    tmp->tmp2       = (double *) malloc(sizeof(double) * p_square);
    tmp->Lambda     = (double *) malloc(sizeof(double) * p);
    tmp->lam_mat    = (double *) malloc(sizeof(double) * p_square);
    tmp->constrained_idx_matrix = (bool *) malloc(sizeof(bool) * p_square);
    *tmp_ = tmp;
}

void free_tmp_vars(tmpvars* tmp){
    if (tmp->Delta != NULL) free(tmp->Delta);
    if (tmp->Delta_old != NULL) free(tmp->Delta_old);
    if (tmp->tmp1 != NULL) free(tmp->tmp1);
    if (tmp->tmp2 != NULL) free(tmp->tmp2);
    if (tmp->Lambda != NULL) free(tmp->Lambda);
    if (tmp->lam_mat  != NULL) free(tmp->lam_mat);
    if (tmp != NULL) free(tmp);
}



bool check_constraint(int p_square, double* omega, const double* tau, const int* bound, tmpvars* tmp)
{
    double constraint = .0; 
    for (int i = 0; i < p_square; ++i)
    {
        if (tmp->constrained_idx_matrix[i]) constraint += fmin(fabs(omega[i]), tau[0]);
    }
    if (constraint > (bound[0]*tau[0]+1e-10))
    {
        // printf("upper bound is: %10.10f, and constraint is: %10.10f. \n", bound[0]*tau[0], constraint);
        return false; 
    } else { return true; }
}


bool check_constraint_convex(int p_square, double* omega, const double* bound, tmpvars* tmp)
{
    double constraint = .0; 
    for (int i = 0; i < p_square; ++i)
    {
        if (tmp->constrained_idx_matrix[i]) constraint += (fabs(omega[i]) * tmp->lam_mat[i]);
    }
    if (constraint > (bound[0]+1e-10))
    {
        // printf("upper bound is: %10.10f, and constraint is: %10.10f. \n", bound[0], constraint);
        return false; 
    } else { return true; }
}


void glasso_gen(const double* Sigma_hat, const int* pp,
                const int* num_iter, const double* eps_abs,
                const double* eps_rel,
                double* Rho, double* Alpha, tmpvars* tmp,
                double* Gamma, double* Omega)
{
    int p = pp[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j;
    double rho = Rho[0], alpha=Alpha[0],primal_res,dual_res,eps_primal,eps_dual;
    double *Delta = tmp->Delta, *Delta_old = tmp->Delta_old, *tmp1 = tmp->tmp1, *tmp2 = tmp->tmp2, *Lambda = tmp->Lambda, *lam_mat=tmp->lam_mat;
    dmat_waxpby(p_square,1.0,Omega,1.0,Gamma,Delta);
    dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
    soft_threshold_vec(Delta,p_square,tmp1);
    
    for (j = 0; j < max_iter; j++) {
        // Omega udpate //
        dmat_waxpby(p_square,rho,Delta,-rho,Gamma,tmp1);
        dmat_waxpby(p_square,1.0,tmp1,-1.0,Sigma_hat,tmp1);
        
        eigen_decomp(p,tmp1,tmp2,Lambda);
        for (i = 0; i < p; i++){ Lambda[i] = (Lambda[i] + sqrt(Lambda[i]*Lambda[i]+4*rho))/(2*rho); }
        // if (isnan_mat(Lambda,p)) { printf("Lambda is nan! \n"); }
        // if (isnan_mat(tmp2,p_square)) { printf("tmp2 is nan! \n"); }
        dmat_B_ATDA(p,p,tmp2,Lambda,Omega);
        // if (isnan_mat(Omega,p_square)) { printf("Omega is nan! \n"); }
        
        // Delta update //
        dmat_vcopy(p_square,Delta,Delta_old);
        dmat_waxpby(p_square,alpha,Omega,1.0-alpha,Delta,Delta);
        dmat_waxpby(p_square,1.0,Delta,1.0,Gamma,Delta);
        dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
        soft_threshold_vec(Delta,p_square,tmp1);
        
        // Gamma update //
        dmat_waxpby(p_square,1.0,Delta_old,-1.0,Delta,tmp1);
        dmat_waxpby(p_square,1.0,Omega,-1.0,Delta,tmp2);
        dmat_waxpby(p_square,1.0,Gamma,1.0-alpha,tmp1,Gamma);
        dmat_waxpby(p_square,1.0,Gamma,alpha,tmp2,Gamma);
        
        // check stopping rule //
        dual_res = rho * dmat_norm2(p_square,tmp1);
        primal_res = dmat_norm2(p_square,tmp2);
        
        eps_primal = eps_abs[0]*p + eps_rel[0]*fmax(dmat_norm2(p_square,Omega),dmat_norm2(p_square,Delta));
        eps_dual = eps_abs[0]*p + eps_rel[0]*dmat_norm2(p_square,Gamma)*rho;
        
        // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
        if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
            dmat_vcopy(p_square,Delta,Omega);
            // printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
            break;
        }
    }
}


void glasso_gen_constrained(const double* Sigma_hat, const int* pp,
                const int* num_iter, const double* eps_abs,
                const double* eps_rel,
                double* Rho, double* Alpha, tmpvars* tmp,
                const double* z, 
                double* Gamma, double* Omega)
{
    int p = pp[0];
    int max_iter = num_iter[0], p_square = p*p;
    int i,j;
    double rho = Rho[0], alpha=Alpha[0],primal_res,dual_res,eps_primal,eps_dual;
    double *Delta = tmp->Delta, *Delta_old = tmp->Delta_old, *tmp1 = tmp->tmp1, *tmp2 = tmp->tmp2, *Lambda = tmp->Lambda, *lam_mat=tmp->lam_mat;
    dmat_waxpby(p_square,1.0,Omega,1.0,Gamma,Delta);
    // dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
    // soft_threshold_vec(Delta,p_square,tmp1);
    proj_L1_gen(Delta, lam_mat, z, &p_square);
    if (tmp->num_of_zero_para[0] != 0) 
    {
        for (int i = 0; i < tmp->num_of_zero_para[0]; ++i)
        {
            Delta[tmp->zero_para_index[2*i] + tmp->zero_para_index[2*i+1] * p] = .0;
            Delta[tmp->zero_para_index[2*i+1] + tmp->zero_para_index[2*i] * p] = .0;
        }
    }
    for (j = 0; j < max_iter; j++) 
    {

        // print_dmatrix(Delta, p, p);
        // print_dmatrix(Omega, p, p);

        // Omega udpate //
        dmat_waxpby(p_square,rho,Delta,-rho,Gamma,tmp1);
        dmat_waxpby(p_square,1.0,tmp1,-1.0,Sigma_hat,tmp1);
        
        eigen_decomp(p,tmp1,tmp2,Lambda);
        for (i = 0; i < p; i++){ Lambda[i] = (Lambda[i] + sqrt(Lambda[i]*Lambda[i]+4*rho))/(2*rho); }
        // if (isnan_mat(Lambda,p)) { printf("Lambda is nan! \n"); }
        // if (isnan_mat(tmp2,p_square)) { printf("tmp2 is nan! \n"); }
        dmat_B_ATDA(p,p,tmp2,Lambda,Omega);
        // if (isnan_mat(Omega,p_square)) { printf("Omega is nan! \n"); }
        
        // Delta update //
        dmat_vcopy(p_square,Delta,Delta_old);
        dmat_waxpby(p_square,alpha,Omega,1.0-alpha,Delta,Delta);
        dmat_waxpby(p_square,1.0,Delta,1.0,Gamma,Delta);
        // dmat_waxpby(p_square,1.0/rho,lam_mat,.0,tmp1,tmp1);
        // soft_threshold_vec(Delta,p_square,tmp1);
        proj_L1_gen(Delta, lam_mat, z, &p_square);

        // print_dmatrix(Delta, p, p);
        // print_dmatrix(Omega, p, p);


        if (tmp->num_of_zero_para[0] != 0) 
        { 
            for (int i = 0; i < tmp->num_of_zero_para[0]; ++i)
            {
                Delta[tmp->zero_para_index[2*i] + tmp->zero_para_index[2*i+1] * p] = .0;
                Delta[tmp->zero_para_index[2*i+1] + tmp->zero_para_index[2*i] * p] = .0;
            }
        }
        
        // Gamma update //
        dmat_waxpby(p_square,1.0,Delta_old,-1.0,Delta,tmp1);
        dmat_waxpby(p_square,1.0,Omega,-1.0,Delta,tmp2);
        dmat_waxpby(p_square,1.0,Gamma,1.0-alpha,tmp1,Gamma);
        dmat_waxpby(p_square,1.0,Gamma,alpha,tmp2,Gamma);
        
        // print_dmatrix(Delta, p, p);
        // print_dmatrix(Omega, p, p);

        // check stopping rule //
        dual_res = rho * dmat_norm2(p_square,tmp1);
        primal_res = dmat_norm2(p_square,tmp2);
        
        eps_primal = eps_abs[0]*p + eps_rel[0]*fmax(dmat_norm2(p_square,Omega),dmat_norm2(p_square,Delta));
        eps_dual = eps_abs[0]*p + eps_rel[0]*dmat_norm2(p_square,Gamma)*rho;
        
        // if (j % 10 == 0) printf("#%d, primal res: %5.4f, dual res: %5.4f. \n", j, primal_res, dual_res);
        // error("stop here");
        if ((primal_res < eps_primal) && (dual_res < eps_dual)) {
            dmat_vcopy(p_square,Delta,Omega);
            // printf("#%d, primal res: %5.10f, dual res: %5.10f. \n", j, primal_res, dual_res);
            break;
        }

        if (j == (max_iter - 1))
        {
            printf("warning: ADMM does not converge! Please increase the number of ADMM iterations. \n");
            printf("current primal res: %5.10f, dual res: %5.10f. \n", primal_res, dual_res);
            dmat_vcopy(p_square,Delta,Omega);
            break;
        }
    }
}




void glasso_L1(const double* Sigma_hat, const double* lam,
               const int* pp,
               const int* num_iter, const double* eps_abs,
               const double* eps_rel,
               double* Rho, double* Alpha, double* Omega)
{
    int i,p_square = pp[0]*pp[0];
    double* Gamma = (double *) malloc(sizeof(double) * p_square);
    dmat_vset(p_square,.0,Gamma);
    tmpvars* tmp;
    create_tmp_vars(&tmp,pp[0]);
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma,Omega);
    if (Gamma != NULL) free(Gamma);
    free_tmp_vars(tmp);
}


// calculate constrained MLE, entries outside the set of free_para_index are set to zero //
void glasso_cMLE(const double* Sigma_hat,
                 const int* free_para_index,
                 const int* num_of_free_para,
                 const int* pp,
                 const int* num_iter, const double* eps_abs,
                 const double* eps_rel,
                 double* Rho, double* Alpha, double* Omega)
{
    int i,p_square = pp[0]*pp[0];
    double* Gamma = (double *) malloc(sizeof(double) * p_square);
    dmat_vset(p_square,.0,Gamma);
    tmpvars* tmp;
    create_tmp_vars(&tmp,pp[0]);
    dmat_vset(p_square,1e10,tmp->lam_mat);
    for (i = 0; i < num_of_free_para[0]; i++) {
        tmp->lam_mat[free_para_index[i]] = .0;
    }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma,Omega);
    if (Gamma != NULL) free(Gamma);
    free_tmp_vars(tmp);
}

void glasso_nonconvex(const double* Sigma_hat, const double* lam,
                      const double* tau,
                      const int* pp,
                      const int* num_iter, const int* dc_max_iter,
                      const double* eps_abs,
                      const double* eps_rel,
                      double* Rho, double* Alpha, tmpvars* tmp,
                      enum penalty_types* pen_type,
                      double* Gamma_L1, double* Gamma_trL1,
                      double* Omega_L1, double* Omega_trL1)
{
    int i,j,p_square = pp[0]*pp[0];
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type[0] == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - abs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (abs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}


// minimize <Omega, Sigma_hat> - log det Omega 
// subject to \sum_{ij} p_{tau}(|Omega_{ij}|) < z. 
// initialized by convex solution // 
void glasso_nonconvex_constrained(const double* Sigma_hat, const int* z,
                      const double* tau,
                      const int* pp,
                      const int* num_iter, const int* dc_max_iter,
                      const double* eps_abs,
                      const double* eps_rel,
                      double* Rho, double* Alpha, tmpvars* tmp,
                      double* Gamma_L1, double* Gamma_trL1,
                      double* Omega_L1, double* Omega_trL1)
{
    int i,j,p_square = pp[0]*pp[0];
    int bound = z[0];
    dmat_vset(p_square,1.0,tmp->lam_mat);
    for (i = 0; i < p_square; i++)
    {
        if (!tmp->constrained_idx_matrix[i]) tmp->lam_mat[i] = .0; 
    }
    double bound_times_tau = bound * tau[0]; 
    glasso_gen_constrained(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&bound_times_tau,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);

    for (j = 0; j < dc_max_iter[0]; j++)
    {
        // printf("DC iteration: %d. \n", j+1);
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1); // store previous lam_mat //

        bound = z[0];
        for (i = 0; i < p_square; i++)
        {
            if (tmp->constrained_idx_matrix[i]) 
            {
                if (fabs(Omega_trL1[i]) <= tau[0]) 
                {
                    tmp->lam_mat[i] = 1.0;
                } else 
                {
                    tmp->lam_mat[i] = .0;
                    bound--;
                }
            } else 
            {
                tmp->lam_mat[i] = .0;
            }
        }
        // printf("upper bound is: %d. \n", bound);
        if (bound < 0) {
            printf("upper bound is: %d. \n", bound);
            error("upper bound becomes negative! \n");
        }
        
        // print_dmatrix(tmp->lam_mat, pp[0], pp[0]);

        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < eps_abs[0]) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        bound_times_tau = bound * tau[0];
        glasso_gen_constrained(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&bound_times_tau,Gamma_trL1,Omega_trL1);
    }
    
    if(!check_constraint(p_square,Omega_trL1, tau, z, tmp)) { 
        printf("tau is %f. \n", tau[0]);
        printf("upper bound is %d. \n", z[0]);
        print_dmatrix(Omega_trL1, pp[0],pp[0]);
        error("nonconvex solution not feasible! \n"); 
    }
}



// single parameter with given initialization // 
void glasso_nonconvex_constrained_single(const double* Sigma_hat, 
                                const int* free_para_index, const int* num_of_free_para, 
                                const int* zero_para_index, const int* num_of_zero_para,
                                const int* bound, const double* tau, const int* pp,
                               const int* num_iter, const int* dc_max_iter, const double* eps_abs,
                               const double* eps_rel, double* Rho, double* Alpha,
                               double* Omega_trL1)
{
    int i,j, p = pp[0], p_square = pp[0]*pp[0]; tmpvars* tmp;
    int bound_tmp = 0;
    // double* Gamma_L1 = (double *) malloc(sizeof(double) * p * p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p * p);
    dmat_vset(p*p,.0,Gamma_trL1);
    create_tmp_vars_constrained(&tmp,p);

    tmp->zero_para_index = zero_para_index; 
    tmp->num_of_zero_para = num_of_zero_para; 

    for (i = 0; i < p*p; ++i) tmp->constrained_idx_matrix[i] = true;
    for (i = 0; i < p; ++i) tmp->constrained_idx_matrix[i*p + i] = false; 
    for (i = 0; i < num_of_free_para[0]; ++i)
    {
        tmp->constrained_idx_matrix[free_para_index[2*i] + p*free_para_index[2*i+1]] = false; 
        tmp->constrained_idx_matrix[free_para_index[2*i+1] + p*free_para_index[2*i]] = false; 
    }

    // print_bmatrix(tmp->constrained_idx_matrix,p,p);
    for (j = 0; j < dc_max_iter[0]; j++)
    {
        // printf("DC iteration: %d. \n", j+1);
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        bound_tmp = bound[0];
        // printf("bound is %d. \n", bound[0]);
        // printf("upper bound (before) is %d. \n", bound_tmp);
        for (i = 0; i < p_square; i++)
        {
            if (tmp->constrained_idx_matrix[i]) 
            {
                if (fabs(Omega_trL1[i]) <= tau[0]) 
                {
                    tmp->lam_mat[i] = 1.0;
                }
                else 
                {
                    tmp->lam_mat[i] = .0;
                    bound_tmp -= 1;
                }
            } else 
            {
                tmp->lam_mat[i] = .0;
            }
        }
        // print_dmatrix(tmp->lam_mat,p,p);
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        // printf("upper bound (after) is %d \n", bound_tmp);
        if (bound_tmp < 0) { 
            // printf("upper bound is: %d. \n", bound_tmp); 
            error("upper bound becomes negative! \n"); 

        }
        if (dmat_norminf(p_square,tmp->tmp1) < eps_abs[0] && j >= 1) 
        {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        double z = bound_tmp*tau[0];
        glasso_gen_constrained(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&z,Gamma_trL1,Omega_trL1);
    }

    if(!check_constraint(p_square,Omega_trL1, tau, bound, tmp)) 
    { 
        printf("tau is %f. \n", tau[0]);
        print_dmatrix(Omega_trL1, pp[0],pp[0]);
        error("nonconvex solution not feasible! \n"); 
    }
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars_constrained(tmp);
}


bool make_sol_feasible(int p, double* omega, const double* tau_prev, const double* tau_cur, const int* bound, tmpvars* tmp)
{
    double constraint = .0; 
    double omega_ij = .0;
    int p_square = p*p;
    /*
    // compute feasibility gap // 
    for (int i = 0; i < p_square; ++i)
    {
        if (tmp->constrained_idx_matrix[i]) constraint += fmin(fabs(omega[i])/tau_cur[0], 1.0);
    }
    constraint = constraint - bound[0]; 
    if (constraint <= .0) { return true; }
    */
    for (int i = 0; i < p; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            if (tmp->constrained_idx_matrix[i*p + j]) 
            {
                omega_ij = fabs(omega[i*p + j]);
                if (omega_ij <= tau_prev[0])
                {
                    omega[i*p + j] = omega[i*p + j] * tau_cur[0] / tau_prev[0];
                    omega[j*p + i] = omega[j*p + i] * tau_cur[0] / tau_prev[0];
                }
            }
        }
    }
    if (check_constraint(p_square, omega, tau_cur, bound, tmp)) { return true; } 
    else 
    { 
        // print_dmatrix(omega,p,p);
        printf("current tau is %f. \n",tau_cur[0]); error("can't make the solution feasible! \n");    
    }
    /*
    for (int i = 0; i < p; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            if (tmp->constrained_idx_matrix[i*p + j]) 
            {
                omega_ij = fabs(omega[i*p + j]);
                if (omega_ij <= tau_prev[0] && omega_ij >= tau_cur[0])
                {
                    // constraint decreasing amount //
                    omega_ij = 2.0 * (fmin(fabs(omega[i*p+j])/tau_cur[0], 1.0))
                    omega[i*p + j] = tau_cur[0];
                    omega[j*p + i] = tau_cur[0];
                }
            }
        }
    }*/
}

void glasso_nonconvex_constrained_path(const double* Sigma_hat, 
                                const int* free_para_index, const int* num_of_free_para, 
                                const int* zero_para_index, const int* num_of_zero_para,
                                const int* bound, const int* bound_length, 
                                const double* tau, const int* tau_length,
                                const int* pp, const int* num_iter, const int* dc_max_iter,
                               const double* eps_abs, const double* eps_rel,
                               double* Rho, double* Alpha,
                               double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p * p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p * p);
    create_tmp_vars_constrained(&tmp,p);
    
    tmp->zero_para_index = zero_para_index; 
    tmp->num_of_zero_para = num_of_zero_para; 
    for (int i = 0; i < p*p; ++i) tmp->constrained_idx_matrix[i] = true;
    for (int i = 0; i < p; ++i) tmp->constrained_idx_matrix[i*p + i] = false; 
    for (int i = 0; i < num_of_free_para[0]; ++i)
    {
        tmp->constrained_idx_matrix[free_para_index[2*i] + p*free_para_index[2*i+1]] = false; 
        tmp->constrained_idx_matrix[free_para_index[2*i+1] + p*free_para_index[2*i]] = false; 
    }
    
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < bound_length[0]; i++)
    {
        // printf("bound is: %d.\n", bound[i]);
        // get solution for maximum tau //
        if (i == 0)
        {
            dmat_vset(p*p,.0,Gamma_L1);
            glasso_nonconvex_constrained(Sigma_hat,bound,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else 
        {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            glasso_nonconvex_constrained(Sigma_hat,bound+i,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_length[0]; j++)
        {
            // printf("tau[%d] is %f. \n", j, tau[j]);
            // Rho[0] = bound[i]; // adaptively varying rho //
            // glasso_nonconvex_constrained(Sigma_hat,bound+i,tau+j,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
            dmat_vcopy(p*p,Omega_trL1+i*tau_length[0]*p*p+(j-1)*p*p,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
            
            // start check feasibility //
            dmat_vcopy(p*p,Omega_trL1+i*tau_length[0]*p*p+(j-1)*p*p,tmp->tmp1);
            
            // if (!check_constraint(p*p, tmp->tmp1, tau+j-1, bound, tmp)) { error("not feasible at previous tau. \n"); }

            // decreasing tau may lead to the event that the previous solution becomes infeasible. 
            if (!check_constraint(p*p, tmp->tmp1, tau+j, bound+i, tmp)) 
            {
                // printf("not feasible at current tau. \n");
                // print_dmatrix(tmp->tmp1, p, p);
                // if (!check_constraint(p*p, tmp->tmp1, tau+j-1, bound+i, tmp)) error("not feasible at previous tau! \n");
                make_sol_feasible(p, tmp->tmp1, tau+j-1, tau+j, bound+i, tmp);
            }
            dmat_vcopy(p*p,tmp->tmp1,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
            // end  checking feasibility //
            // printf("tau[%d] is %f. \n", j, tau[j]);
            // print_dmatrix(tmp->tmp1, p, p);
            // printf("bound is %d. \n", bound[0]);
            glasso_nonconvex_constrained_single(Sigma_hat,free_para_index,num_of_free_para,zero_para_index,num_of_zero_para,bound+i, tau+j, pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars_constrained(tmp);
}




void glasso_nonconvex_path(const double* Sigma_hat, const double* lam,
                           const int* lam_length,const double* tau,
                           const int* pp,
                           const int* num_iter, const int* dc_max_iter,
                           const double* eps_abs,
                           const double* eps_rel,
                           double* Rho, double* Alpha,
                           const int* method,
                           double* Omega_L1, double* Omega_trL1)
{
    int i, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; printf("TLP is used. \n");}
    else if (method[0] == 2) { pen_type = MCP; printf("MCP is used. \n");}
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++)
    {
        // Rho[0] = lam[i]; // adaptively varying rho //
        printf("lambda is: %f\n", lam[i]);
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            glasso_nonconvex(Sigma_hat,lam,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            glasso_nonconvex(Sigma_hat,lam+i,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}


void glasso_nonconvex_path_new(const double* Sigma_hat, const double* lam,
                               const int* lam_length,const double* tau,
                               const int* tau_length,const int* pp,
                               const int* num_iter, const int* dc_max_iter,
                               const double* eps_abs,
                               const double* eps_rel,
                               double* Rho, double* Alpha,
                               const int* method,
                               double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) {
        pen_type = Truncated_L1;
        // printf("TLP is used. \n");
    }
    else if (method[0] == 2) {
        pen_type = MCP;
        // printf("MCP is used. \n");
    }
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        // printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            glasso_nonconvex(Sigma_hat,lam,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            glasso_nonconvex(Sigma_hat,lam+i,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            glasso_nonconvex(Sigma_hat,lam+i,tau+j,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}


void glasso_nonconvex_path_v2(const double* Sigma_hat, const double* lam,
                              const int* lam_length,const double* tau_div_lam,
                              const int* tau_div_lam_length,const int* pp,
                              const int* num_iter, const int* dc_max_iter,
                              const double* eps_abs,
                              const double* eps_rel,
                              double* Rho, double* Alpha,
                              const int* method,
                              double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; printf("TLP is used. \n");}
    else if (method[0] == 2) { pen_type = MCP; printf("MCP is used. \n");}
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    double tau;
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            tau = tau_div_lam[0] * lam[0];
            glasso_nonconvex(Sigma_hat,lam,&tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            tau = tau_div_lam[0] * lam[i];
            glasso_nonconvex(Sigma_hat,lam+i,&tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_div_lam_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            tau = tau_div_lam[j] * lam[i];
            glasso_nonconvex(Sigma_hat,lam+i,&tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_div_lam_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}




// set some parameters to zero //
void glasso_nonconvex_zero(const double* Sigma_hat, const double* lam,
                                const int* zero_para_index, const int* num_of_zero_para,
                                const double* tau,
                                const int* pp,
                                const int* num_iter, const int* dc_max_iter,
                                const double* eps_abs,
                                const double* eps_rel,
                                double* Rho, double* Alpha, const int* method,
                                double* Gamma_L1, double* Gamma_trL1,
                                double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0],p_square = pp[0]*pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) {
        pen_type = Truncated_L1;
        // printf("TLP is used. \n");
    }
    else if (method[0] == 2) {
        pen_type = MCP;
        // printf("MCP is used. \n");
    } else { error("nonconvex penlaty not supported! \n"); }
    create_tmp_vars(&tmp,p);
    
    
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    for (i = 0; i < num_of_zero_para[0]; i++) {
        tmp->lam_mat[zero_para_index[2*i]*pp[0]+zero_para_index[2*i+1]] = 1e10;
        tmp->lam_mat[zero_para_index[2*i+1]*pp[0]+zero_para_index[2*i]] = 1e10;
    }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - abs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (abs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        for (i = 0; i < num_of_zero_para[0]; i++) {
            tmp->lam_mat[zero_para_index[2*i]*pp[0]+zero_para_index[2*i+1]] = 1e10;
            tmp->lam_mat[zero_para_index[2*i+1]*pp[0]+zero_para_index[2*i]] = 1e10;
        }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
    free_tmp_vars(tmp);
}


// do not penalize some free parameters //
void glasso_nonconvex_free(const double* Sigma_hat, const double* lam,
                                const int* free_para_index, const int* num_of_free_para,
                                const double* tau,
                                const int* pp,
                                const int* num_iter, const int* dc_max_iter,
                                const double* eps_abs,
                                const double* eps_rel,
                                double* Rho, double* Alpha,const int* method,
                                double* Gamma_L1, double* Gamma_trL1,
                                double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0],p_square = pp[0]*pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) {
        pen_type = Truncated_L1;
        // printf("TLP is used. \n");
    }
    else if (method[0] == 2) {
        pen_type = MCP;
        // printf("MCP is used. \n");
    } else { error("nonconvex penlaty not supported! \n"); }
    create_tmp_vars(&tmp,p);
    
    
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    for (i = 0; i < num_of_free_para[0]; i++) {
        tmp->lam_mat[free_para_index[2*i]*pp[0]+free_para_index[2*i+1]] = .0;
        tmp->lam_mat[free_para_index[2*i+1]*pp[0]+free_para_index[2*i]] = .0;
    }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - abs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (abs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        for (i = 0; i < num_of_free_para[0]; i++) {
            tmp->lam_mat[free_para_index[2*i]*pp[0]+free_para_index[2*i+1]] = .0;
            tmp->lam_mat[free_para_index[2*i+1]*pp[0]+free_para_index[2*i]] = .0;
        }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
    free_tmp_vars(tmp);
}
                                


// inference part: do not penalize some free parameters //


// asymptotic normality based CI //
void glasso_nonconvex_inference(const double* Sigma_hat, const double* lam,
                                const int* free_para_index, const int* num_of_free_para,
                                const double* tau,
                                const int* pp,
                                const int* num_iter, const int* dc_max_iter,
                                const double* eps_abs,
                                const double* eps_rel,
                                double* Rho, double* Alpha, tmpvars* tmp,
                                enum penalty_types* pen_type,
                                double* Gamma_L1, double* Gamma_trL1,
                                double* Omega_L1, double* Omega_trL1)
{
    int i,j,p_square = pp[0]*pp[0];
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    for (i = 0; i < num_of_free_para[0]; i++) {
        tmp->lam_mat[free_para_index[2*i]*pp[0]+free_para_index[2*i+1]] = .0;
        tmp->lam_mat[free_para_index[2*i+1]*pp[0]+free_para_index[2*i]] = .0;
    }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type[0] == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - abs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (abs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        for (i = 0; i < num_of_free_para[0]; i++) {
            tmp->lam_mat[free_para_index[2*i]*pp[0]+free_para_index[2*i+1]] = .0;
            tmp->lam_mat[free_para_index[2*i+1]*pp[0]+free_para_index[2*i]] = .0;
        }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}

// asymptotic normality based CI //
void glasso_nonconvex_path_inference(const double* Sigma_hat, const double* lam,
                                     const int* lam_length,
                                     const int* free_para_index, const int* num_of_free_para,
                                     const double* tau_div_lam,
                                     const int* tau_div_lam_length,const int* pp,
                                     const int* num_iter, const int* dc_max_iter,
                                     const double* eps_abs,
                                     const double* eps_rel,
                                     double* Rho, double* Alpha,
                                     const int* method,
                                     double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; printf("TLP is used. \n");}
    else if (method[0] == 2) { pen_type = MCP; printf("MCP is used. \n");}
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    double tau;
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            tau = tau_div_lam[0] * lam[0];
            glasso_nonconvex_inference(Sigma_hat,lam,free_para_index,num_of_free_para,&tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            tau = tau_div_lam[0] * lam[i];
            glasso_nonconvex_inference(Sigma_hat,lam+i,free_para_index, num_of_free_para, &tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_div_lam_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            tau = tau_div_lam[j] * lam[i];
            glasso_nonconvex_inference(Sigma_hat,lam+i,free_para_index,num_of_free_para, &tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_div_lam_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}



// Likelihood ratio based CI //
void glasso_nonconvex_inference_LR(const double* Sigma_hat, const double* lam,
                                const int* zero_para_index, const int* num_of_zero_para,
                                const double* tau,
                                const int* pp,
                                const int* num_iter, const int* dc_max_iter,
                                const double* eps_abs,
                                const double* eps_rel,
                                double* Rho, double* Alpha, tmpvars* tmp,
                                enum penalty_types* pen_type,
                                double* Gamma_L1, double* Gamma_trL1,
                                double* Omega_L1, double* Omega_trL1)
{
    int i,j,p_square = pp[0]*pp[0];
    dmat_vset(p_square,lam[0],tmp->lam_mat);
    for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
    for (i = 0; i < num_of_zero_para[0]; i++) {
        tmp->lam_mat[zero_para_index[2*i]*pp[0]+zero_para_index[2*i+1]] = 1e10;
        tmp->lam_mat[zero_para_index[2*i+1]*pp[0]+zero_para_index[2*i]] = 1e10;
    }
    glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_L1,Omega_L1);
    dmat_vcopy(p_square,Omega_L1,Omega_trL1);
    dmat_vcopy(p_square,Gamma_L1,Gamma_trL1);
    for (j = 0; j < dc_max_iter[0]; j++) {
        dmat_vcopy(p_square,tmp->lam_mat,tmp->tmp1);
        if (pen_type[0] == MCP) { for (i = 0; i < p_square; i++) { tmp->lam_mat[i] = lam[0]*fmax(.0, 1.0 - abs(Omega_trL1[i])/(2*tau[0])); } }
        if (pen_type[0] == Truncated_L1)
        {
            for (i = 0; i < p_square; i++)
            {
                if (abs(Omega_trL1[i]) < tau[0]) tmp->lam_mat[i] = lam[0];
                else tmp->lam_mat[i] = .0;
            }
        }
        for (i = 0; i < pp[0]; i++) { tmp->lam_mat[i*pp[0]+i] = .0; }
        for (i = 0; i < num_of_zero_para[0]; i++) {
            tmp->lam_mat[zero_para_index[2*i]*pp[0]+zero_para_index[2*i+1]] = 1e10;
            tmp->lam_mat[zero_para_index[2*i+1]*pp[0]+zero_para_index[2*i]] = 1e10;
        }
        dmat_waxpby(p_square,1.0,tmp->lam_mat,-1.0,tmp->tmp1,tmp->tmp1);
        // print_dmatrix(tmp->lam_mat,pp[0],pp[0]);
        // print_dmatrix(Omega_trL1,pp[0],pp[0]);
        // printf("DC gap: %5.20f. \n",dmat_norminf(p_square,tmp->tmp1));
        if (dmat_norminf(p_square,tmp->tmp1) < (lam[0]*eps_abs[0])) {
            // printf("DC iteration used: %d. \n", j);
            break;
        }
        glasso_gen(Sigma_hat,pp,num_iter,eps_abs,eps_rel,Rho,Alpha,tmp,Gamma_trL1,Omega_trL1);
    }
}

// Likelihood ratio based CI //
void glasso_nonconvex_path_inference_LR(const double* Sigma_hat, const double* lam,
                                     const int* lam_length,
                                     const int* zero_para_index, const int* num_of_zero_para,
                                     const double* tau_div_lam,
                                     const int* tau_div_lam_length,const int* pp,
                                     const int* num_iter, const int* dc_max_iter,
                                     const double* eps_abs,
                                     const double* eps_rel,
                                     double* Rho, double* Alpha,
                                     const int* method,
                                     double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; }
    else if (method[0] == 2) { pen_type = MCP; }
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    double tau;
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        // printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            tau = tau_div_lam[0] * lam[0];
            glasso_nonconvex_inference_LR(Sigma_hat,lam,zero_para_index,num_of_zero_para,&tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            tau = tau_div_lam[0] * lam[i];
            glasso_nonconvex_inference_LR(Sigma_hat,lam+i,zero_para_index,num_of_zero_para, &tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_div_lam_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            tau = tau_div_lam[j] * lam[i];
            glasso_nonconvex_inference_LR(Sigma_hat,lam+i,zero_para_index,num_of_zero_para, &tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_div_lam_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}


void glasso_nonconvex_path_inference_LR_v2(const double* Sigma_hat, const double* lam,
                                        const int* lam_length,
                                        const int* zero_para_index, const int* num_of_zero_para,
                                        const double* tau,
                                        const int* tau_length,const int* pp,
                                        const int* num_iter, const int* dc_max_iter,
                                        const double* eps_abs,
                                        const double* eps_rel,
                                        double* Rho, double* Alpha,
                                        const int* method,
                                        double* Omega_L1, double* Omega_trL1)
{
    int i,j, p = pp[0]; tmpvars* tmp;
    enum penalty_types pen_type;
    if (method[0] == 1) { pen_type = Truncated_L1; }
    else if (method[0] == 2) { pen_type = MCP; }
    else { error("nonconvex penlaty not supported! \n"); }
    double* Gamma_L1 = (double *) malloc(sizeof(double) * p*p);
    double* Gamma_trL1 = (double *) malloc(sizeof(double) * p*p);
    create_tmp_vars(&tmp,p);
    // double* tmp = (double *) malloc(sizeof(double) * p * p);
    for (i = 0; i < lam_length[0]; i++){
        // printf("lambda is: %f.\n", lam[i]);
        // get solution for maximum tau //
        if (i == 0) {
            dmat_vset(p*p,.0,Gamma_L1);
            glasso_nonconvex_inference_LR(Sigma_hat,lam,zero_para_index,num_of_zero_para,tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1,Omega_trL1);
        } else {
            dmat_vcopy(p*p,Omega_L1+(i-1)*p*p,Omega_L1+i*p*p);
            glasso_nonconvex_inference_LR(Sigma_hat,lam+i,zero_para_index,num_of_zero_para, tau,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*p*p);
        }
        for (j = 1; j < tau_length[0]; j++){
            // Rho[0] = lam[i]; // adaptively varying rho //
            glasso_nonconvex_inference_LR(Sigma_hat,lam+i,zero_para_index,num_of_zero_para, tau+j,pp,num_iter,dc_max_iter,eps_abs,eps_rel,Rho,Alpha,tmp,&pen_type,Gamma_L1,Gamma_trL1,Omega_L1+i*p*p,Omega_trL1+i*tau_length[0]*p*p+j*p*p);
        }
    }
    if (Gamma_L1 != NULL) free(Gamma_L1);
    if (Gamma_trL1 != NULL) free(Gamma_trL1);
    free_tmp_vars(tmp);
}

