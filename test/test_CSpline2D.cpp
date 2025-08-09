#include<stdio.h>
#include<math.h>
#include<string>
#ifdef _OPENMP
#include<omp.h>
#endif
#include"InterCSpline2D.h"

/* number of bins in x
 * where the function is interpolated */
int nbin_x = 16;
/* range in x
 * xmin < x < xmax */
double xmin = 0.;
double xmax = 1.;

/* number of bins in y
 * where the function is interpolated */
int nbin_y = 16;
/* range in y
 * ymin < y < ymax */
double ymin = 0.;
double ymax = 1.;

/* number of points in x and y
 * at which the interpolated function is evaluated */
int n_pt_x = 128;
int n_pt_y = 128;

// function to be interpolated
double func_test(double x, double y,
                 double *df_dx = NULL,
                 double *df_dy = NULL);

int main(int argc, char *argv[]) {
    // x bin
    double *tab_x = new double[nbin_x + 1];
    for (int ix = 0; ix <= nbin_x; ix++) {
        tab_x[ix] = xmin +
            (xmax - xmin) * static_cast<double>(ix) /
                            static_cast<double>(nbin_x);
    }

    // y bin
    double *tab_y = new double[nbin_y + 1];
    for (int iy = 0; iy <= nbin_y; iy++) {
        tab_y[iy] = ymin +
            (ymax - ymin) * static_cast<double>(iy) /
                            static_cast<double>(nbin_y);
    }

    // tabulated function
    double **tab_f =
        InterCSpline2D::new_array_func(nbin_x, nbin_y);
    // boundary condition for the first derivative
    double **tab_bc_df_dx =
        InterCSpline2D::new_array_bc_df_dx(nbin_y);
    double **tab_bc_df_dy =
        InterCSpline2D::new_array_bc_df_dy(nbin_x);

    for (int ix = 0; ix <= nbin_x; ix++) {
        for (int iy = 0; iy <= nbin_y; iy++) {
            double *ptr_df_dx;
            if (ix == 0) {
                ptr_df_dx = &tab_bc_df_dx[iy][0];
            } else if (ix == nbin_x) {
                ptr_df_dx = &tab_bc_df_dx[iy][1];
            } else {
                ptr_df_dx = NULL;
            }

            double *ptr_df_dy;
            if (iy == 0) {
                ptr_df_dy = &tab_bc_df_dy[ix][0];
            } else if (iy == nbin_y) {
                ptr_df_dy = &tab_bc_df_dy[ix][1];
            } else {
                ptr_df_dy = NULL;
            }

            tab_f[ix][iy] =
                func_test(tab_x[ix], tab_y[iy],
                          ptr_df_dx,
                          ptr_df_dy);
        }
    }

    // initialize the interpolator
    InterCSpline2D csp_test;
    csp_test.init(nbin_x, nbin_y,
                  tab_x, tab_y,
                  tab_f,
                  tab_bc_df_dx,
                  tab_bc_df_dy);

    InterCSpline2D::del_array_func(nbin_x, nbin_y,
                                   tab_f);
    InterCSpline2D::del_array_bc_df_dx(nbin_y,
                                       tab_bc_df_dx);
    InterCSpline2D::del_array_bc_df_dy(nbin_x,
                                       tab_bc_df_dy);

    delete [] tab_x;
    delete [] tab_y;

    double *pt_x = new double[n_pt_x + 1];
    for (int ix = 0; ix <= n_pt_x; ix++) {
        pt_x[ix] = xmin +
            (xmax - xmin) * static_cast<double>(ix) /
                            static_cast<double>(n_pt_x);
    }
    double *pt_y = new double[n_pt_y + 1];
    for (int iy = 0; iy <= n_pt_y; iy++) {
        pt_y[iy] = ymin +
            (ymax - ymin) * static_cast<double>(iy) /
                            static_cast<double>(n_pt_y);
    }

    double **pt_f_fin =
        InterCSpline2D::new_array_func(n_pt_x, n_pt_y);
    double **pt_f_lin =
        InterCSpline2D::new_array_func(n_pt_x, n_pt_y);

    // evaluate the interpolated function
    #ifdef _OPENMP
    #pragma omp parallel
    {  // parallel code begins
    #endif
        #ifdef _OPENMP
        int n_thread = omp_get_num_threads();
        int tid = omp_get_thread_num();
        #endif

        for (int ix = 0; ix <= n_pt_x; ix++) {
            #ifdef _OPENMP
            if (ix % n_thread != tid) {
                continue;
            }
            #endif

            for (int iy = 0; iy <= n_pt_y; iy++) {
                pt_f_fin[ix][iy] =
                    csp_test.get_func(pt_x[ix], pt_y[iy]);
                pt_f_lin[ix][iy] =
                    csp_test.get_func_lin(pt_x[ix], pt_y[iy]);
            }
        }
    #ifdef _OPENMP
    }  // parallel code ends
    #endif

    FILE *ptr_fout = stderr;
    bool have_out_ext = argc > 1;
    if (have_out_ext) {
        std::string name_file_out = argv[1];
        ptr_fout = fopen(name_file_out.c_str(), "w");
    }

    if (ptr_fout == NULL) {
        have_out_ext = false;
        ptr_fout = stderr;
    }

    for (int ix = 0; ix <= n_pt_x; ix++) {
        for (int iy = 0; iy <= n_pt_y; iy++) {
            double f_ini = func_test(pt_x[ix], pt_y[iy]);

            double f_fin = pt_f_fin[ix][iy];
            double f_lin = pt_f_lin[ix][iy];

            fprintf(ptr_fout, "    %e    %e    %e",
                    pt_x[ix], pt_y[iy], f_ini);
            fprintf(ptr_fout, "    %e    %e",
                    f_fin, f_fin - f_ini);
            fprintf(ptr_fout, "    %e    %e\n",
                    f_lin, f_lin - f_ini);
        }
        fprintf(ptr_fout, "\n");
    }

    if (have_out_ext) {
        fclose(ptr_fout);
    }

    delete [] pt_x;
    delete [] pt_y;

    InterCSpline2D::del_array_func(n_pt_x, n_pt_y,
                                   pt_f_fin);
    InterCSpline2D::del_array_func(n_pt_x, n_pt_y,
                                   pt_f_lin);

    return 0;
}

double func_test(double x, double y,
                 double *df_dx,
                 double *df_dy) {
    double fac_gauss =
        exp(-0.5 * (x - 0.5) * (x - 0.5) / 0.64) *
        exp(-0.5 * (y - 0.5) * (y - 0.5) / 0.36);
    double f = fac_gauss *
        cos(4. * M_PI * x) * sin(2. * M_PI * y);

    if (df_dx != NULL) {
        *df_dx = fac_gauss *
            (((0.5 - x) / 0.64) * cos(4. * M_PI * x) -
             4. * M_PI * sin(4. * M_PI * x)) *
            sin(2. * M_PI * y);
    }

    if (df_dy != NULL) {
        *df_dy = fac_gauss *
            (((0.5 - y) / 0.36) * sin(2. * M_PI * y) +
             2. * M_PI * cos(2. * M_PI * y)) *
            cos(4. * M_PI * x);
    }

    return f;
}
