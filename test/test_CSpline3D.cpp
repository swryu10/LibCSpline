#include<stdio.h>
#include<math.h>
#include<string>
#ifdef _OPENMP
#include<omp.h>
#endif
#include"InterCSpline2D.h"
#include"InterCSpline3D.h"

/* number of bins in z
 * where the function is interpolated */
int nbin_z = 16;
/* range in z
 * zmin < x < zmax */
double zmin = 0.;
double zmax = 1.;

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

/* number of points in z, x and y
 * at which the interpolated function is evaluated */
int n_pt_z = 64;
int n_pt_x = 64;
int n_pt_y = 64;

/* values of z, x and y
 * at which 2D outputs are produced */
double z_mid = 0.5;
double x_mid = 0.5;
double y_mid = 0.25;

double func_test(double z, double x, double y,
                 double *df_dz = NULL,
                 double *df_dx = NULL,
                 double *df_dy = NULL);

int main(int argc, char *argv[]) {
    // z bin
    double *tab_z = new double[nbin_z + 1];
    for (int iz = 0; iz <= nbin_z; iz++) {
        tab_z[iz] = zmin +
            (zmax - zmin) * static_cast<double>(iz) /
                            static_cast<double>(nbin_z);
    }

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
    double ***tab_f =
        InterCSpline3D::new_array_func(nbin_z, nbin_x, nbin_y);
    // boundary condition for the first derivative
    double ***tab_bc_df_dz =
        InterCSpline3D::new_array_bc_df_dz(nbin_x, nbin_y);
    double ***tab_bc_df_dx =
        InterCSpline3D::new_array_bc_df_dx(nbin_y, nbin_z);
    double ***tab_bc_df_dy =
        InterCSpline3D::new_array_bc_df_dy(nbin_z, nbin_x);

    for (int iz = 0; iz <= nbin_z; iz++) {
        for (int ix = 0; ix <= nbin_x; ix++) {
            for (int iy = 0; iy <= nbin_y; iy++) {
                double *ptr_df_dz;
                if (iz == 0) {
                    ptr_df_dz = &tab_bc_df_dz[ix][iy][0];
                } else if (iz == nbin_z) {
                    ptr_df_dz = &tab_bc_df_dz[ix][iy][1];
                } else {
                    ptr_df_dz = NULL;
                }

                double *ptr_df_dx;
                if (ix == 0) {
                    ptr_df_dx = &tab_bc_df_dx[iy][iz][0];
                } else if (ix == nbin_x) {
                    ptr_df_dx = &tab_bc_df_dx[iy][iz][1];
                } else {
                    ptr_df_dx = NULL;
                }

                double *ptr_df_dy;
                if (iy == 0) {
                    ptr_df_dy = &tab_bc_df_dy[iz][ix][0];
                } else if (iy == nbin_y) {
                    ptr_df_dy = &tab_bc_df_dy[iz][ix][1];
                } else {
                    ptr_df_dy = NULL;
                }

                tab_f[iz][ix][iy] =
                    func_test(tab_z[iz], tab_x[ix], tab_y[iy],
                              ptr_df_dz,
                              ptr_df_dx,
                              ptr_df_dy);
            }
        }
    }

    // initialize the interpolator
    InterCSpline3D csp_test;
    csp_test.init(nbin_z, nbin_x, nbin_y,
                  tab_z, tab_x, tab_y,
                  tab_f,
                  tab_bc_df_dz,
                  tab_bc_df_dx,
                  tab_bc_df_dy);

    InterCSpline3D::del_array_func(nbin_z, nbin_x, nbin_y,
                                   tab_f);
    InterCSpline3D::del_array_bc_df_dz(nbin_x, nbin_y,
                                       tab_bc_df_dz);
    InterCSpline3D::del_array_bc_df_dx(nbin_y, nbin_z,
                                       tab_bc_df_dx);
    InterCSpline3D::del_array_bc_df_dy(nbin_z, nbin_x,
                                       tab_bc_df_dy);

    delete [] tab_z;
    delete [] tab_x;
    delete [] tab_y;

    double *pt_z = new double[n_pt_z + 1];
    for (int iz = 0; iz <= n_pt_z; iz++) {
        pt_z[iz] = zmin +
            (zmax - zmin) * static_cast<double>(iz) /
                            static_cast<double>(n_pt_z);
    }
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

    double **pt_f_xy_fin =
        InterCSpline2D::new_array_func(n_pt_x, n_pt_y);
    double **pt_f_xy_lin =
        InterCSpline2D::new_array_func(n_pt_x, n_pt_y);
    double **pt_f_yz_fin =
        InterCSpline2D::new_array_func(n_pt_y, n_pt_z);
    double **pt_f_yz_lin =
        InterCSpline2D::new_array_func(n_pt_y, n_pt_z);
    double **pt_f_zx_fin =
        InterCSpline2D::new_array_func(n_pt_z, n_pt_x);
    double **pt_f_zx_lin =
        InterCSpline2D::new_array_func(n_pt_z, n_pt_x);

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
                pt_f_xy_fin[ix][iy] =
                    csp_test.get_func(z_mid, pt_x[ix], pt_y[iy]);
                pt_f_xy_lin[ix][iy] =
                    csp_test.get_func_lin(z_mid, pt_x[ix], pt_y[iy]);
            }
        }

        #ifdef _OPENMP
        // syncronize threads
        #pragma omp barrier
        #endif

        for (int iy = 0; iy <= n_pt_y; iy++) {
            #ifdef _OPENMP
            if (iy % n_thread != tid) {
                continue;
            }
            #endif

            for (int iz = 0; iz <= n_pt_z; iz++) {
                pt_f_yz_fin[iy][iz] =
                    csp_test.get_func(pt_z[iz], x_mid, pt_y[iy]);
                pt_f_yz_lin[iy][iz] =
                    csp_test.get_func_lin(pt_z[iz], x_mid, pt_y[iy]);
            }
        }

        #ifdef _OPENMP
        // syncronize threads
        #pragma omp barrier
        #endif

        for (int iz = 0; iz <= n_pt_z; iz++) {
            #ifdef _OPENMP
            if (iz % n_thread != tid) {
                continue;
            }
            #endif

            for (int ix = 0; ix <= n_pt_x; ix++) {
                pt_f_zx_fin[iz][ix] =
                    csp_test.get_func(pt_z[iz], pt_x[ix], y_mid);
                pt_f_zx_lin[iz][ix] =
                    csp_test.get_func_lin(pt_z[iz], pt_x[ix], y_mid);
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

    std::string axis = "z";
    if (argc > 2) {
        std::string axis_in = argv[2];
        if (axis_in == "z" ||
            axis_in == "Z") {
            axis = "z";
        } else if (axis_in == "x" ||
                   axis_in == "X") {
            axis = "x";
        } else if (axis_in == "y" ||
                   axis_in == "Y") {
            axis = "y";
        } else {
            axis = "z";
        }
    }

    if (axis == "z") {
        for (int ix = 0; ix <= n_pt_x; ix++) {
            for (int iy = 0; iy <= n_pt_y; iy++) {
                double f_ini =
                    func_test(z_mid, pt_x[ix], pt_y[iy]);

                fprintf(ptr_fout, "    %e    %e    %e",
                        pt_x[ix], pt_y[iy], f_ini);
                fprintf(ptr_fout, "    %e    %e",
                        pt_f_xy_fin[ix][iy],
                        pt_f_xy_fin[ix][iy] -
                        f_ini);
                fprintf(ptr_fout, "    %e    %e\n",
                        pt_f_xy_lin[ix][iy],
                        pt_f_xy_lin[ix][iy] -
                        f_ini);
            }
        }
    }
    if (axis == "x") {
        for (int iy = 0; iy <= n_pt_y; iy++) {
            for (int iz = 0; iz <= n_pt_z; iz++) {
                double f_ini =
                    func_test(pt_z[iz], x_mid, pt_y[iy]);

                fprintf(ptr_fout, "    %e    %e    %e",
                        pt_y[iy], pt_z[iz], f_ini);
                fprintf(ptr_fout, "    %e    %e",
                        pt_f_yz_fin[iy][iz],
                        pt_f_yz_fin[iy][iz] -
                        f_ini);
                fprintf(ptr_fout, "    %e    %e\n",
                        pt_f_yz_lin[iy][iz],
                        pt_f_yz_lin[iy][iz] -
                        f_ini);
            }
        }
    }
    if (axis == "y") {
        for (int iz = 0; iz <= n_pt_z; iz++) {
            for (int ix = 0; ix <= n_pt_x; ix++) {
                double f_ini =
                    func_test(pt_z[iz], pt_x[ix], y_mid);

                fprintf(ptr_fout, "    %e    %e    %e",
                        pt_z[iz], pt_x[ix], f_ini);
                fprintf(ptr_fout, "    %e    %e",
                        pt_f_zx_fin[iz][ix],
                        pt_f_zx_fin[iz][ix] -
                        f_ini);
                fprintf(ptr_fout, "    %e    %e\n",
                        pt_f_zx_lin[iz][ix],
                        pt_f_zx_lin[iz][ix] -
                        f_ini);
            }
        }
    }

    if (have_out_ext) {
        fclose(ptr_fout);
    }

    delete [] pt_z;
    delete [] pt_x;
    delete [] pt_y;

    InterCSpline2D::del_array_func(n_pt_x, n_pt_y,
                                   pt_f_xy_fin);
    InterCSpline2D::del_array_func(n_pt_x, n_pt_y,
                                   pt_f_xy_lin);
    InterCSpline2D::del_array_func(n_pt_y, n_pt_z,
                                   pt_f_yz_fin);
    InterCSpline2D::del_array_func(n_pt_y, n_pt_z,
                                   pt_f_yz_lin);
    InterCSpline2D::del_array_func(n_pt_z, n_pt_x,
                                   pt_f_zx_fin);
    InterCSpline2D::del_array_func(n_pt_z, n_pt_x,
                                   pt_f_zx_lin);

    return 0;
}

double func_test(double z, double x, double y,
                 double *df_dz,
                 double *df_dx,
                 double *df_dy) {
    double fac_gauss =
        exp(-0.5 * (x - 0.5) * (x - 0.5) / 0.64) *
        exp(-0.5 * (y - 0.5) * (y - 0.5) / 0.36);
    double f = fac_gauss *
        cos(4. * M_PI * x) * sin(2. * M_PI * y) *
        4. * z * (1. - z);

    if (df_dz != NULL) {
        *df_dz = fac_gauss *
            cos(4. * M_PI * x) * sin(2. * M_PI * y) *
            4. * (1. - 2. * z);
    }

    if (df_dx != NULL) {
        *df_dx = fac_gauss *
            (((0.5 - x) / 0.64) * cos(4. * M_PI * x) -
             4. * M_PI * sin(4. * M_PI * x)) *
            sin(2. * M_PI * y) *
            4. * z * (1. - z);
    }

    if (df_dy != NULL) {
        *df_dy = fac_gauss *
            (((0.5 - y) / 0.36) * sin(2. * M_PI * y) +
             2. * M_PI * cos(2. * M_PI * y)) *
            cos(4. * M_PI * x) *
            4. * z * (1. - z);
    }

    return f;
}
