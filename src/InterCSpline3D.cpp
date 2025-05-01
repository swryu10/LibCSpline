#ifdef _OPENMP
#include<omp.h>
#endif
#include"InterCSpline1D.h"
#include"InterCSpline2D.h"
#include"InterCSpline3D.h"

void InterCSpline3D::init(int nbin_in_z,
                          int nbin_in_x,
                          int nbin_in_y,
                          double *z_in,
                          double *x_in,
                          double *y_in,
                          double ***f_in,
                          double ***bc_df_dz,
                          double ***bc_df_dx,
                          double ***bc_df_dy) {
    reset();

    if (nbin_in_z < 1 ||
        nbin_in_x < 1 ||
        nbin_in_y < 1) {
        return;
    }

    nbin_z_ = nbin_in_z;
    tab_z_ = new double[nbin_z_ + 1];
    for (int iz = 0; iz <= nbin_z_; iz++) {
        tab_z_[iz] = z_in[iz];
    }

    nbin_x_ = nbin_in_x;
    tab_x_ = new double[nbin_x_ + 1];
    for (int ix = 0; ix <= nbin_x_; ix++) {
        tab_x_[ix] = x_in[ix];
    }

    nbin_y_ = nbin_in_y;
    tab_y_ = new double[nbin_y_ + 1];
    for (int iy = 0; iy <= nbin_y_; iy++) {
        tab_y_[iy] = y_in[iy];
    }

    tab_f_ = new_array_func(nbin_z_, nbin_x_, nbin_y_);
    tab_d2f_dz_dz_ = new_array_func(nbin_z_, nbin_x_, nbin_y_);
    tab_d2f_dx_dx_ = new_array_func(nbin_z_, nbin_x_, nbin_y_);
    tab_d2f_dy_dy_ = new_array_func(nbin_z_, nbin_x_, nbin_y_);
    for (int iz = 0; iz <= nbin_z_; iz++) {
        for (int ix = 0; ix <= nbin_x_; ix++) {
            for (int iy = 0; iy <= nbin_y_; iy++) {
                tab_f_[iz][ix][iy] = f_in[iz][ix][iy];
            }
        }
    }

    bool sorted_z = false;
    while (!sorted_z) {
        double **f_l = new double *[nbin_x_ + 1];
        double **f_u = new double *[nbin_x_ + 1];
        for (int ix = 0; ix <= nbin_x_; ix++) {
            f_l[ix] = new double[nbin_y_ + 1];
            f_u[ix] = new double[nbin_y_ + 1];
        }

        bool swapped = false;
        for (int iz = 0; iz < nbin_z_; iz++) {
            if (tab_z_[iz] > tab_z_[iz + 1]) {
                double z_l = tab_z_[iz + 1];
                double z_u = tab_z_[iz];

                for (int ix = 0; ix <= nbin_x_; ix++) {
                    for (int iy = 0; iy <= nbin_y_; iy++) {
                        f_l[ix][iy] = tab_f_[iz + 1][ix][iy];
                        f_u[ix][iy] = tab_f_[iz][ix][iy];
                    }
                }

                tab_z_[iz] = z_l;
                tab_z_[iz + 1] = z_u;

                for (int ix = 0; ix <= nbin_x_; ix++) {
                    for (int iy = 0; iy <= nbin_y_; iy++) {
                        tab_f_[iz][ix][iy] = f_l[ix][iy];
                        tab_f_[iz + 1][ix][iy] = f_u[ix][iy];
                    }
                }

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted_z = true;
        }

        for (int iy = 0; iy <= nbin_y_; iy++) {
            delete [] f_l[iy];
            delete [] f_u[iy];
        }
        delete [] f_l;
        delete [] f_u;
    }

    zmin_ = tab_z_[0];
    zmax_ = tab_z_[nbin_z_];

    bool sorted_x = false;
    while (!sorted_x) {
        double **f_l = new double *[nbin_y_ + 1];
        double **f_u = new double *[nbin_y_ + 1];
        for (int iy = 0; iy <= nbin_y_; iy++) {
            f_l[iy] = new double[nbin_z_ + 1];
            f_u[iy] = new double[nbin_z_ + 1];
        }

        bool swapped = false;
        for (int ix = 0; ix < nbin_x_; ix++) {
            if (tab_x_[ix] > tab_x_[ix + 1]) {
                double x_l = tab_x_[ix + 1];
                double x_u = tab_x_[ix];

                for (int iy = 0; iy <= nbin_y_; iy++) {
                    for (int iz = 0; iz <= nbin_z_; iz++) {
                        f_l[iy][iz] = tab_f_[iz][ix + 1][iy];
                        f_u[iy][iz] = tab_f_[iz][ix][iy];
                    }
                }

                tab_x_[ix] = x_l;
                tab_x_[ix + 1] = x_u;

                for (int iy = 0; iy <= nbin_y_; iy++) {
                    for (int iz = 0; iz <= nbin_z_; iz++) {
                        tab_f_[iz][ix][iy] = f_l[iy][iz];
                        tab_f_[iz][ix + 1][iy] = f_u[iy][iz];
                    }
                }

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted_x = true;
        }

        for (int iy = 0; iy <= nbin_y_; iy++) {
            delete [] f_l[iy];
            delete [] f_u[iy];
        }
        delete [] f_l;
        delete [] f_u;
    }

    xmin_ = tab_x_[0];
    xmax_ = tab_x_[nbin_x_];

    bool sorted_y = false;
    while (!sorted_y) {
        double **f_l = new double *[nbin_z_ + 1];
        double **f_u = new double *[nbin_z_ + 1];
        for (int iz = 0; iz <= nbin_z_; iz++) {
            f_l[iz] = new double[nbin_x_ + 1];
            f_u[iz] = new double[nbin_x_ + 1];
        }

        bool swapped = false;
        for (int iy = 0; iy < nbin_y_; iy++) {
            if (tab_y_[iy] > tab_y_[iy + 1]) {
                double y_l = tab_y_[iy + 1];
                double y_u = tab_y_[iy];

                for (int iz = 0; iz <= nbin_z_; iz++) {
                    for (int ix = 0; ix <= nbin_x_; ix++) {
                        f_l[iz][ix] = tab_f_[iz][ix][iy + 1];
                        f_u[iz][ix] = tab_f_[iz][ix][iy];
                    }
                }

                tab_y_[iy] = y_l;
                tab_y_[iy + 1] = y_u;

                for (int iz = 0; iz <= nbin_z_; iz++) {
                    for (int ix = 0; ix <= nbin_x_; ix++) {
                        tab_f_[iz][ix][iy] = f_l[iz][ix];
                        tab_f_[iz][ix][iy + 1] = f_u[iz][ix];
                    }
                }

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted_y = true;
        }

        for (int iz = 0; iz <= nbin_z_; iz++) {
            delete [] f_l[iz];
            delete [] f_u[iz];
        }
        delete [] f_l;
        delete [] f_u;
    }

    ymin_ = tab_y_[0];
    ymax_ = tab_y_[nbin_y_];

    InterCSpline2D *list_csp_xy_f =
        new InterCSpline2D [nbin_z_ + 1]();

    InterCSpline1D **list_csp_z_f =
        new InterCSpline1D *[nbin_x_ + 1];
    for (int ix = 0; ix <= nbin_x_; ix++) {
        list_csp_z_f[ix] = new InterCSpline1D [nbin_y_ + 1]();
    }

    for (int iz = 0; iz <= nbin_z_; iz++) {
        double **ptr_bc_df_dx = NULL;
        if (bc_df_dx != NULL) {
            ptr_bc_df_dx = new double *[nbin_y_ + 1];
            for (int iy = 0; iy <= nbin_y_; iy++) {
                ptr_bc_df_dx[iy] = new double[2];
                ptr_bc_df_dx[iy][0] = bc_df_dx[iy][iz][0];
                ptr_bc_df_dx[iy][1] = bc_df_dx[iy][iz][1];
            }
        }

        double **ptr_bc_df_dy = NULL;
        if (bc_df_dy != NULL) {
            ptr_bc_df_dy = new double *[nbin_x_ + 1];
            for (int ix = 0; ix <= nbin_x_; ix++) {
                ptr_bc_df_dy[ix] = new double[2];
                ptr_bc_df_dy[ix][0] = bc_df_dy[iz][ix][0];
                ptr_bc_df_dy[ix][1] = bc_df_dy[iz][ix][1];
            }
        }

        list_csp_xy_f[iz].init(nbin_x_, nbin_y_,
                               tab_x_, tab_y_,
                               tab_f_[iz],
                               ptr_bc_df_dx,
                               ptr_bc_df_dy);

        if (bc_df_dx != NULL) {
            for (int iy = 0; iy <= nbin_y_; iy++) {
                delete [] ptr_bc_df_dx[iy];
            }
            delete [] ptr_bc_df_dx;
        }

        if (bc_df_dy != NULL) {
            for (int ix = 0; ix <= nbin_x_; ix++) {
                delete [] ptr_bc_df_dy[ix];
            }
            delete [] ptr_bc_df_dy;
        }

        for (int ix = 0; ix <= nbin_x_; ix++) {
            for (int iy = 0; iy <= nbin_y_; iy++) {
                double *ptr_df_dx = NULL;
                double *ptr_df_dy = NULL;
                double *ptr_d2f_dx_dx =
                    &tab_d2f_dx_dx_[iz][ix][iy];
                double *ptr_d2f_dx_dy =
                    NULL;
                double *ptr_d2f_dy_dy =
                    &tab_d2f_dy_dy_[iz][ix][iy];
                double f_tmp =
                    list_csp_xy_f[iz].get_func(tab_x_[ix],
                                               tab_y_[iy],
                                               ptr_df_dx,
                                               ptr_df_dy,
                                               ptr_d2f_dx_dx,
                                               ptr_d2f_dx_dy,
                                               ptr_d2f_dy_dy);
            }
        }
    }

    #ifdef _OPENMP
    #pragma omp parallel
    {  // parallel code begins
    #endif
        #ifdef _OPENMP
        int n_thread = omp_get_num_threads();
        int tid = omp_get_thread_num();
        #endif

        for (int ix = 0; ix <= nbin_x_; ix++) {
            #ifdef _OPENMP
            if (ix % n_thread != tid) {
                continue;
            }
            #endif

            for (int iy = 0; iy <= nbin_y_; iy++) {
                double *tab_tmp_f = new double[nbin_z_ + 1];

                for (int iz = 0; iz <= nbin_z_; iz++) {
                    tab_tmp_f[iz] = tab_f_[iz][ix][iy];
                }

                double *ptr_bc_df_dz = NULL;
                if (bc_df_dz != NULL) {
                    ptr_bc_df_dz = bc_df_dz[ix][iy];
                }
                list_csp_z_f[ix][iy].init(nbin_z_,
                                          tab_z_,
                                          tab_tmp_f,
                                          ptr_bc_df_dz);

                for (int iz = 0; iz <= nbin_z_; iz++) {
                    double *ptr_df_dz = NULL;
                    double *ptr_d2f_dz_dz =
                        &tab_d2f_dz_dz_[iz][ix][iy];
                    double f_tmp =
                        list_csp_z_f[ix][iy].get_func(tab_z_[iz],
                                                      ptr_df_dz,
                                                      ptr_d2f_dz_dz);
                }

                delete [] tab_tmp_f;
            }
        }
    #ifdef _OPENMP
    }  // parallel code ends
    #endif

    for (int ix = 0; ix <= nbin_x_; ix++) {
        delete [] list_csp_z_f[ix];
    }
    delete [] list_csp_z_f;

    delete [] list_csp_xy_f;

    initialized_ = true;

    return;
}

double InterCSpline3D::get_func(double z_in,
                                double x_in,
                                double y_in,
                                double *ptr_df_dz_out,
                                double *ptr_df_dx_out,
                                double *ptr_df_dy_out,
                                double *ptr_d2f_dz_dz_out,
                                double *ptr_d2f_dz_dx_out,
                                double *ptr_d2f_dx_dx_out,
                                double *ptr_d2f_dx_dy_out,
                                double *ptr_d2f_dy_dy_out,
                                double *ptr_d2f_dy_dz_out) {
    if (!initialized_) {
        if (ptr_df_dz_out != NULL) {
            *ptr_df_dz_out = 0.;
        }

        if (ptr_df_dx_out != NULL) {
            *ptr_df_dx_out = 0.;
        }

        if (ptr_df_dy_out != NULL) {
            *ptr_df_dy_out = 0.;
        }

        if (ptr_d2f_dz_dz_out != NULL) {
            *ptr_d2f_dz_dz_out = 0.;
        }

        if (ptr_d2f_dz_dx_out != NULL) {
            *ptr_d2f_dz_dx_out = 0.;
        }

        if (ptr_d2f_dx_dx_out != NULL) {
            *ptr_d2f_dx_dx_out = 0.;
        }

        if (ptr_d2f_dx_dy_out != NULL) {
            *ptr_d2f_dx_dy_out = 0.;
        }

        if (ptr_d2f_dy_dy_out != NULL) {
            *ptr_d2f_dy_dy_out = 0.;
        }

        if (ptr_d2f_dy_dz_out != NULL) {
            *ptr_d2f_dy_dz_out = 0.;
        }

        return 0.;
    }

    int iz = get_index_z(z_in);
    int ix = get_index_x(x_in);
    int iy = get_index_y(y_in);

    double delta_z = tab_z_[iz + 1] - tab_z_[iz];
    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_y = tab_y_[iy + 1] - tab_y_[iy];

    double *frac_z_d0f = new double [2];
    frac_z_d0f[0] = (tab_z_[iz + 1] - z_in) / delta_z;
    frac_z_d0f[1] = 1. - frac_z_d0f[0];

    double **coeff_z_f = new double *[3];
    coeff_z_f[0] = new double[2];
    coeff_z_f[1] = new double[2];
    coeff_z_f[2] = new double[2];
    coeff_z_f[0][0] = frac_z_d0f[0];
    coeff_z_f[0][1] = frac_z_d0f[1];
    coeff_z_f[1][0] =
        -(3. * frac_z_d0f[0] * frac_z_d0f[0] - 1.) / 6.;
    coeff_z_f[1][1] =
        (3. * frac_z_d0f[1] * frac_z_d0f[1] - 1.) / 6.;
    coeff_z_f[2][0] =
        frac_z_d0f[0] * (frac_z_d0f[0] * frac_z_d0f[0] - 1.) / 6.;
    coeff_z_f[2][1] =
        frac_z_d0f[1] * (frac_z_d0f[1] * frac_z_d0f[1] - 1.) / 6.;

    double *frac_x_d0f = new double [2];
    frac_x_d0f[0] = (tab_x_[ix + 1] - x_in) / delta_x;
    frac_x_d0f[1] = 1. - frac_x_d0f[0];

    double **coeff_x_f = new double *[3];
    coeff_x_f[0] = new double[2];
    coeff_x_f[1] = new double[2];
    coeff_x_f[2] = new double[2];
    coeff_x_f[0][0] = frac_x_d0f[0];
    coeff_x_f[0][1] = frac_x_d0f[1];
    coeff_x_f[1][0] =
        -(3. * frac_x_d0f[0] * frac_x_d0f[0] - 1.) / 6.;
    coeff_x_f[1][1] =
        (3. * frac_x_d0f[1] * frac_x_d0f[1] - 1.) / 6.;
    coeff_x_f[2][0] =
        frac_x_d0f[0] * (frac_x_d0f[0] * frac_x_d0f[0] - 1.) / 6.;
    coeff_x_f[2][1] =
        frac_x_d0f[1] * (frac_x_d0f[1] * frac_x_d0f[1] - 1.) / 6.;

    double *frac_y_d0f = new double [2];
    frac_y_d0f[0] = (tab_y_[iy + 1] - y_in) / delta_y;
    frac_y_d0f[1] = 1. - frac_y_d0f[0];

    double **coeff_y_f = new double *[3];
    coeff_y_f[0] = new double[2];
    coeff_y_f[1] = new double[2];
    coeff_y_f[2] = new double[2];
    coeff_y_f[0][0] = frac_y_d0f[0];
    coeff_y_f[0][1] = frac_y_d0f[1];
    coeff_y_f[1][0] =
        -(3. * frac_y_d0f[0] * frac_y_d0f[0] - 1.) / 6.;
    coeff_y_f[1][1] =
        (3. * frac_y_d0f[1] * frac_y_d0f[1] - 1.) / 6.;
    coeff_y_f[2][0] =
        frac_y_d0f[0] * (frac_y_d0f[0] * frac_y_d0f[0] - 1.) / 6.;
    coeff_y_f[2][1] =
        frac_y_d0f[1] * (frac_y_d0f[1] * frac_y_d0f[1] - 1.) / 6.;

    double f_out = 0.;
    for (int k = 0; k < 8; k++) {
        int kz = k % 2;
        int kxy = (k - kz) / 2;
        int kx = kxy % 2;
        int ky = (kxy - kx) / 2;

        f_out +=
            coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
            tab_f_[iz + kz][ix + kx][iy + ky] +
            delta_z * delta_z *
            coeff_z_f[2][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
            tab_d2f_dz_dz_[iz + kz][ix + kx][iy + ky] +
            delta_x * delta_x *
            coeff_z_f[0][kz] * coeff_x_f[2][kx] * coeff_y_f[0][ky] *
            tab_d2f_dx_dx_[iz + kz][ix + kx][iy + ky] +
            delta_y * delta_y *
            coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[2][ky] *
            tab_d2f_dy_dy_[iz + kz][ix + kx][iy + ky];
    }

    if (ptr_df_dz_out != NULL) {
        *ptr_df_dz_out = 0.;
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_df_dz_out +=
                coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                    (tab_f_[iz + 1][ix + kx][iy + ky] -
                     tab_f_[iz][ix + kx][iy + ky]) / delta_z +
                coeff_x_f[2][kx] * coeff_y_f[0][ky] *
                    (tab_d2f_dx_dx_[iz + 1][ix + kx][iy + ky] -
                     tab_d2f_dx_dx_[iz][ix + kx][iy + ky]) *
                    delta_x * delta_x / delta_z +
                coeff_x_f[0][kx] * coeff_y_f[2][ky] *
                    (tab_d2f_dy_dy_[iz + 1][ix + kx][iy + ky] -
                     tab_d2f_dy_dy_[iz][ix + kx][iy + ky]) *
                    delta_y * delta_y / delta_z;
        }
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_df_dz_out +=
                delta_z *
                coeff_z_f[1][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dz_dz_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_df_dx_out != NULL) {
        *ptr_df_dx_out = 0.;
        for (int k = 0; k < 4; k++) {
            int ky = k % 2;
            int kz = (k - ky) / 2;

            *ptr_df_dx_out +=
                coeff_y_f[0][ky] * coeff_z_f[0][kz] *
                    (tab_f_[iz + kz][ix + 1][iy + ky] -
                     tab_f_[iz + kz][ix][iy + ky]) / delta_x +
                coeff_y_f[2][ky] * coeff_z_f[0][kz] *
                    (tab_d2f_dy_dy_[iz + kz][ix + 1][iy + ky] -
                     tab_d2f_dy_dy_[iz + kz][ix][iy + ky]) *
                    delta_y * delta_y / delta_x +
                coeff_y_f[0][ky] * coeff_z_f[2][kz] *
                    (tab_d2f_dz_dz_[iz + kz][ix + 1][iy + ky] -
                     tab_d2f_dz_dz_[iz + kz][ix][iy + ky]) *
                    delta_z * delta_z / delta_x;
        }
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_df_dx_out +=
                delta_x *
                coeff_z_f[0][kz] * coeff_x_f[1][kx] * coeff_y_f[0][ky] *
                tab_d2f_dx_dx_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_df_dy_out != NULL) {
        *ptr_df_dy_out = 0.;
        for (int k = 0; k < 4; k++) {
            int kz = k % 2;
            int kx = (k - kz) / 2;

            *ptr_df_dy_out +=
                coeff_z_f[0][kz] * coeff_x_f[0][kx] *
                    (tab_f_[iz + kz][ix + kx][iy + 1] -
                     tab_f_[iz + kz][ix + kx][iy]) / delta_y +
                coeff_z_f[2][kz] * coeff_x_f[0][kx] *
                    (tab_d2f_dz_dz_[iz + kz][ix + kx][iy + 1] -
                     tab_d2f_dz_dz_[iz + kz][ix + kx][iy]) *
                    delta_z * delta_z / delta_y +
                coeff_z_f[0][kz] * coeff_x_f[2][kx] *
                    (tab_d2f_dx_dx_[iz + kz][ix + kx][iy + 1] -
                     tab_d2f_dx_dx_[iz + kz][ix + kx][iy]) *
                    delta_x * delta_x / delta_y;
        }
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_df_dy_out +=
                delta_y *
                coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[1][ky] *
                tab_d2f_dy_dy_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dz_dz_out != NULL) {
        *ptr_d2f_dz_dz_out = 0.;
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_d2f_dz_dz_out +=
                coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dz_dz_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dz_dx_out != NULL) {
        *ptr_d2f_dz_dx_out = 0.;
        for (int ky = 0; ky < 2; ky++) {
            *ptr_d2f_dz_dx_out +=
                coeff_y_f[0][ky] *
                    (tab_f_[iz + 1][ix + 1][iy + ky] -
                     tab_f_[iz][ix + 1][iy + ky] -
                     tab_f_[iz + 1][ix][iy + ky] +
                     tab_f_[iz][ix][iy + ky]) /
                    (delta_z * delta_x) +
                coeff_y_f[2][ky] *
                    (tab_d2f_dy_dy_[iz + 1][ix + 1][iy + ky] -
                     tab_d2f_dy_dy_[iz][ix + 1][iy + ky] -
                     tab_d2f_dy_dy_[iz + 1][ix][iy + ky] +
                     tab_d2f_dy_dy_[iz][ix][iy + ky]) *
                    delta_y * delta_y / (delta_z * delta_x);
        }
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_d2f_dz_dx_out +=
                coeff_x_f[1][kx] * coeff_y_f[0][ky] *
                    (tab_d2f_dx_dx_[iz + 1][ix + kx][iy + ky] -
                     tab_d2f_dx_dx_[iz][ix + kx][iy + ky]) *
                    delta_x / delta_z;
        }
        for (int k = 0; k < 4; k++) {
            int ky = k % 2;
            int kz = (k - ky) / 2;

            *ptr_d2f_dz_dx_out +=
                coeff_z_f[1][kz] * coeff_y_f[0][ky] *
                    (tab_d2f_dz_dz_[iz + kz][ix + 1][iy + ky] -
                     tab_d2f_dz_dz_[iz + kz][ix][iy + ky]) *
                    delta_z / delta_x;
        }
    }

    if (ptr_d2f_dx_dx_out != NULL) {
        *ptr_d2f_dx_dx_out = 0.;
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_d2f_dx_dx_out +=
                coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dx_dx_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dx_dy_out != NULL) {
        *ptr_d2f_dx_dy_out = 0.;
        for (int kz = 0; kz < 2; kz++) {
            *ptr_d2f_dx_dy_out +=
                coeff_z_f[0][kz] *
                    (tab_f_[iz + kz][ix + 1][iy + 1] -
                     tab_f_[iz + kz][ix + 1][iy] -
                     tab_f_[iz + kz][ix][iy + 1] +
                     tab_f_[iz + kz][ix][iy]) /
                    (delta_x * delta_y) +
                coeff_z_f[2][kz] *
                    (tab_d2f_dz_dz_[iz + kz][ix + 1][iy + 1] -
                     tab_d2f_dz_dz_[iz + kz][ix + 1][iy] -
                     tab_d2f_dz_dz_[iz + kz][ix][iy + 1] +
                     tab_d2f_dz_dz_[iz + kz][ix][iy]) *
                    delta_z * delta_z / (delta_x * delta_y);
        }
        for (int k = 0; k < 4; k++) {
            int ky = k % 2;
            int kz = (k - ky) / 2;

            *ptr_d2f_dx_dy_out +=
                coeff_y_f[1][ky] * coeff_z_f[0][kz] *
                    (tab_d2f_dy_dy_[iz + kz][ix + 1][iy + ky] -
                     tab_d2f_dy_dy_[iz + kz][ix][iy + ky]) *
                    delta_y / delta_x;
        }
        for (int k = 0; k < 4; k++) {
            int kz = k % 2;
            int kx = (k - kz) / 2;

            *ptr_d2f_dx_dy_out +=
                coeff_x_f[1][kx] * coeff_z_f[0][kz] *
                    (tab_d2f_dx_dx_[iz + kz][ix + kx][iy + 1] -
                     tab_d2f_dx_dx_[iz + kz][ix + kx][iy]) *
                    delta_x / delta_y;
        }
    }

    if (ptr_d2f_dy_dy_out != NULL) {
        *ptr_d2f_dy_dy_out = 0.;
        for (int k = 0; k < 8; k++) {
            int kz = k % 2;
            int kxy = (k - kz) / 2;
            int kx = kxy % 2;
            int ky = (kxy - kx) / 2;

            *ptr_d2f_dy_dy_out +=
                coeff_z_f[0][kz] * coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dy_dy_[iz + kz][ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dy_dz_out != NULL) {
        *ptr_d2f_dy_dz_out = 0.;
        for (int kx = 0; kx < 2; kx++) {
            *ptr_d2f_dy_dz_out +=
                coeff_x_f[0][kx] *
                    (tab_f_[iz + 1][ix + kx][iy + 1] -
                     tab_f_[iz + 1][ix + kx][iy] -
                     tab_f_[iz][ix + kx][iy + 1] +
                     tab_f_[iz][ix + kx][iy]) /
                    (delta_y * delta_z) +
                coeff_x_f[2][kx] *
                    (tab_d2f_dx_dx_[iz + 1][ix + kx][iy + 1] -
                     tab_d2f_dx_dx_[iz + 1][ix + kx][iy] -
                     tab_d2f_dx_dx_[iz][ix + kx][iy + 1] +
                     tab_d2f_dx_dx_[iz][ix + kx][iy]) *
                    delta_x * delta_x / (delta_y * delta_z);
        }
        for (int k = 0; k < 4; k++) {
            int kz = k % 2;
            int kx = (k - kz) / 2;

            *ptr_d2f_dy_dz_out +=
                coeff_z_f[1][kz] * coeff_x_f[0][kx] *
                    (tab_d2f_dz_dz_[iz + kz][ix + kx][iy + 1] -
                     tab_d2f_dz_dz_[iz + kz][ix + kx][iy]) *
                    delta_z / delta_y;
        }
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_d2f_dy_dz_out +=
                coeff_y_f[1][ky] * coeff_x_f[0][kx] *
                    (tab_d2f_dy_dy_[iz + 1][ix + kx][iy + ky] -
                     tab_d2f_dy_dy_[iz][ix + kx][iy + ky]) *
                    delta_y / delta_z;
        }
    }

    delete [] coeff_z_f[0];
    delete [] coeff_z_f[1];
    delete [] coeff_z_f[2];
    delete [] coeff_z_f;

    delete [] coeff_x_f[0];
    delete [] coeff_x_f[1];
    delete [] coeff_x_f[2];
    delete [] coeff_x_f;

    delete [] coeff_y_f[0];
    delete [] coeff_y_f[1];
    delete [] coeff_y_f[2];
    delete [] coeff_y_f;

    delete [] frac_z_d0f;
    delete [] frac_x_d0f;
    delete [] frac_y_d0f;

    return f_out;
}

double InterCSpline3D::get_func_lin(double z_in,
                                    double x_in,
                                    double y_in) {
    if (!initialized_) {
        return 0.;
    }

    int iz = get_index_z(z_in);
    int ix = get_index_x(x_in);
    int iy = get_index_y(y_in);

    double delta_z = tab_z_[iz + 1] - tab_z_[iz];
    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_y = tab_y_[iy + 1] - tab_y_[iy];

    double *frac_z_d0f = new double [2];
    frac_z_d0f[0] = (tab_z_[iz + 1] - z_in) / delta_z;
    frac_z_d0f[1] = 1. - frac_z_d0f[0];

    double *frac_x_d0f = new double [2];
    frac_x_d0f[0] = (tab_x_[ix + 1] - x_in) / delta_x;
    frac_x_d0f[1] = 1. - frac_x_d0f[0];

    double *frac_y_d0f = new double [2];
    frac_y_d0f[0] = (tab_y_[iy + 1] - y_in) / delta_y;
    frac_y_d0f[1] = 1. - frac_y_d0f[0];

    double f_out = 0.;
    for (int k = 0; k < 8; k++) {
        int kz = k % 2;
        int kxy = (k - kz) / 2;
        int kx = kxy % 2;
        int ky = (kxy - kx) / 2;

        f_out +=
            frac_z_d0f[kz] * frac_x_d0f[kx] * frac_y_d0f[ky] *
            tab_f_[iz + kz][ix + kx][iy + ky];
    }

    delete [] frac_z_d0f;
    delete [] frac_x_d0f;
    delete [] frac_y_d0f;

    return f_out;
}

int InterCSpline3D::get_index_z(double z_in) {
    int iz = 0;
    if (z_in < zmin_) {
        iz = 0;
    } else if (z_in >= zmax_) {
        iz = nbin_z_ - 1;
    } else {
        for (int jz = 0; jz < nbin_z_; jz++) {
            if (z_in >= tab_z_[jz] && z_in < tab_z_[jz + 1]) {
                iz = jz;
                break;
            }
        }
    }

    return iz;
}

int InterCSpline3D::get_index_x(double x_in) {
    int ix = 0;
    if (x_in < xmin_) {
        ix = 0;
    } else if (x_in >= xmax_) {
        ix = nbin_x_ - 1;
    } else {
        for (int jx = 0; jx < nbin_x_; jx++) {
            if (x_in >= tab_x_[jx] && x_in < tab_x_[jx + 1]) {
                ix = jx;
                break;
            }
        }
    }

    return ix;
}

int InterCSpline3D::get_index_y(double y_in) {
    int iy = 0;
    if (y_in < ymin_) {
        iy = 0;
    } else if (y_in >= ymax_) {
        iy = nbin_y_ - 1;
    } else {
        for (int jy = 0; jy < nbin_y_; jy++) {
            if (y_in >= tab_y_[jy] && y_in < tab_y_[jy + 1]) {
                iy = jy;
                break;
            }
        }
    }

    return iy;
}
