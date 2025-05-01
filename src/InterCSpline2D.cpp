#ifdef _OPENMP
#include<omp.h>
#endif
#include"InterCSpline1D.h"
#include"InterCSpline2D.h"

void InterCSpline2D::init(int nbin_in_x,
                          int nbin_in_y,
                          double *x_in,
                          double *y_in,
                          double **f_in,
                          double **bc_df_dx,
                          double **bc_df_dy) {
    reset();

    if (nbin_in_x < 1 ||
        nbin_in_y < 1) {
        return;
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

    tab_f_ = new double *[nbin_x_ + 1];
    tab_d2f_dx_dx_ = new double *[nbin_x_ + 1];
    tab_d2f_dy_dy_ = new double *[nbin_x_ + 1];

    for (int ix = 0; ix <= nbin_x_; ix++) {
        tab_f_[ix] = new double[nbin_y_ + 1];
        tab_d2f_dx_dx_[ix] = new double[nbin_y_ + 1];
        tab_d2f_dy_dy_[ix] = new double[nbin_y_ + 1];

        for (int iy = 0; iy <= nbin_y_; iy++) {
            tab_f_[ix][iy] = f_in[ix][iy];

            tab_d2f_dx_dx_[ix][iy] = 0.;
            tab_d2f_dy_dy_[ix][iy] = 0.;
        }
    }

    double **tab_f_tr = new double *[nbin_y_ + 1];
    for (int iy = 0; iy <= nbin_y_; iy++) {
        tab_f_tr[iy] = new double[nbin_x_ + 1];

        for (int ix = 0; ix <= nbin_x_; ix++) {
            tab_f_tr[iy][ix] = tab_f_[ix][iy];
        }
    }

    bool sorted_x = false;
    while (!sorted_x) {
        double *f_l = new double[nbin_y_ + 1];
        double *f_u = new double[nbin_y_ + 1];

        bool swapped = false;
        for (int ix = 0; ix < nbin_x_; ix++) {
            if (tab_x_[ix] > tab_x_[ix + 1]) {
                double x_l = tab_x_[ix + 1];
                double x_u = tab_x_[ix];

                for (int iy = 0; iy <= nbin_y_; iy++) {
                    f_l[iy] = tab_f_[ix + 1][iy];
                    f_u[iy] = tab_f_[ix][iy];
                }

                tab_x_[ix] = x_l;
                tab_x_[ix + 1] = x_u;

                for (int iy = 0; iy <= nbin_y_; iy++) {
                    tab_f_[ix][iy] = f_l[iy];
                    tab_f_[ix + 1][iy] = f_u[iy];
                }

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted_x = true;
        }

        delete [] f_l;
        delete [] f_u;
    }

    xmin_ = tab_x_[0];
    xmax_ = tab_x_[nbin_x_];

    bool sorted_y = false;
    while (!sorted_y) {
        double *f_l = new double[nbin_x_ + 1];
        double *f_u = new double[nbin_x_ + 1];

        bool swapped = false;
        for (int iy = 0; iy < nbin_y_; iy++) {
            if (tab_y_[iy] > tab_y_[iy + 1]) {
                double y_l = tab_y_[iy + 1];
                double y_u = tab_y_[iy];

                for (int ix = 0; ix <= nbin_x_; ix++) {
                    f_l[ix] = tab_f_[ix][iy + 1];
                    f_u[ix] = tab_f_[ix][iy];
                }

                tab_y_[iy] = y_l;
                tab_y_[iy + 1] = y_u;

                for (int ix = 0; ix <= nbin_x_; ix++) {
                    tab_f_[ix][iy] = f_l[ix];
                    tab_f_[ix][iy + 1] = f_u[ix];
                }

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted_y = true;
        }

        delete [] f_l;
        delete [] f_u;
    }

    ymin_ = tab_y_[0];
    ymax_ = tab_y_[nbin_y_];

    InterCSpline1D *list_csp_y_f =
        new InterCSpline1D [nbin_x_ + 1]();
    InterCSpline1D *list_csp_x_f =
        new InterCSpline1D [nbin_y_ + 1]();

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

            double *ptr_bc_df_dy = NULL;
            if (bc_df_dy != NULL) {
                ptr_bc_df_dy = bc_df_dy[ix];
            }
            list_csp_y_f[ix].init(nbin_y_,
                                  tab_y_,
                                  tab_f_[ix],
                                  ptr_bc_df_dy);

            for (int iy = 0; iy <= nbin_y_; iy++) {
                double *ptr_df_dy = NULL;
                double *ptr_d2f_dy_dy =
                    &tab_d2f_dy_dy_[ix][iy];
                double f_tmp =
                    list_csp_y_f[ix].get_func(tab_y_[iy],
                                              ptr_df_dy,
                                              ptr_d2f_dy_dy);
            }
        }

        #ifdef _OPENMP
        // syncronize threads
        #pragma omp barrier
        #endif

        for (int iy = 0; iy <= nbin_y_; iy++) {
            #ifdef _OPENMP
            if (iy % n_thread != tid) {
                continue;
            }
            #endif

            double *ptr_bc_df_dx = NULL;
            if (bc_df_dx != NULL) {
                ptr_bc_df_dx = bc_df_dx[iy];
            }
            list_csp_x_f[iy].init(nbin_x_,
                                  tab_x_,
                                  tab_f_tr[iy],
                                  ptr_bc_df_dx);

            for (int ix = 0; ix <= nbin_x_; ix++) {
                double *ptr_df_dx = NULL;
                double *ptr_d2f_dx_dx =
                    &tab_d2f_dx_dx_[ix][iy];
                double f_tmp =
                    list_csp_x_f[iy].get_func(tab_x_[ix],
                                              ptr_df_dx,
                                              ptr_d2f_dx_dx);
            }
        }
    #ifdef _OPENMP
    }  // parallel code ends
    #endif

    delete [] list_csp_y_f;
    delete [] list_csp_x_f;

    for (int iy = 0; iy <= nbin_y_; iy++) {
        delete [] tab_f_tr[iy];
    }
    delete [] tab_f_tr;

    initialized_ = true;

    return;
}

double InterCSpline2D::get_func(double x_in,
                                double y_in,
                                double *ptr_df_dx_out,
                                double *ptr_df_dy_out,
                                double *ptr_d2f_dx_dx_out,
                                double *ptr_d2f_dx_dy_out,
                                double *ptr_d2f_dy_dy_out) {
    if (!initialized_) {
        if (ptr_df_dx_out != NULL) {
            *ptr_df_dx_out = 0.;
        }

        if (ptr_df_dy_out != NULL) {
            *ptr_df_dy_out = 0.;
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

        return 0.;
    }

    int ix = get_index_x(x_in);
    int iy = get_index_y(y_in);

    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_y = tab_y_[iy + 1] - tab_y_[iy];

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
    for (int k = 0; k < 4; k++) {
        int kx = k % 2;
        int ky = (k - kx) / 2;

        f_out +=
            coeff_x_f[0][kx] * coeff_y_f[0][ky] *
            tab_f_[ix + kx][iy + ky] +
            delta_x * delta_x *
            coeff_x_f[2][kx] * coeff_y_f[0][ky] *
            tab_d2f_dx_dx_[ix + kx][iy + ky] +
            delta_y * delta_y *
            coeff_x_f[0][kx] * coeff_y_f[2][ky] *
            tab_d2f_dy_dy_[ix + kx][iy + ky];
    }

    if (ptr_df_dx_out != NULL) {
        *ptr_df_dx_out =
            coeff_y_f[0][0] *
            (tab_f_[ix + 1][iy] - tab_f_[ix][iy]) / delta_x +
            coeff_y_f[0][1] *
            (tab_f_[ix + 1][iy + 1] - tab_f_[ix][iy + 1]) / delta_x +
            coeff_y_f[2][0] * (delta_y * delta_y / delta_x) *
            (tab_d2f_dy_dy_[ix + 1][iy] - tab_d2f_dy_dy_[ix][iy]) +
            coeff_y_f[2][1] * (delta_y * delta_y / delta_x) *
            (tab_d2f_dy_dy_[ix + 1][iy + 1] - tab_d2f_dy_dy_[ix][iy + 1]);
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_df_dx_out +=
                delta_x *
                coeff_x_f[1][kx] * coeff_y_f[0][ky] *
                tab_d2f_dx_dx_[ix + kx][iy + ky];
        }
    }

    if (ptr_df_dy_out != NULL) {
        *ptr_df_dy_out =
            coeff_x_f[0][0] *
            (tab_f_[ix][iy + 1] - tab_f_[ix][iy]) / delta_y +
            coeff_x_f[0][1] *
            (tab_f_[ix + 1][iy + 1] - tab_f_[ix + 1][iy]) / delta_y +
            coeff_x_f[2][0] * (delta_x * delta_x / delta_y) *
            (tab_d2f_dx_dx_[ix][iy + 1] - tab_d2f_dx_dx_[ix][iy]) +
            coeff_x_f[2][1] * (delta_x * delta_x / delta_y) *
            (tab_d2f_dx_dx_[ix + 1][iy + 1] - tab_d2f_dx_dx_[ix + 1][iy]);

        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_df_dy_out +=
                delta_y *
                coeff_x_f[0][kx] * coeff_y_f[1][ky] *
                tab_d2f_dy_dy_[ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dx_dx_out != NULL) {
        *ptr_d2f_dx_dx_out = 0.;
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_d2f_dx_dx_out +=
                coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dx_dx_[ix + kx][iy + ky];
        }
    }

    if (ptr_d2f_dx_dy_out != NULL) {
        *ptr_d2f_dx_dy_out =
            (tab_f_[ix + 1][iy + 1] - tab_f_[ix][iy + 1] -
             tab_f_[ix + 1][iy] + tab_f_[ix][iy]) /
                (delta_x * delta_y) +
            coeff_x_f[1][0] * (delta_x / delta_y) *
                (tab_d2f_dx_dx_[ix][iy + 1] -
                 tab_d2f_dx_dx_[ix][iy]) +
            coeff_x_f[1][1] * (delta_x / delta_y) *
                (tab_d2f_dx_dx_[ix + 1][iy + 1] -
                 tab_d2f_dx_dx_[ix + 1][iy]) +
            coeff_y_f[1][0] * (delta_y / delta_x) *
                (tab_d2f_dy_dy_[ix + 1][iy] -
                 tab_d2f_dy_dy_[ix][iy]) +
            coeff_y_f[1][1] * (delta_y / delta_x) *
                (tab_d2f_dy_dy_[ix + 1][iy + 1] -
                 tab_d2f_dy_dy_[ix][iy + 1]);            
    }

    if (ptr_d2f_dy_dy_out != NULL) {
        *ptr_d2f_dy_dy_out = 0.;
        for (int k = 0; k < 4; k++) {
            int kx = k % 2;
            int ky = (k - kx) / 2;

            *ptr_d2f_dy_dy_out +=
                coeff_x_f[0][kx] * coeff_y_f[0][ky] *
                tab_d2f_dy_dy_[ix + kx][iy + ky];
        }
    }

    delete [] coeff_x_f[0];
    delete [] coeff_x_f[1];
    delete [] coeff_x_f[2];
    delete [] coeff_x_f;

    delete [] coeff_y_f[0];
    delete [] coeff_y_f[1];
    delete [] coeff_y_f[2];
    delete [] coeff_y_f;

    delete [] frac_x_d0f;
    delete [] frac_y_d0f;

    return f_out;
}

double InterCSpline2D::get_func_lin(double x_in,
                                    double y_in) {
    if (!initialized_) {
        return 0.;
    }

    int ix = get_index_x(x_in);
    int iy = get_index_y(y_in);

    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_y = tab_y_[iy + 1] - tab_y_[iy];

    double *frac_x_d0f = new double [2];
    frac_x_d0f[0] = (tab_x_[ix + 1] - x_in) / delta_x;
    frac_x_d0f[1] = 1. - frac_x_d0f[0];

    double *frac_y_d0f = new double [2];
    frac_y_d0f[0] = (tab_y_[iy + 1] - y_in) / delta_y;
    frac_y_d0f[1] = 1. - frac_y_d0f[0];

    double f_out =
        frac_x_d0f[0] * frac_y_d0f[0] * tab_f_[ix][iy] +
        frac_x_d0f[0] * frac_y_d0f[1] * tab_f_[ix][iy + 1] +
        frac_x_d0f[1] * frac_y_d0f[0] * tab_f_[ix + 1][iy] +
        frac_x_d0f[1] * frac_y_d0f[1] * tab_f_[ix + 1][iy + 1];

    delete [] frac_x_d0f;
    delete [] frac_y_d0f;

    return f_out;
}

int InterCSpline2D::get_index_x(double x_in) {
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

int InterCSpline2D::get_index_y(double y_in) {
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
