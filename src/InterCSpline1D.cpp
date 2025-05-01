#include<stdio.h>
#include"InterCSpline1D.h"

void InterCSpline1D::init(int nbin_in,
                          double *x_in,
                          double *f_in,
                          double *bc_df_dx) {
    reset();

    if (nbin_in < 1) {
        return;
    }

    nbin_ = nbin_in;
    tab_x_ = new double[nbin_ + 1];
    tab_f_ = new double[nbin_ + 1];
    tab_d2f_dx_ = new double[nbin_ + 1];

    for (int ix = 0; ix <= nbin_; ix++) {
        tab_x_[ix] = x_in[ix];
        tab_f_[ix] = f_in[ix];
    }

    bool sorted = false;
    while (!sorted) {
        bool swapped = false;
        for (int ix = 0; ix < nbin_; ix++) {
            if (tab_x_[ix] > tab_x_[ix + 1]) {
                double x_l = tab_x_[ix + 1];
                double f_l = tab_f_[ix + 1];

                double x_u = tab_x_[ix];
                double f_u = tab_f_[ix];

                tab_x_[ix] = x_l;
                tab_f_[ix] = f_l;

                tab_x_[ix + 1] = x_u;
                tab_f_[ix + 1] = f_u;

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted = true;
        }
    }

    xmin_ = tab_x_[0];
    xmax_ = tab_x_[nbin_];

    double *coeff_a = new double[nbin_ - 1];
    double *coeff_b = new double[nbin_ - 1];
    double *coeff_c = new double[nbin_ - 1];
    double *coeff_d = new double[nbin_ - 1];

    for (int ix = 1; ix < nbin_; ix++) {
        coeff_a[ix - 1] = (tab_x_[ix] - tab_x_[ix - 1]) / 6.;
        coeff_b[ix - 1] = (tab_x_[ix + 1] - tab_x_[ix - 1]) / 3.;
        coeff_c[ix - 1] = (tab_x_[ix + 1] - tab_x_[ix]) / 6.;
        coeff_d[ix - 1] =
            (tab_f_[ix + 1] - tab_f_[ix]) /
            (tab_x_[ix + 1] - tab_x_[ix]) -
            (tab_f_[ix] - tab_f_[ix - 1]) /
            (tab_x_[ix] - tab_x_[ix - 1]);
    }

    if (bc_df_dx != NULL) {
        /*
        coeff_b[0] * tab_d2f_dx_[1] +
        coeff_c[0] * tab_d2f_dx_[2]
            = coeff_d[0] -
              coeff_a[0] * tab_d2f_dx_[0]
        bc_df_dx[0] =
            (tab_f_[1] - tab_f_[0]) /
            (tab_x_[1] - tab_x_[0]) -
            tab_d2f_dx_[1] *
                (tab_x_[1] - tab_x_[0]) / 6. -
            tab_d2f_dx_[0] *
                (tab_x_[1] - tab_x_[0]) / 3.
        */
        coeff_b[0] -= 0.5 * coeff_a[0];
        coeff_d[0] +=
            3. * coeff_a[0] *
            (bc_df_dx[0] - (tab_f_[1] - tab_f_[0]) /
                           (tab_x_[1] - tab_x_[0])) /
            (tab_x_[1] - tab_x_[0]);

        /*
        coeff_a[nbin_ - 2] * tab_d2f_dx_[nbin_ - 2] +
        coeff_b[nbin_ - 2] * tab_d2f_dx_[nbin_ - 1]
            = coeff_d[nbin_ - 2] -
              coeff_c[nbin_ - 2] * tab_d2f_dx_[nbin_]
        bc_df_dx[1] =
            (tab_f_[nbin_] - tab_f_[nbin_ - 1]) /
            (tab_x_[nbin_] - tab_x_[nbin_ - 1]) +
            tab_d2f_dx_[nbin_] *
                (tab_x_[nbin_] - tab_x_[nbin_ - 1]) / 3. +
            tab_d2f_dx_[nbin_ - 1] *
                (tab_x_[nbin_] - tab_x_[nbin_ - 1]) / 6.
        */
        coeff_b[nbin_ - 2] -= 0.5 * coeff_c[nbin_ - 2];
        coeff_d[nbin_ - 2] -=
            3. * coeff_c[nbin_ - 2] *
            (bc_df_dx[1] - (tab_f_[nbin_] - tab_f_[nbin_ - 1]) /
                           (tab_x_[nbin_] - tab_x_[nbin_ - 1])) /
            (tab_x_[nbin_] - tab_x_[nbin_ - 1]);
    }

    double *coeff_ap = new double[nbin_ - 1];
    double *coeff_bp = new double[nbin_ - 1];
    double *coeff_cp = new double[nbin_ - 1];
    double *coeff_dp = new double[nbin_ - 1];

    coeff_ap[0] = 0.;
    coeff_bp[0] = 1.;
    coeff_cp[0] = coeff_c[0] / coeff_b[0];
    coeff_dp[0] = coeff_d[0] / coeff_b[0];
    for (int ix = 1; ix < nbin_ - 1; ix++) {
        coeff_ap[ix] = 0.;
        coeff_bp[ix] = 1.;
        coeff_cp[ix] =
            coeff_c[ix] /
            (coeff_b[ix] - coeff_cp[ix - 1] * coeff_a[ix]);
        coeff_dp[ix] =
            (coeff_d[ix] - coeff_dp[ix - 1] * coeff_a[ix]) /
            (coeff_b[ix] - coeff_cp[ix - 1] * coeff_a[ix]);
    }

    tab_d2f_dx_[nbin_ - 1] = coeff_dp[nbin_ - 2];
    if (bc_df_dx == NULL) {
        tab_d2f_dx_[nbin_] = 0.;
    } else {
        tab_d2f_dx_[nbin_] =
            3. * (bc_df_dx[1] -
                    (tab_f_[nbin_] - tab_f_[nbin_ - 1]) /
                    (tab_x_[nbin_] - tab_x_[nbin_ - 1])) /
            (tab_x_[nbin_] - tab_x_[nbin_ - 1]) -
            0.5 * tab_d2f_dx_[nbin_ - 1];
    }

    for (int ix = nbin_ - 2; ix > 0; ix--) {
        tab_d2f_dx_[ix] =
            coeff_dp[ix - 1] -
            coeff_cp[ix - 1] * tab_d2f_dx_[ix + 1];
    }

    if (bc_df_dx == NULL) {
        tab_d2f_dx_[0] = 0.;
    } else {
        tab_d2f_dx_[0] =
            -3. * (bc_df_dx[0] -
                    (tab_f_[1] - tab_f_[0]) /
                    (tab_x_[1] - tab_x_[0])) /
            (tab_x_[1] - tab_x_[0]) -
            0.5 * tab_d2f_dx_[1];
    }

    initialized_ = true;

    return;
}

double InterCSpline1D::get_func(double x_in,
                                double *ptr_df_dx_out,
                                double *ptr_d2f_dx_dx_out) {
    if (!initialized_) {
        if (ptr_df_dx_out != NULL) {
            *ptr_df_dx_out = 0.;
        }

        if (ptr_d2f_dx_dx_out != NULL) {
            *ptr_d2f_dx_dx_out = 0.;
        }

        return 0.;
    }

    int ix = get_index_x(x_in);

    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_f = tab_f_[ix + 1] - tab_f_[ix];

    double *frac_d0f = new double[2];
    frac_d0f[0] =
        (tab_x_[ix + 1] - x_in) / delta_x;
    frac_d0f[1] = 1. - frac_d0f[0];

    double *frac_d2f = new double[2];
    frac_d2f[0] =
        frac_d0f[0] *
        (frac_d0f[0] * frac_d0f[0] - 1.) *
        delta_x * delta_x / 6.;
    frac_d2f[1] =
        frac_d0f[1] *
        (frac_d0f[1] * frac_d0f[1] - 1.) *
        delta_x * delta_x / 6.;

    double f_out =
        frac_d0f[0] * tab_f_[ix] +
        frac_d0f[1] * tab_f_[ix + 1] +
        frac_d2f[0] * tab_d2f_dx_[ix] +
        frac_d2f[1] * tab_d2f_dx_[ix + 1];

    if (ptr_df_dx_out != NULL) {
        *ptr_df_dx_out =
            delta_f / delta_x -
            (3. * frac_d0f[0] * frac_d0f[0] - 1.) *
                delta_x * tab_d2f_dx_[ix] / 6. +
            (3. * frac_d0f[1] * frac_d0f[1] - 1.) *
                delta_x * tab_d2f_dx_[ix + 1] / 6.;
    }

    if (ptr_d2f_dx_dx_out != NULL) {
        *ptr_d2f_dx_dx_out =
            frac_d0f[0] * tab_d2f_dx_[ix] +
            frac_d0f[1] * tab_d2f_dx_[ix + 1];
    }

    return f_out;
}

double InterCSpline1D::get_func_lin(double x_in) {
    if (!initialized_) {
        return 0.;
    }

    int ix = get_index_x(x_in);

    double delta_x = tab_x_[ix + 1] - tab_x_[ix];

    double *frac_x_d0f = new double [2];
    frac_x_d0f[0] = (tab_x_[ix + 1] - x_in) / delta_x;
    frac_x_d0f[1] = 1. - frac_x_d0f[0];

    double f_out =
        frac_x_d0f[0] * tab_f_[ix] +
        frac_x_d0f[1] * tab_f_[ix + 1];

    delete [] frac_x_d0f;

    return f_out;
}

int InterCSpline1D::get_index_x(double x_in) {
    int ix = 0;
    if (x_in < xmin_) {
        ix = 0;
    } else if (x_in >= xmax_) {
        ix = nbin_ - 1;
    } else {
        for (int jx = 0; jx < nbin_; jx++) {
            if (x_in >= tab_x_[jx] && x_in < tab_x_[jx + 1]) {
                ix = jx;
                break;
            }
        }
    }

    return ix;
}
