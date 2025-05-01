#ifndef _INTERCSPLINE2D_H_
#define _INTERCSPLINE2D_H_

#include<stdlib.h>

/* implementation of
 * natural cubic spline interpolation
 * for functions in 2-dimension */
class InterCSpline2D {
  private :

    int nbin_x_;
    double xmin_;
    double xmax_;

    int nbin_y_;
    double ymin_;
    double ymax_;

    double *tab_x_;
    double *tab_y_;

    double **tab_f_;
    double **tab_d2f_dx_dx_;
    double **tab_d2f_dy_dy_;

    bool initialized_;

  public :

    InterCSpline2D() {
        initialized_ = false;

        return;
    }

    ~InterCSpline2D() {
        reset();

        return;
    }

    void reset() {
        if (!initialized_) {
            return;
        }

        delete [] tab_x_;
        delete [] tab_y_;

        del_array_func(nbin_x_, nbin_y_,
                       tab_f_);
        del_array_func(nbin_x_, nbin_y_,
                       tab_d2f_dx_dx_);
        del_array_func(nbin_x_, nbin_y_,
                       tab_d2f_dy_dy_);

        initialized_ = false;

        return;
    }

    void init(int nbin_in_x,
              int nbin_in_y,
              double *x_in,
              double *y_in,
              double **f_in,
              double **bc_df_dx = NULL,
              double **bc_df_dy = NULL);

    double get_func(double x_in,
                    double y_in,
                    double *ptr_df_dx_out = NULL,
                    double *ptr_df_dy_out = NULL,
                    double *ptr_d2f_dx_dx_out = NULL,
                    double *ptr_d2f_dx_dy_out = NULL,
                    double *ptr_d2f_dy_dy_out = NULL);

    double get_func_lin(double x_in,
                        double y_in);

    int get_index_x(double x_in);
    int get_index_y(double y_in);

    static double **new_array_func(int nbin_in_x,
                                   int nbin_in_y) {
        double **ptr_out = new double *[nbin_in_x + 1];
        for (int ix = 0; ix <= nbin_in_x; ix++) {
            ptr_out[ix] = new double [nbin_in_y + 1];
            for (int iy = 0; iy <= nbin_in_y; iy++) {
                ptr_out[ix][iy] = 0.;
            }
        }

        return ptr_out;
    }
    static void del_array_func(int nbin_in_x,
                               int nbin_in_y,
                               double **array_in) {
        for (int ix = 0; ix <= nbin_in_x; ix++) {
            delete [] array_in[ix];
        }
        delete [] array_in;

        return;
    }

    static double **new_array_bc_df_dx(int nbin_in_y) {
        double **ptr_out = new double *[nbin_in_y + 1];
        for (int iy = 0; iy <= nbin_in_y; iy++) {
            ptr_out[iy] = new double[2];
            ptr_out[iy][0] = 0.;
            ptr_out[iy][1] = 0.;
        }

        return ptr_out;
    }
    static void del_array_bc_df_dx(int nbin_in_y,
                                   double **array_in) {
        for (int iy = 0; iy <= nbin_in_y; iy++) {
            delete [] array_in[iy];
        }
        delete [] array_in;

        return;
    }

    static double **new_array_bc_df_dy(int nbin_in_x) {
        double **ptr_out = new double *[nbin_in_x + 1];
        for (int ix = 0; ix <= nbin_in_x; ix++) {
            ptr_out[ix] = new double[2];
            ptr_out[ix][0] = 0.;
            ptr_out[ix][1] = 0.;
        }

        return ptr_out;
    }
    static void del_array_bc_df_dy(int nbin_in_x,
                                   double **array_in) {
        for (int ix = 0; ix <= nbin_in_x; ix++) {
            delete [] array_in[ix];
        }
        delete [] array_in;

        return;
    }
};

#endif
