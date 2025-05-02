#ifndef _INTERCSPLINE2D_H_
#define _INTERCSPLINE2D_H_

#include<stdlib.h>

/* implementation of
 * natural cubic spline interpolation
 * for functions in 2-dimension */
class InterCSpline2D {
  private :

    /* number of bins in x
     * For ix = 0 ... nbin_x_,
     *   tab_x_[ix] = x
     * xmin_ = tab_x_[0] and
     * xmax_ = tab_x_[nbin_x_], respectively. */
    int nbin_x_;
    double xmin_;
    double xmax_;

    /* number of bins in y
     * For iy = 0 ... nbin_y_,
     *   tab_y_[iy] = y
     * ymin_ = tab_y_[0] and
     * ymax_ = tab_y_[nbin_y_], respectively. */
    int nbin_y_;
    double ymin_;
    double ymax_;

    double *tab_x_;
    double *tab_y_;

    /* tabulated function and its second derivatives
     * For ix = 0 ... nbin_x_ and
     *     iy = 0 ... nbin_y_,
     *   tab_f_[ix][iy] = f
     *   tab_d2f_dx_dx_[ix][iy] = d^{2}f / dx^{2}
     *   tab_d2f_dy_dy_[ix][iy] = d^{2}f / dy^{2}
     * at x = tab_x_[ix], y = tab_y_[iy] */
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

    /* initialized the 2D interpolator
     * with tabulated function at discrete set of points
     *
     * nbin_in_x : number of bins in x
     * nbin_in_y : number of bins in y
     * x_in : array for x
     * y_in : array for y
     * f_in : array for function f
     * For ix = 0 ... nbin_in_x and
     *     iy = 0 ... nbin_in_y,
     *   nbin_x_ = nbin_in_x
     *   nbin_y_ = nbin_in_y
     *   tab_x_[ix] = x_in[ix]
     *   tab_y_[iy] = y_in[iy]
     *   tab_f_[ix][iy] = f_in[ix][iy]
     * Note that the array will be sorted
     * in ascending order of x and y
     *   tab_x_[ix] < tab_x_[ix + 1]
     *   tab_y_[iy] < tab_y_[iy + 1]
     *
     * bc_df_dx : boundary condition for the first derivative in x
     *   bc_df_dx[iy][0] = df / dx at x = xmin_ (= tab_x_[0]),
     *                                y = tab_y_[iy]
     *   bc_df_dx[iy][1] = df / dx at x = xmax_ (= tab_x_[nbin_x_]),
     *                                y = tab_y_[iy]
     * bc_df_dy : boundary condition for the first derivative in y
     *   bc_df_dy[ix][0] = df / dy at x = tab_x_[ix],
     *                                y = ymin_ (= tab_y_[0])
     *   bc_df_dy[ix][1] = df / dy at x = tab_x_[ix],
     *                                y = ymax_ (= tab_y_[nbin_y_])
     * If a NULL pointer is given
     * for the value of bc_df_dx and/or bc_df_dy,
     * it performs a natural cubic spline.
     *   d^{2}f / dx^{2} = 0 at x = xmin_ and x = xmax_ and/or
     *   d^{2}f / dy^{2} = 0 at y = ymin_ and y = ymax_ */
    void init(int nbin_in_x,
              int nbin_in_y,
              double *x_in,
              double *y_in,
              double **f_in,
              double **bc_df_dx = NULL,
              double **bc_df_dy = NULL);

    /* get the interpolated function
     *
     * x_in : value of x at which the function is evaluated
     * y_in : value of y at which the function is evaluated
     * ptr_df_dx_out : pointer to the first derivative in x
     *                 df / dx
     * ptr_df_dy_out : pointer to the first derivative in y
     *                 df / dy
     * ptr_d2f_dx_dx_out : pointer to the second derivative
     *                     d^{2}f / dx^{2}
     * ptr_d2f_dx_dy_out : pointer to the second derivative
     *                     d^{2}f / dx dy
     * ptr_d2f_dy_dy_out : pointer to the second derivative
     *                     d^{2}f / dy^{2}
     * If a NULL pointer is provided,
     * it does not calculate the derivative.
     *
     * returns value of the interpolated function f
     * at x = x_in, y = y_in */
    double get_func(double x_in,
                    double y_in,
                    double *ptr_df_dx_out = NULL,
                    double *ptr_df_dy_out = NULL,
                    double *ptr_d2f_dx_dx_out = NULL,
                    double *ptr_d2f_dx_dy_out = NULL,
                    double *ptr_d2f_dy_dy_out = NULL);

    /* returns value of the function f
     * given by a linear interpolation */
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
