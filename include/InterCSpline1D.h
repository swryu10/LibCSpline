#ifndef _INTERCSPLINE1D_H_
#define _INTERCSPLINE1D_H_

#include<stdlib.h>

/* implementation of
 * natural cubic spline interpolation
 * for functions in 1-dimension */
class InterCSpline1D {
  private :

    /* number of bins and tabulated functions
     * For ix = 0 ... nbin_,
     *   tab_x_[ix] = x
     *   tab_f_[ix] = f at x
     *   tab_d2f_dx_[ix] = d^{2}f / dx^{2} at x
     * xmin_ = tab_x_[0] and
     * xmax_ = tab_x_[nbin_], respectively. */
    int nbin_;
    double xmin_;
    double xmax_;
    double *tab_x_;
    double *tab_f_;
    double *tab_d2f_dx_;

    bool initialized_;

  public :

    InterCSpline1D() {
        initialized_ = false;

        return;
    }

    ~InterCSpline1D() {
        reset();

        return;
    }

    void reset() {
        if (!initialized_) {
            return;
        }

        delete [] tab_x_;
        delete [] tab_f_;
        delete [] tab_d2f_dx_;

        initialized_ = false;

        return;
    }

    /* initialized the 1D interpolator
     * with tabulated function at discrete set of points
     *
     * nbin_in : number of bins
     * x_in : array for x
     * f_in : array for function f
     * For ix = 0 ... nbin_in
     *   nbin_ = nbin_in,
     *   tab_x_[ix] = x_in[ix]
     *   tab_f_[ix] = f_in[ix]
     * Note that the array will be sorted
     * in ascending order of x
     *   tab_x_[ix] < tab_x_[ix + 1]
     *
     * bc_df_dx : boundary condition for the first derivative
     *   bc_df_dx[0] = df / dx at x = xmin_ (= tab_x_[0])
     *   bc_df_dx[1] = df / dx at x = xmax_ (= tab_x_[nbin_])
     * If a NULL pointer is given for the value of bc_df_dx,
     * it performs a natural cubic spline.
     *   d^{2}f / dx^{2} = 0 at x = xmin_ and x = xmax_ */
    void init(int nbin_in,
              double *x_in,
              double *f_in,
              double *bc_df_dx = NULL);

    /* get the interpolated function
     *
     * x_in : value of x at which the function is evaluated
     * ptr_df_dx_out : pointer to the first derivative
     * ptr_d2f_dx_dx_out : pointer to the second derivative
     * If a NULL pointer is provided,
     * it does not calculate the derivative.
     *
     * returns value of the interpolated function f
     * at x = x_in */
    double get_func(double x_in,
                    double *ptr_df_dx_out = NULL,
                    double *ptr_d2f_dx_dx_out = NULL);

    /* returns value of the function f
     * given by a linear interpolation */
    double get_func_lin(double x_in);

    int get_index_x(double x_in);
};

#endif
