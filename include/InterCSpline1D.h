#ifndef _INTERCSPLINE1D_H_
#define _INTERCSPLINE1D_H_

#include<stdlib.h>

/* implementation of
 * natural cubic spline interpolation
 * for functions in 1-dimension */
class InterCSpline1D {
  private :

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

    void init(int nbin_in,
              double *x_in,
              double *f_in,
              double *bc_df_dx = NULL);

    double get_func(double x_in,
                    double *ptr_df_dx_out = NULL,
                    double *ptr_d2f_dx_dx_out = NULL);

    double get_func_lin(double x_in);

    int get_index_x(double x_in);
};

#endif
