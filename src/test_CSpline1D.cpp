#include<stdio.h>
#include<math.h>
#include<string>
#include"InterCSpline1D.h"

int nbin_x = 8;
double xmin = 0.;
double xmax = 2. * M_PI;

double n_pt_x = 128;

int main(int argc, char *argv[]) {
    double *tab_x = new double[nbin_x + 1];
    for (int ix = 0; ix <= nbin_x; ix++) {
        tab_x[ix] = xmin +
            (xmax - xmin) * static_cast<double>(ix) /
                            static_cast<double>(nbin_x);
    }

    double *tab_f = new double [nbin_x + 1];
    for (int ix = 0; ix <= nbin_x; ix++) {
            tab_f[ix] = cos(tab_x[ix]);
    }

    double *tab_bc_df_dx = new double[2];
    tab_bc_df_dx[0] = 0.;
    tab_bc_df_dx[1] = 0.;

    InterCSpline1D csp_test;
    csp_test.init(nbin_x, tab_x,
                  tab_f,
                  tab_bc_df_dx);

    delete [] tab_f;
    delete [] tab_x;

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
        double x_now = xmin +
            (xmax - xmin) * static_cast<double>(ix) /
                            static_cast<double>(n_pt_x);

        double df_dx_out = 0.;
        double d2f_dx_dx_out = 0.;
        double f_out =
            csp_test.get_func(x_now,
                              &df_dx_out, &d2f_dx_dx_out);

        fprintf(ptr_fout, "    %e    %e    %e    %e\n",
                x_now, f_out, df_dx_out, d2f_dx_dx_out);
    }

    if (have_out_ext) {
        fclose(ptr_fout);
    }

    return 0;
}
