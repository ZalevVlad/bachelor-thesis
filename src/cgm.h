#ifndef CGM_H_
#define CGM_H_

#include <cmath>
#include <vector>
#include <iostream>
using namespace std;


double vec_vec(vector<double>& X, vector<double>& Y);
double nev(vector<double>& p, vector<double>& pr);
void matrix_vector_i0(vector<double>& ggl, vector<double>& ggu, vector<long>& ig, vector<long>& jg, vector<double>& di, vector<double>& x, vector<double>& y);
void matrix_vector(vector<double>& ggl, vector<double>& ggu, vector<long>& ig, vector<long>& jg, vector<double>& di, vector<double>& x, vector<double>& y);



class CGM
{
public:
    void init(vector< long>& gi_s, vector< long>& gj_s, vector< double>& di_s, vector< double>& gg_s, vector< double>& rp_s);
    void solve(vector<double>& solution, int max_iter, double eps, double relax);
    void solve_msg(vector<double>& ggl, vector<long>& ia, vector<long>& ja, vector<double>& di, vector<double>& pr,
        vector<double>& q, int max_iter, double eps, double relax);

private:
    void make_LLT_decomposition();
    void mul_matrix(vector<double>& f, vector<double>& x);
    void solve_L(vector<double>& f, vector<double>& x);
    void solve_LT(vector<double>& f, vector<double>& x);
    void solve_LLT(vector<double>& f, vector<double>& x);
    double dot_prod(vector<double>& a, vector<double>& b);

    int n;
    vector< long> ia, ja;
    vector<double> di, gg, pr, r, x0, z, p, s;
    vector<double> L_di, L_gg;
};

#endif // CGM_H_
