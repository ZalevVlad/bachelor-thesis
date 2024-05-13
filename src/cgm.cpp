#include "cgm.h"

// Скалярное умножение векторов
double vec_vec(vector<double>& X, vector<double>& Y) {
  double result = 0;
  for (int i = 0; i < X.size(); i++)
    result += X[i] * Y[i];
  return result;
}

// Вычисление невязки
double nev(vector<double>& p, vector<double>& pr) {
  double neva = sqrt(vec_vec(p, p)) / sqrt(vec_vec(pr, pr));
  return neva;
}

// Индексация с 0
void matrix_vector_i0(vector<double>& ggl, vector<double>& ggu, vector<long>& ig, vector<long>& jg, vector<double>& di, vector<double>& x, vector<double>& y) {
  vector<double> Z;
  Z = x;
  int n = x.size();
  for (int i = 0; i < n; i++)
    Z[i] = x[i] * di[i];
  for (int i = 0; i < n; i++)
    for (int j = ig[i]; j < ig[i + 1]; j++) {
      Z[i] += x[jg[j]] * ggl[j];
      Z[jg[j]] += x[i] * ggu[j];
    }
  y = Z;
}

// Индексация с 1
void matrix_vector(vector<double>& ggl, vector<double>& ggu, vector<long>& ig, vector<long>& jg, vector<double>& di, vector<double>& x, vector<double>& y) {
  vector<double> Z;
  Z = x;
  int n = x.size();
  for (int i = 0; i < n; i++)
    Z[i] = x[i] * di[i];
  for (int i = 0; i < n; i++)
    for (int j = ig[i] - 1; j < ig[i + 1] - 1; j++) {
      Z[i] += x[jg[j] - 1] * ggl[j];
      Z[jg[j] - 1] += x[i] * ggu[j];
    }
  y = Z;
}
void matrix_vector(vector<double>& ggl, vector<double>& ggu, vector<int>& ig, vector<int>& jg, vector<double>& di, vector<double>& x, vector<double>& y) {
  vector<double> Z;
  Z = x;
  int n = x.size();
  for (int i = 0; i < n; i++)
    Z[i] = x[i] * di[i];
  for (int i = 0; i < n; i++)
    for (int j = ig[i] - 1; j < ig[i + 1] - 1; j++) {
      Z[i] += x[jg[j] - 1] * ggl[j];
      Z[jg[j] - 1] += x[i] * ggu[j];
    }
  y = Z;
}

void CGM::init(vector<long>& gi_s, vector<long>& gj_s, vector<double>& di_s, vector<double>& gg_s, vector<double>& rp_s) {
  ia = gi_s;
  ja = gj_s;
  // Перейдем с нумерации с 1 на нумерацию с 0
  for (int i = 0; i < ia.size(); i++) {
    ia[i]--;
  }
  for (int i = 0; i < ja.size(); i++) {
    ja[i]--;
  }

  di = di_s;
  gg = gg_s;
  pr = rp_s;
  n = di.size();

  unsigned int m = ia[n];
  r.resize(n);
  x0.resize(n);
  z.resize(n);
  p.resize(n);
  s.resize(n);

  L_di.resize(n);
  L_gg.resize(m);

  for (unsigned int i = 0; i < n; i++) {
    L_di[i] = di[i];
    x0[i] = 0;  // Начальное приближение
  }

  for (unsigned int i = 0; i < m; i++)
    L_gg[i] = gg[i];
}

void CGM::make_LLT_decomposition() {
  double sum_d, sum_l;

  for (unsigned int k = 0; k < n; k++) {
    sum_d = 0;
    unsigned int i_s = ia[k], i_e = ia[k + 1];

    for (unsigned int i = i_s; i < i_e; i++) {
      sum_l = 0;
      unsigned int j_s = ia[ja[i]], j_e = ia[ja[i] + 1];

      for (unsigned int m = i_s; m < i; m++) {
        for (unsigned int j = j_s; j < j_e; j++) {
          if (ja[m] == ja[j]) {
            sum_l += L_gg[m] * L_gg[j];
            j_s++;
          }
        }
      }
      L_gg[i] = (L_gg[i] - sum_l) / L_di[ja[i]];

      sum_d += L_gg[i] * L_gg[i];
    }
    L_di[k] = sqrt(L_di[k] - sum_d);
  }
}

double CGM::dot_prod(vector<double>& a, vector<double>& b) {
  double d_p = 0;
  for (unsigned int i = 0; i < n; i++)
    d_p += a[i] * b[i];
  return d_p;
}

void CGM::mul_matrix(vector<double>& f, vector<double>& x) {
  for (unsigned int i = 0; i < n; i++) {
    double v_el = f[i];
    x[i] = di[i] * v_el;
    for (unsigned int k = ia[i], k1 = ia[i + 1]; k < k1; k++) {
      unsigned int j = ja[k];
      x[i] += gg[k] * f[j];
      x[j] += gg[k] * v_el;
    }
  }
}

void CGM::solve_L(vector<double>& f, vector<double>& x) {
  for (unsigned int k = 1, k1 = 0; k <= n; k++, k1++) {
    double sum = 0;

    for (unsigned int i = ia[k1]; i < ia[k]; i++)
      sum += L_gg[i] * x[ja[i]];

    x[k1] = (f[k1] - sum) / L_di[k1];
  }
}

void CGM::solve_LT(vector<double>& f, vector<double>& x) {
  for (unsigned int k = n, k1 = n - 1; k > 0; k--, k1--) {
    x[k1] = f[k1] / L_di[k1];
    double v_el = x[k1];

    for (unsigned int i = ia[k1]; i < ia[k]; i++)
      f[ja[i]] -= L_gg[i] * v_el;
  }
}

void CGM::solve_LLT(vector<double>& f, vector<double>& x) {
  solve_L(f, x);
  solve_LT(x, x);
}

void CGM::solve(vector<double>& solution, int max_iter, double eps, double relax) {
  // int a=0;
  mul_matrix(x0, r);
  make_LLT_decomposition();

  for (unsigned int i = 0; i < n; i++) {
    r[i] = pr[i] - r[i];
    /* if (r[i] != 0) {
         a++;
     }*/
  }

  solve_LLT(r, z);
  for (unsigned int i = 0; i < n; i++)
    p[i] = z[i];

  double alpha, betta, prod_1, prod_2;
  double discr, rp_norm;

  rp_norm = sqrt(dot_prod(pr, pr));

  prod_1 = dot_prod(p, r);

  bool end = false;

  int max_iters = 0;
  double nevazka_0 = 1, nevazka_1 = 2;
  vector<double> buff(di.size());
  vector<double> x_past(x0.size());
  for (unsigned int iter = 0; iter < max_iter && !end; iter++) {
    max_iters++;
    discr = sqrt(dot_prod(r, r));

    if (discr != discr || rp_norm != rp_norm)
      cerr << "Error: NaN detected!" << endl;

    // Способ вычисления невязки который был в решателе
    /*nevazka_0 = nevazka_1;
    nevazka_1 = discr / rp_norm;*/

    if (nevazka_1 > eps && nevazka_0 != nevazka_1) {
      mul_matrix(z, s);

      alpha = prod_1 / dot_prod(s, z);

      x_past = x0;
      for (unsigned int i = 0; i < n; i++) {
        x0[i] += alpha * z[i];
        r[i] -= alpha * s[i];
      }
      ////Релаксация
      for (int i = 0; i < x0.size(); i++) {
        x0[i] = relax * x0[i] + (1. - relax) * x_past[i];
      }

      solve_LLT(r, p);
      prod_2 = dot_prod(p, r);

      betta = prod_2 / prod_1;

      prod_1 = prod_2;

      for (unsigned int i = 0; i < n; i++)
        z[i] = p[i] + betta * z[i];

      // Вычисление невязки
      // Используется файл solver.h
      matrix_vector_i0(gg, gg, ia, ja, di, x0, buff);  // A*x
      for (int i = 0; i < di.size(); i++) {
        buff[i] -= pr[i];
      }
      nevazka_0 = nevazka_1;
      nevazka_1 = nev(buff, pr);
    } else
      end = true;

    // printf("CGM %d : nev= %e    \n", iter, nevazka_1);
  }
  printf("CGM : nev= %e ; iterations= %d\n", nevazka_1, max_iters);

  solution.resize(x0.size());
  for (unsigned int i = 0; i < n; i++)
    solution[i] = x0[i];
}

void CGM::solve_msg(vector<double>& ggl, vector<long>& ia, vector<long>& ja, vector<double>& di, vector<double>& pr,
                    vector<double>& q, int max_iter, double eps, double relax) {
  int n = di.size();
  // Начальное приближение
  vector<double> x0(n), x1(n);
  for (int i = 0; i < n; x1[i++] = 1);
  x0 = x1;

  vector<double> buff(di.size());
  matrix_vector(ggl, ggl, ia, ja, di, x1, buff);
  vector<double> r(di.size());
  for (int i = 0; i < r.size(); i++) {
    r[i] = pr[i] - buff[i];
  }
  vector<double> z = r;

  int iter = 0;
  double alpha, betta, nevazka_1 = 1., nevazka_0 = 0;
  while (iter < max_iter && nevazka_1 > eps /*&& nevazka_0!=nevazka_1*/) {
    nevazka_0 = nevazka_1;
    matrix_vector(ggl, ggl, ia, ja, di, z, buff);  // A*z
    alpha = vec_vec(r, r) / vec_vec(buff, z);
    for (int i = 0; i < n; i++) {
      x1[i] = x1[i] + alpha * z[i];
    }
    betta = 1 / vec_vec(r, r);
    for (int i = 0; i < n; i++) {
      r[i] = r[i] - alpha * buff[i];
    }
    betta *= vec_vec(r, r);
    for (int i = 0; i < n; i++) {
      z[i] = r[i] + betta * z[i];
    }

    // релаксация
    /*for (int i = 0; i < n; i++) { x1[i] = w_relax * x1[i] + (1. - w_relax) * x0[i]; }
    x0 = x1;*/

    matrix_vector(ggl, ggl, ia, ja, di, x1, buff);  // A*x
    for (int i = 0; i < n; i++) {
      buff[i] -= pr[i];
    }

    nevazka_1 = nev(buff, pr);
    iter++;
    // printf("nev = %e, iterations = %d\n", nevazka_1, iter);
  }
  q = x1;
  printf("MSG : nev = %e, iterations = %d\n", nevazka_1, iter);
}
