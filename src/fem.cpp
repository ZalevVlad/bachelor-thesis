#include "fem.h"

#define I 1.
#define P_X 0
#define P_Y 0
#define P_Z 15.

double FEM::u(double x, double y, double z) {
  // return x * y + y * z + z * x;
  return x * y * z;
}
double FEM::f(double x, double y, double z) {
  return 0;
}
double FEM::u_rz(double r, double z) {
  return 10.;
}
double FEM::f_rz(double r, double z) {
  // if (r == 0) { return 0; }
  return 0.;
}

// Внесение элемента в разреженную матрицу
// ia, ja - профиль матрицы
// ggl - вектор куда вставляется значение
// str - номер строки начиная с 1
// col - номер столбца начиная с 1
// x - вносимое значение
void FEM::add_to_sparse(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x) {
  long start = ia[str - 1] - 1, end = ia[str] - 1;
  for (long i = start; i < end; i++) {
    if (ja[i] == col) {
      ggl[i] += x;
      break;
    }
  }
}

// Замена значения в разреженной матрице
// ia, ja - профиль матрицы
// ggl - вектор куда вставляется значение
// str - номер строки начиная с 1
// col - номер столбца начиная с 1
// x - вносимое значение
double FEM::replace_in_sparse_r(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x) {
  double el = 0;
  long start = ia[str - 1] - 1, end = ia[str] - 1;
  for (long i = start; i < end; i++) {
    if (ja[i] == col) {
      el = ggl[i];
      ggl[i] = x;
      break;
    }
  }
  return el;
}
// Заменит значение в разреженной матрице и вернет его
// ia, ja - профиль матрицы
// ggl - вектор куда вставляется значение
// str - номер строки начиная с 1
// col - номер столбца начиная с 1
// x - вносимое значение
void FEM::replace_in_sparse(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x) {
  long start = ia[str - 1] - 1, end = ia[str] - 1;
  for (long i = start; i < end; i++) {
    if (ja[i] == col) {
      ggl[i] = x;
      break;
    }
  }
}

void FEM::vec_vec_sum(std::vector<double>& a, std::vector<double>& b, std::vector<double>& x) {
  x.resize(a.size());
  for (int i = 0; i < a.size(); i++) {
    x[i] = a[i] + b[i];
  }
}
void FEM::matrix_vector_mul(std::vector<double>& q, std::vector<double>& x) {
  int n = q.size();
  x.resize(0);
  x.resize(n);
  for (int i = 0; i < n; i++) {
    x[i] = di[i] * q[i];
    for (unsigned int k = ia[i] - 1, k1 = ia[i + 1] - 1; k < k1; k++) {
      unsigned int j = ja[k] - 1;
      x[i] += ggl[k] * q[j];
      x[j] += ggl[k] * q[i];
    }
  }
}

// Посчитает и запишет в G[8][8] локальную матрицу жесткости
void FEM::local_G(vector<vector<double>>& G, double h_x, double h_y, double h_z, double lambda) {
  double k1 = lambda * h_y * h_z / (h_x * 36.), k2 = lambda * h_x * h_z / (h_y * 36.), k3 = lambda * h_x * h_y / (h_z * 36.);
  G[0] = {+k1 * 4 + k2 * 4 + k3 * 4};
  G[1] = {-k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4};
  G[2] = {+k1 * 2 - k2 * 4 + k3 * 2, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4};
  G[3] = {-k1 * 2 - k2 * 2 + k3 * 1, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4};
  G[4] = {+k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4};
  G[5] = {-k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4};
  G[6] = {+k1 * 1 - k2 * 2 - k3 * 2, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4};
  G[7] = {-k1 * 1 - k2 * 1 - k3 * 1, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4};
}
// Посчитает и запишет в G[4][4] локальную матрицу жесткости
void FEM::local_G_rzn(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i) {
  double k1 = h_z / h_r * lambda, k2 = h_r / h_z * lambda;
  G[0] = {2. / 6. * k1 + 2. / 6. * k2, -2. / 6. * k1 + 1. / 6. * k2, 1. / 6. * k1 - 2. / 6. * k2, -1. / 6. * k1 - 1. / 6. * k2};
  G[1] = {-2. / 6. * k1 + 1. / 6. * k2, 2. / 6. * k1 + 2. / 6. * k2, -1. / 6. * k1 - 1. / 6. * k2, 1. / 6. * k1 - 2. / 6. * k2};
  G[2] = {1. / 6. * k1 - 2. / 6. * k2, -1. / 6. * k1 - 1. / 6. * k2, 2. / 6. * k1 + 2. / 6. * k2, -2. / 6. * k1 + 1. / 6. * k2};
  G[3] = {-1. / 6. * k1 - 1. / 6. * k2, 1. / 6. * k1 - 2. / 6. * k2, -2. / 6. * k1 + 1. / 6. * k2, 2. / 6. * k1 + 2. / 6. * k2};
}
void FEM::local_G_rz(vector<vector<double>>& G, double hr, double hz, double lambda, int i) {
  double k1 = lambda * hz / (hr * 12.), k2 = lambda * hr / (hz * 12.);
  double r0 = ssrz.coord[ssrz.nvtr[i][0] - 1][0];
  double G24 = k1 * (2 * hr + 4 * r0), G14 = k2 * (hr + 4 * r0);
  double G112 = k1 * (hr + 2 * r0), G212 = k2 * (hr + 2 * r0);
  G[0] = {G24 + G14, -G24 + G212, G112 - G14, -G112 - G212};
  G[1] = {G[0][1], G[0][0], -G112 - G212, G112 - G14};
  G[2] = {G[0][2], G[1][2], G[0][0], -G24 + G212};
  G[3] = {G[0][3], G[1][3], G[2][3], G[0][0]};
}
void FEM::local_b_p(vector<double>& b, double h_x, double h_y, double h_z, int i) {
  double x0, x1, y0, y1, z0, z1;
  double f0, f1, f2, f3, f4, f5, f6, f7;
  double k = h_x * h_y * h_z / 216.;
  x0 = ssxyz.coord[ssxyz.nvtr[i][0] - 1][0];
  x1 = ssxyz.coord[ssxyz.nvtr[i][1] - 1][0];
  y0 = ssxyz.coord[ssxyz.nvtr[i][0] - 1][1];
  y1 = ssxyz.coord[ssxyz.nvtr[i][2] - 1][1];
  z0 = ssxyz.coord[ssxyz.nvtr[i][0] - 1][2];
  z1 = ssxyz.coord[ssxyz.nvtr[i][4] - 1][2];
  f0 = f(x0, y0, z0);
  f1 = f(x1, y0, z0);
  f2 = f(x0, y1, z0);
  f3 = f(x1, y1, z0);
  f4 = f(x0, y0, z1);
  f5 = f(x1, y0, z1);
  f6 = f(x0, y1, z1);
  f7 = f(x1, y1, z1);
  b[0] = k * (8 * f0 + 4 * f1 + 4 * f2 + 4 * f3 + 2 * f4 + 2 * f5 + 2 * f6 + 1 * f7);
  b[1] = k * (4 * f0 + 8 * f1 + 2 * f2 + 2 * f3 + 4 * f4 + 1 * f5 + 4 * f6 + 2 * f7);
  b[2] = k * (4 * f0 + 2 * f1 + 8 * f2 + 2 * f3 + 4 * f4 + 4 * f5 + 1 * f6 + 2 * f7);
  b[3] = k * (4 * f0 + 2 * f1 + 2 * f2 + 8 * f3 + 1 * f4 + 4 * f5 + 4 * f6 + 2 * f7);
  b[4] = k * (2 * f0 + 4 * f1 + 4 * f2 + 1 * f3 + 8 * f4 + 2 * f5 + 2 * f6 + 4 * f7);
  b[5] = k * (2 * f0 + 1 * f1 + 4 * f2 + 4 * f3 + 2 * f4 + 8 * f5 + 2 * f6 + 4 * f7);
  b[6] = k * (2 * f0 + 4 * f1 + 1 * f2 + 4 * f3 + 2 * f4 + 2 * f5 + 8 * f6 + 4 * f7);
  b[7] = k * (1 * f0 + 2 * f1 + 2 * f2 + 2 * f3 + 4 * f4 + 4 * f5 + 4 * f6 + 8 * f7);
}
void FEM::local_b_rz_pn(vector<double>& b, double h_r, double h_z, int i) {
  double k = h_r * h_z / 36.;
  double r0, r1, z0, z1;
  r0 = ssrz.coord[ssrz.nvtr[i][0] - 1][0];
  r1 = ssrz.coord[ssrz.nvtr[i][1] - 1][0];
  z0 = ssrz.coord[ssrz.nvtr[i][0] - 1][1];
  z1 = ssrz.coord[ssrz.nvtr[i][2] - 1][1];
  b[0] = k * (4 * f_rz(r0, z0) + 2 * f_rz(r1, z0) + 2 * f_rz(r0, z1) + 1 * f_rz(r1, z1));
  b[1] = k * (2 * f_rz(r0, z0) + 4 * f_rz(r1, z0) + 1 * f_rz(r0, z1) + 2 * f_rz(r1, z1));
  b[2] = k * (2 * f_rz(r0, z0) + 1 * f_rz(r1, z0) + 4 * f_rz(r0, z1) + 2 * f_rz(r1, z1));
  b[3] = k * (1 * f_rz(r0, z0) + 2 * f_rz(r1, z0) + 2 * f_rz(r0, z1) + 4 * f_rz(r1, z1));
}
void FEM::local_b_rz_p(vector<double>& b, double hr, double hz, int i) {
  double k = hr * hz / 72.;
  double r0, r1, z0, z1;
  r0 = ssrz.coord[ssrz.nvtr[i][0] - 1][0];
  r1 = ssrz.coord[ssrz.nvtr[i][1] - 1][0];
  z0 = ssrz.coord[ssrz.nvtr[i][0] - 1][1];
  z1 = ssrz.coord[ssrz.nvtr[i][2] - 1][1];
  b[0] = k * ((2 * hr + 8 * r0) * f_rz(r0, z0) + (2 * hr + 4 * r0) * f_rz(r1, z0) + (1 * hr + 4 * r0) * f_rz(r0, z1) + (1 * hr + 2 * r0) * f_rz(r1, z1));
  b[1] = k * ((2 * hr + 4 * r0) * f_rz(r0, z0) + (6 * hr + 8 * r0) * f_rz(r1, z0) + (1 * hr + 2 * r0) * f_rz(r0, z1) + (3 * hr + 4 * r0) * f_rz(r1, z1));
  b[2] = k * ((1 * hr + 4 * r0) * f_rz(r0, z0) + (1 * hr + 2 * r0) * f_rz(r1, z0) + (2 * hr + 8 * r0) * f_rz(r0, z1) + (2 * hr + 4 * r0) * f_rz(r1, z1));
  b[3] = k * ((1 * hr + 2 * r0) * f_rz(r0, z0) + (3 * hr + 4 * r0) * f_rz(r1, z0) + (2 * hr + 4 * r0) * f_rz(r0, z1) + (6 * hr + 8 * r0) * f_rz(r1, z1));
}

void FEM::build_matrix_profile() {
  A.resize(0);
  A.resize(mesh.kuzlov);

  // Пройдем по элементам и занесем узлы и их соседей в соответствующие позиции
  for (int i = 0; i < mesh.kel; i++) {
    for (int j = 0; j < FEM_XYZ_NODES_NUM; j++) {
      A[ssxyz.nvtr[i][j] - 1].insert(ssxyz.nvtr[i].begin(), ssxyz.nvtr[i].end());
    }
  }

  ia.resize(mesh.kuzlov + 1);
  ia[0] = 1;
  ia[1] = 1;
  list<long> ja_list;

  // Построим по списку узлов А профиль матрицы
  for (int i = 1; i < mesh.kuzlov; i++) {
    // Складываем в ja все элементы строки которые находятся под диагональю
    for (set<long>::iterator it = A[i].begin(); *it < i + 1; it++) {
      ja_list.push_back(*it);
    }
    // Добавляем в ia количество считанных элементов на момент текущей строки
    ia[i + 1] = ja_list.size() + 1;
  }

  // Скопируем значения с листа в вектор
  ja.resize(ja_list.size());
  int i = 0;
  for (list<long>::iterator it = ja_list.begin(); i < ja_list.size(); it++) {
    ja[i++] = *it;
  }
}
void FEM::build_matrix_profile_rz() {
  A.resize(0);
  A.resize(mesh_rz.kuzlov);

  // Пройдем по элементам и занесем узлы и их соседи в соответствующие позиции
  for (int i = 0; i < mesh_rz.kel; i++) {
    for (int j = 0; j < FEM_RZ_NODES_NUM; j++) {
      A[ssrz.nvtr[i][j] - 1].insert(ssrz.nvtr[i].begin(), ssrz.nvtr[i].end());
    }
  }

  ia.resize(mesh_rz.kuzlov + 1);
  ia[0] = 1;
  ia[1] = 1;
  list<long> ja_list;

  // Построим по списку узлов А профиль матрицы
  for (int i = 1; i < mesh_rz.kuzlov; i++) {
    // Складываем в ja все элементы строки которые находятся под диагональю
    for (set<long>::iterator it = A[i].begin(); *it < i + 1; it++) {
      ja_list.push_back(*it);
    }
    // Добавляем в ia количество считанных элементов на момент текущей строки
    ia[i + 1] = ja_list.size() + 1;
  }

  // Скопируем значения с листа в вектор
  ja.resize(ja_list.size());
  int i = 0;
  for (list<long>::iterator it = ja_list.begin(); i < ja_list.size(); it++) {
    ja[i++] = *it;
  }
}

void FEM::slae_init() {
  ggl.resize(0);
  di.resize(0);
  pr.resize(0);
  q.resize(0);

  ggl.resize(ja.size());
  di.resize(ia.size() - 1);
  pr.resize(di.size());
  q.resize(di.size());
  // Зададим начальное решение единицами
  for (int i = 0; i < q.size(); i++) {
    q[i] = 1.;
  }
}

void FEM::add_G() {
  int ggls = ggl.size();
  ggl.resize(0);
  di.resize(0);
  di.resize(mesh.kuzlov);
  ggl.resize(ggls);
  // Инициализация локального элемента
  vector<vector<double>> local_g;
  local_g.resize(FEM_XYZ_NODES_NUM);
  for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
    local_g[i].resize(FEM_XYZ_NODES_NUM);
  }

  // Основной цикл по элементам
  double h_x, h_y, h_z, lambda;
  vector<int> buff(FEM_XYZ_NODES_NUM);
  for (int i = 0; i < mesh.kel; i++) {
    // Сортируем номера локального элемента
    buff = ssxyz.nvtr[i];
    sort(buff.begin(), buff.end());

    // Посчитаем длину шага по x, y и z
    // Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
    h_x = ssxyz.coord[buff[1] - 1][0] - ssxyz.coord[buff[0] - 1][0];
    h_y = ssxyz.coord[buff[2] - 1][1] - ssxyz.coord[buff[0] - 1][1];
    h_z = ssxyz.coord[buff[4] - 1][2] - ssxyz.coord[buff[0] - 1][2];
    lambda = sreda.sigma[ssxyz.nvkat2d[i] - 1][SIGMA];

    // Посчитаем соответствующую локальную матрицу
    local_G(local_g, h_x, h_y, h_z, lambda);

    // Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
    // Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
    for (int l = 0; l < FEM_XYZ_NODES_NUM; l++) {
      for (int m = 0; m < FEM_XYZ_NODES_NUM; m++) {
        if (buff[m] < buff[l]) {
          add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_g[l][m]);
        } else if (buff[m] == buff[l]) {
          di[buff[m] - 1] += local_g[m][m];
        }
      }
    }
  }
}
void FEM::add_b() {
  // Номер узла источника
  int n_uzl = -1;
  // Ищем узел совпадающий с координатами источника
  for (int i = 0; i < mesh.kuzlov; i++) {
    // Сравниваем координаты узлов и источника (подразумевается что источник расположен в координате {0,0,0})
    if (ssxyz.coord[i][0] == 0 && ssxyz.coord[i][1] == 0 && ssxyz.coord[i][2] == 0) {
      n_uzl = i;
      break;
    }
  }
  if (n_uzl == -1) {
    printf("Source node not found\n");
  } else {
    pr[n_uzl] += sreda.current_sources[0][SOURCE_POW];
    printf("Source coord: %f %f %f, Node num = %d\n",
           ssxyz.coord[n_uzl][0], ssxyz.coord[n_uzl][1], ssxyz.coord[n_uzl][2], n_uzl + 1);
  }
  /* int n_el = n_el_xyz(0, 0, 0);
   for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
       pr[ssxyz.nvtr[n_el][i]-1]+= sreda.current_sources[0][SOURCE_POW]/FEM_XYZ_NODES_NUM;
   }*/
}
void FEM::add_G_sigma_diff() {
  // Инициализация локального элемента
  vector<vector<double>> local_g;
  local_g.resize(FEM_XYZ_NODES_NUM);
  for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
    local_g[i].resize(FEM_XYZ_NODES_NUM);
  }

  // Основной цикл по элементам
  double h_x, h_y, h_z, lambda;
  vector<int> buff(FEM_XYZ_NODES_NUM);
  for (int i = 0; i < mesh.kel; i++) {
    // Сортируем номера локального элемента
    buff = ssxyz.nvtr[i];
    sort(buff.begin(), buff.end());

    // Посчитаем длину шага по x, y и z
    // Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
    h_x = ssxyz.coord[buff[1] - 1][0] - ssxyz.coord[buff[0] - 1][0];
    h_y = ssxyz.coord[buff[2] - 1][1] - ssxyz.coord[buff[0] - 1][1];
    h_z = ssxyz.coord[buff[4] - 1][2] - ssxyz.coord[buff[0] - 1][2];
    lambda = sreda.sigma[ssxyz.nvkat2dr[i] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[i] - 1][SIGMA];

    if (lambda) {
      // Посчитаем соответствующую локальную матрицу
      local_G(local_g, h_x, h_y, h_z, lambda);
      // Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
      // Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
      for (int l = 0; l < FEM_XYZ_NODES_NUM; l++) {
        for (int m = 0; m < FEM_XYZ_NODES_NUM; m++) {
          if (buff[m] < buff[l]) {
            add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_g[l][m]);
          } else if (buff[m] == buff[l]) {
            di[buff[m] - 1] += local_g[m][m];
          }
        }
      }
    }
  }
}
void FEM::add_b_segregation_1() {
  /*vector<double> rz_to_xyz_solution;
  q_rz_to_xyz1(rz_to_xyz_solution);
  V1 = rz_to_xyz_solution;
  matrix_vector_mul(rz_to_xyz_solution,pr);*/
}
void FEM::add_b_segregation_2() {
  /*vector<double> rz_to_xyz_solution;
  q_rz_to_xyz2(rz_to_xyz_solution);
  V2 = rz_to_xyz_solution;
  matrix_vector_mul(rz_to_xyz_solution, pr);*/
}

void FEM::add_b_segregation_3() {
  const double Xi[3] = {-sqrt(0.6), 0, sqrt(0.6)};
  const double Ci[3] = {5. / 9., 8. / 9., 5. / 9.};
  std::vector<vector<vector<vector<double>>>> gradV(3);
  double sigma_diff, counter = 0;
  // Обьявление для производных в точках Гаусса трехмерного элемента
  for (int i = 0; i < 3; i++) {
    gradV[i].resize(3);
    for (int j = 0; j < 3; j++) {
      gradV[i][j].resize(3);
      for (int k = 0; k < 3; k++) {
        gradV[i][j][k].resize(3);
      }
    }
  }

  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // double dotx, doty, dotz;
      // sigma_diff = 1.;

      // Считаем производные в точках Гаусса
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            // dotx = cx + Xi[i] * hx / 2.; doty = cy + Xi[j] * hy / 2.; dotz = cz + Xi[k] * hz / 2.;
            grad_rz_in_point(cx + Xi[i] * hx / 2., cy + Xi[j] * hy / 2., cz + Xi[k] * hz / 2., gradV[i][j][k]);
            // gradV[i][j][k] = { doty * dotz,dotx * dotz,dotx * doty };
          }
        }
      }

      double dpsidx, dpsidy, dpsidz, sx, sy, sz, Gaus, f;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        Gaus = 0;
        // Неконстантная часть частных производных по базисным функциям
        if (nuzl % 2 == 0) {
          dpsidx = -1. / hx;
          sx = -1.;
        } else {
          dpsidx = 1. / hx;
          sx = 1.;
        }
        if (nuzl / 2 % 2 == 0) {
          dpsidy = -1. / hy;
          sy = -1.;
        } else {
          dpsidy = 1. / hy;
          sy = 1.;
        }
        if (nuzl / 4 % 2 == 0) {
          dpsidz = -1. / hz;
          sz = -1.;
        } else {
          dpsidz = 1. / hz;
          sz = 1.;
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              f = gradV[i][j][k][0] * dpsidx * (1 + sy * Xi[j]) * (1 + sz * Xi[k]) / 4. +
                  gradV[i][j][k][1] * dpsidy * (1 + sx * Xi[i]) * (1 + sz * Xi[k]) / 4. +
                  gradV[i][j][k][2] * dpsidz * (1 + sx * Xi[i]) * (1 + sy * Xi[j]) / 4.;
              /*f = dpsidx * (1 + sy * Xi[j]) * (1 + sz * Xi[k]) / 4. +
                  dpsidy * (1 + sx * Xi[i]) * (1 + sz * Xi[k]) / 4. +
                  dpsidz * (1 + sx * Xi[i]) * (1 + sy * Xi[j]) / 4.;*/
              Gaus += Ci[i] * Ci[j] * Ci[k] * f;
            }
          }
        }
        Gaus *= hx * hy * hz / 8.;
        pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus * sigma_diff;
        // pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus*3;
      }
    }
  }
  printf("non zero elements = %f\n", counter);
}
void FEM::add_b_segregation_4() {
  double sigma_diff, counter = 0;
  vector<double> gradV(3);
  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // Считаем производную в центре элемента
      grad_rz_in_point(cx, cy, cz, gradV);
      // gradV = { 1.00100025,1.00100025,1.00100025};
      double sx, sy, sz, Integral;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
        (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
        (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;

        Integral = sigma_diff * (gradV[0] * sx * hy * hz / 4. + gradV[1] * sy * hx * hz / 4. + gradV[2] * sz * hx * hy / 4.);
        pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
      }
    }
  }
  printf("non zero elements = %f\n", counter);
}
void FEM::add_b_segregation_5() {
  double sigma_diff, counter = 0;
  vector<vector<double>> gradV(FEM_XYZ_NODES_NUM);  // Вектор произвожных по x,y и z в каждом узле элемента
  for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
    gradV[i].resize(3);
  }
  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // Считаем производныe
      grad_rz_in_point(x1, y1, z1, gradV[0]);
      grad_rz_in_point(x2, y1, z1, gradV[1]);
      grad_rz_in_point(x1, y2, z1, gradV[2]);
      grad_rz_in_point(x2, y2, z1, gradV[3]);
      grad_rz_in_point(x1, y1, z2, gradV[4]);
      grad_rz_in_point(x2, y1, z2, gradV[5]);
      grad_rz_in_point(x1, y2, z2, gradV[6]);
      grad_rz_in_point(x2, y2, z2, gradV[7]);

      /*gradV[0] = { 1,1,1 };
      gradV[1] = { 1,1.001,1.001 };
      gradV[2] = { 1.001,1.,1.001};
      gradV[3] = { 1.001,1.,1.002001};
      gradV[4] = { 1.001,1.001,1.};
      gradV[5] = { 1.001,1.002001,1.001};
      gradV[6] = { 1.002001,1.001,1.001};
      gradV[7] = { 1.002001,1.001,1.002001};*/
      double sx, sy, sz, Integral;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
        (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
        (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;
        Integral = sigma_diff * (gradV[nuzl][0] * sx * hy * hz / 4. +
                                 gradV[nuzl][1] * sy * hx * hz / 4. +
                                 gradV[nuzl][2] * sz * hx * hy / 4.);
        if (Integral != Integral) {
          printf("Nan %f %f %f : %f %f %f\n", x1, y1, z1, x2, y2, z2);
        }
        pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
      }
    }
  }
  printf("non zero elements = %f\n", counter);
}

void FEM::add_b_segregation_3(int NOMEL) {
  const double Xi[3] = {-sqrt(0.6), 0, sqrt(0.6)};
  const double Ci[3] = {5. / 9., 8. / 9., 5. / 9.};
  std::vector<vector<vector<vector<double>>>> gradV(3);
  double sigma_diff, counter = 0;
  // Обьявление для производных в точках Гаусса трехмерного элемента
  for (int i = 0; i < 3; i++) {
    gradV[i].resize(3);
    for (int j = 0; j < 3; j++) {
      gradV[i][j].resize(3);
      for (int k = 0; k < 3; k++) {
        gradV[i][j][k].resize(3);
      }
    }
  }

  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // Считаем производные в точках Гаусса
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            grad_rz_in_point(cx + Xi[i] * hx / 2., cy + Xi[j] * hy / 2., cz + Xi[k] * hz / 2., gradV[i][j][k]);
            if (nel == NOMEL) {
              printf("i - %d, j - %d, k - %d  || dx = %e | dy = %e | dz = %e\n", i, j, k, gradV[i][j][k][0], gradV[i][j][k][1], gradV[i][j][k][2]);
            }
          }
        }
      }

      double dpsidx, dpsidy, dpsidz, sx, sy, sz, Gaus, f;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        Gaus = 0;
        // Неконстантная часть частных производных по базисным функциям
        if (nuzl % 2 == 0) {
          dpsidx = -1. / hx;
          sx = -1.;
        } else {
          dpsidx = 1. / hx;
          sx = 1.;
        }
        if (nuzl / 2 % 2 == 0) {
          dpsidy = -1. / hy;
          sy = -1.;
        } else {
          dpsidy = 1. / hy;
          sy = 1.;
        }
        if (nuzl / 4 % 2 == 0) {
          dpsidz = -1. / hz;
          sz = -1.;
        } else {
          dpsidz = 1. / hz;
          sz = 1.;
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              f = gradV[i][j][k][0] * dpsidx * (1 + sy * Xi[j]) * (1 + sz * Xi[k]) / 4. +
                  gradV[i][j][k][1] * dpsidy * (1 + sx * Xi[i]) * (1 + sz * Xi[k]) / 4. +
                  gradV[i][j][k][2] * dpsidz * (1 + sx * Xi[i]) * (1 + sy * Xi[j]) / 4.;
              /*f = dpsidx * (1 + sy * Xi[j]) * (1 + sz * Xi[k]) / 4. +
                  dpsidy * (1 + sx * Xi[i]) * (1 + sz * Xi[k]) / 4. +
                  dpsidz * (1 + sx * Xi[i]) * (1 + sy * Xi[j]) / 4.;*/
              Gaus += Ci[i] * Ci[j] * Ci[k] * f;
            }
          }
        }
        Gaus *= hx * hy * hz / 8.;
        pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus * sigma_diff;
        // pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus*3;
      }
    }
  }
  printf("non zero elements = %f\n", counter);
}
void FEM::add_b_segregation_4(int NOMEL) {
  double sigma_diff, counter = 0;
  vector<double> gradV(3);
  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // Считаем производную в центре элемента
      grad_rz_in_point(cx, cy, cz, gradV);
      if (nel == NOMEL) {
        printf("dx = %e | dy = %e | dz = %e\n", gradV[0], gradV[1], gradV[2]);
      }

      double sx, sy, sz, Integral;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
        (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
        (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;

        Integral = sigma_diff * (gradV[0] * sx * hy * hz / 4. + gradV[1] * sy * hx * hz / 4. + gradV[2] * sz * hx * hy / 4.);
        if (nel == NOMEL) {
          printf("local b %d = %e\n", nuzl, Integral);
        }
        pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
      }
    }
  }
  printf("non zero elements = %f\n", counter);
}
void FEM::add_b_segregation_5(int NOMEL) {
  double sigma_diff, counter = 0;
  vector<vector<double>> gradV(FEM_XYZ_NODES_NUM);  // Вектор произвожных по x,y и z в каждом узле элемента
  for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
    gradV[i].resize(3);
  }
  for (int nel = 0; nel < mesh.kel; nel++) {
    sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
    if (sigma_diff) {
      counter++;
      double
          x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0],
          x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
          y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
          z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
      double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
      double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

      // Считаем производныe
      grad_rz_in_point(x1, y1, z1, gradV[0]);
      grad_rz_in_point(x2, y1, z1, gradV[1]);
      grad_rz_in_point(x1, y2, z1, gradV[2]);
      grad_rz_in_point(x2, y2, z1, gradV[3]);
      grad_rz_in_point(x1, y1, z2, gradV[4]);
      grad_rz_in_point(x2, y1, z2, gradV[5]);
      grad_rz_in_point(x1, y2, z2, gradV[6]);
      grad_rz_in_point(x2, y2, z2, gradV[7]);

      double sx, sy, sz, Integral;
      for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
        if (nel == NOMEL) {
          printf("%d -  || dx = %e | dy = %e | dz = %e\n", nuzl, gradV[nuzl][0], gradV[nuzl][1], gradV[nuzl][2]);
        }

        (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
        (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
        (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;
        Integral = sigma_diff * (gradV[nuzl][0] * sx * hy * hz / 4. +
                                 gradV[nuzl][1] * sy * hx * hz / 4. +
                                 gradV[nuzl][2] * sz * hx * hy / 4.);
        if (Integral != Integral) { printf("Nan %f %f %f : %f %f %f\n", x1, y1, z1,