#include "fem.h"

#define I 1.
#define P_X 0
#define P_Y 0
#define P_Z 15.

double FEM::u(double x, double y, double z) {
    //return x * y + y * z + z * x;
    return x*y*z;
}
double FEM::f(double x, double y, double z) {
    return 0;
}
double FEM::u_rz(double r, double z) {
    return 10.;
}
double FEM::f_rz(double r, double z) {
    //if (r == 0) { return 0; }
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
        if (ja[i] == col) { ggl[i] += x; break; }
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
        if (ja[i] == col) { el = ggl[i];  ggl[i] = x; break; }
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
        if (ja[i] == col) { ggl[i] = x; break; }
    }
}

void FEM::vec_vec_sum(std::vector<double>& a, std::vector<double>& b, std::vector<double>& x) {
    x.resize(a.size());
    for (int i = 0; i < a.size(); i++) { x[i] = a[i] + b[i]; }
}
void FEM::matrix_vector_mul(std::vector<double>& q, std::vector<double>& x) {
    int n = q.size();
    x.resize(0);
    x.resize(n);
    for (int i = 0; i < n; i++)
    {
        x[i] = di[i] * q[i];
        for (unsigned int k = ia[i] - 1, k1 = ia[i + 1] - 1; k < k1; k++)
        {
            unsigned int j = ja[k] - 1;
            x[i] += ggl[k] * q[j];
            x[j] += ggl[k] * q[i];
        }
    }
}

// Посчитает и запишет в G[8][8] локальную матрицу жесткости
void FEM::local_G(vector<vector<double>>& G, double h_x, double h_y, double h_z, double lambda) {
    double k1 = lambda * h_y * h_z / (h_x*36.), k2 = lambda * h_x * h_z / (h_y*36.), k3 = lambda * h_x * h_y / (h_z*36.);
    G[0] = { +k1 * 4 + k2 * 4 + k3 * 4 };
    G[1] = { -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[2] = { +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[3] = { -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[4] = { +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[5] = { -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[6] = { +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 4 + k2 * 4 + k3 * 4 };
    G[7] = { -k1 * 1 - k2 * 1 - k3 * 1, +k1 * 1 - k2 * 2 - k3 * 2, -k1 * 2 + k2 * 1 - k3 * 2, +k1 * 2 + k2 * 2 - k3 * 4, -k1 * 2 - k2 * 2 + k3 * 1, +k1 * 2 - k2 * 4 + k3 * 2, -k1 * 4 + k2 * 2 + k3 * 2, +k1 * 4 + k2 * 4 + k3 * 4 };
}
// Посчитает и запишет в G[4][4] локальную матрицу жесткости
void FEM::local_G_rzn(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i) {
    double k1 = h_z / h_r * lambda, k2 = h_r / h_z * lambda;
    G[0] = { 2. / 6. * k1 + 2. / 6. * k2,-2. / 6. * k1 + 1. / 6. * k2,1. / 6. * k1 - 2. / 6. * k2,-1. / 6. * k1 - 1. / 6. * k2 };
    G[1] = { -2. / 6. * k1 + 1. / 6. * k2,2. / 6. * k1 + 2. / 6. * k2,-1. / 6. * k1 - 1. / 6. * k2,1. / 6. * k1 - 2. / 6. * k2 };
    G[2] = { 1. / 6. * k1 - 2. / 6. * k2,-1. / 6. * k1 - 1. / 6. * k2,2. / 6. * k1 + 2. / 6. * k2,-2. / 6. * k1 + 1. / 6. * k2 };
    G[3] = { -1. / 6. * k1 - 1. / 6. * k2,1. / 6. * k1 - 2. / 6. * k2,-2. / 6. * k1 + 1. / 6. * k2,2. / 6. * k1 + 2. / 6. * k2 };
}
void FEM::local_G_rz(vector<vector<double>>& G, double hr, double hz, double lambda, int i) {
    double k1 = lambda * hz / (hr * 12.), k2 = lambda * hr / (hz * 12.);
    double r0 = ssrz.coord[ssrz.nvtr[i][0] - 1][0];
    double G24 = k1 * (2 * hr + 4 * r0), G14 = k2 * (hr + 4 * r0);
    double G112 = k1 * (hr + 2 * r0), G212 = k2 * (hr + 2 * r0);
    G[0] = { G24 + G14,-G24 + G212,G112 - G14,-G112 - G212 };
    G[1] = { G[0][1],G[0][0],-G112 - G212,G112 - G14 };
    G[2] = { G[0][2],G[1][2],G[0][0],-G24+G212 };
    G[3] = { G[0][3],G[1][3],G[2][3],G[0][0] };
}
void FEM::local_b_p(vector<double>& b, double h_x, double h_y, double h_z, int i) {
    double x0, x1, y0, y1, z0,z1;
    double f0, f1, f2, f3, f4, f5, f6, f7;
    double k = h_x*h_y*h_z / 216.;
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
    r0 = ssrz.coord[ssrz.nvtr[i][0]-1][0];
    r1 = ssrz.coord[ssrz.nvtr[i][1]-1][0];
    z0 = ssrz.coord[ssrz.nvtr[i][0]-1][1];
    z1 = ssrz.coord[ssrz.nvtr[i][2]-1][1];
    b[0] = k*(4 * f_rz(r0, z0) + 2 * f_rz(r1, z0) + 2 * f_rz(r0, z1) + 1 * f_rz(r1, z1));
    b[1] = k*(2 * f_rz(r0, z0) + 4 * f_rz(r1, z0) + 1 * f_rz(r0, z1) + 2 * f_rz(r1, z1));
    b[2] = k*(2 * f_rz(r0, z0) + 1 * f_rz(r1, z0) + 4 * f_rz(r0, z1) + 2 * f_rz(r1, z1));
    b[3] = k*(1 * f_rz(r0, z0) + 2 * f_rz(r1, z0) + 2 * f_rz(r0, z1) + 4 * f_rz(r1, z1));
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

    //Построим по списку узлов А профиль матрицы
    for (int i = 1; i < mesh.kuzlov; i++) {
        //Складываем в ja все элементы строки которые находятся под диагональю
        for (set<long>::iterator it = A[i].begin(); *it < i + 1; it++) {
            ja_list.push_back(*it);
        }
        //Добавляем в ia количество считанных элементов на момент текущей строки
        ia[i + 1] = ja_list.size() + 1;
    }

    //Скопируем значения с листа в вектор
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

    //Построим по списку узлов А профиль матрицы
    for (int i = 1; i < mesh_rz.kuzlov; i++) {
        //Складываем в ja все элементы строки которые находятся под диагональю
        for (set<long>::iterator it = A[i].begin(); *it < i + 1; it++) {
            ja_list.push_back(*it);
        }
        //Добавляем в ia количество считанных элементов на момент текущей строки
        ia[i + 1] = ja_list.size() + 1;
    }

    //Скопируем значения с листа в вектор
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
    for (int i = 0; i < q.size(); i++) { q[i] = 1.; }
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
    for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { local_g[i].resize(FEM_XYZ_NODES_NUM); }

    // Основной цикл по элементам
    double h_x, h_y, h_z, lambda;
    vector<int> buff(FEM_XYZ_NODES_NUM);
    for (int i = 0; i < mesh.kel; i++) {
        //Сортируем номера локального элемента
        buff = ssxyz.nvtr[i];
        sort(buff.begin(), buff.end());

        //Посчитаем длину шага по x, y и z
        //Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
        h_x = ssxyz.coord[buff[1] - 1][0] - ssxyz.coord[buff[0] - 1][0];
        h_y = ssxyz.coord[buff[2] - 1][1] - ssxyz.coord[buff[0] - 1][1];
        h_z = ssxyz.coord[buff[4] - 1][2] - ssxyz.coord[buff[0] - 1][2];
        lambda = sreda.sigma[ssxyz.nvkat2d[i] - 1][SIGMA];

        //Посчитаем соответствующую локальную матрицу
        local_G(local_g, h_x, h_y, h_z, lambda);

        //Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
        //Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
        for (int l = 0; l < FEM_XYZ_NODES_NUM; l++) {
            for (int m = 0; m < FEM_XYZ_NODES_NUM; m++) {
                if (buff[m] < buff[l]) { add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_g[l][m]); }
                else if (buff[m] == buff[l]) {
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
        if (ssxyz.coord[i][0] == 0
            && ssxyz.coord[i][1] == 0
            && ssxyz.coord[i][2] == 0) {
            n_uzl = i;
            break;
        }
    }
    if (n_uzl == -1) {
        printf("Source node not found\n");
    }
    else {
        pr[n_uzl] += sreda.current_sources[0][SOURCE_POW];
        printf("Source coord: %f %f %f, Node num = %d\n",
            ssxyz.coord[n_uzl][0], ssxyz.coord[n_uzl][1], ssxyz.coord[n_uzl][2], n_uzl+1);
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
    for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { local_g[i].resize(FEM_XYZ_NODES_NUM); }

    // Основной цикл по элементам
    double h_x, h_y, h_z, lambda;
    vector<int> buff(FEM_XYZ_NODES_NUM);
    for (int i = 0; i < mesh.kel; i++) {
        //Сортируем номера локального элемента
        buff = ssxyz.nvtr[i];
        sort(buff.begin(), buff.end());

        //Посчитаем длину шага по x, y и z
        //Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
        h_x = ssxyz.coord[buff[1] - 1][0] - ssxyz.coord[buff[0] - 1][0];
        h_y = ssxyz.coord[buff[2] - 1][1] - ssxyz.coord[buff[0] - 1][1];
        h_z = ssxyz.coord[buff[4] - 1][2] - ssxyz.coord[buff[0] - 1][2];
        lambda = sreda.sigma[ssxyz.nvkat2dr[i] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[i]-1][SIGMA];
        
        if (lambda) {
            //Посчитаем соответствующую локальную матрицу
                local_G(local_g, h_x, h_y, h_z, lambda);
            //Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
            //Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
            for (int l = 0; l < FEM_XYZ_NODES_NUM; l++) {
                for (int m = 0; m < FEM_XYZ_NODES_NUM; m++) {
                    if (buff[m] < buff[l]) { add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_g[l][m]); }
                    else if (buff[m] == buff[l]) {
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
    const double Xi[3] = { -sqrt(0.6), 0,sqrt(0.6) };
    const double Ci[3] = { 5. / 9., 8. / 9.,5. / 9. };
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
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
                y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
                z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
            double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
            double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;
            
            //double dotx, doty, dotz;
            //sigma_diff = 1.;

            // Считаем производные в точках Гаусса
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        //dotx = cx + Xi[i] * hx / 2.; doty = cy + Xi[j] * hy / 2.; dotz = cz + Xi[k] * hz / 2.;
                        grad_rz_in_point(cx + Xi[i] * hx / 2., cy + Xi[j] * hy / 2., cz + Xi[k] * hz / 2., gradV[i][j][k]);
                        //gradV[i][j][k] = { doty * dotz,dotx * dotz,dotx * doty };
                    }
                }
            }

            double dpsidx, dpsidy, dpsidz, sx, sy, sz, Gaus, f;
            for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
                Gaus = 0;
                // Неконстантная часть частных производных по базисным функциям
                if (nuzl % 2 == 0) { dpsidx = -1. / hx; sx = -1.; }
                else { dpsidx = 1. / hx; sx = 1.; }
                if (nuzl / 2 % 2 == 0) { dpsidy = -1. / hy; sy = -1.; }
                else { dpsidy = 1. / hy; sy = 1.; }
                if (nuzl / 4 % 2 == 0) { dpsidz = -1. / hz; sz = -1.; }
                else { dpsidz = 1. / hz; sz = 1.; }

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
                //pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus*3;
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
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
                y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
                z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
            double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
            double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

            // Считаем производную в центре элемента
            grad_rz_in_point(cx, cy, cz, gradV);
            //gradV = { 1.00100025,1.00100025,1.00100025};
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
    vector<vector<double>> gradV(FEM_XYZ_NODES_NUM); // Вектор произвожных по x,y и z в каждом узле элемента
    for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { gradV[i].resize(3); }
    for (int nel = 0; nel < mesh.kel; nel++) {
        sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
        if (sigma_diff) {
            counter++;
            double
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
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
                Integral = sigma_diff * (
                    gradV[nuzl][0] * sx * hy * hz / 4. +
                    gradV[nuzl][1] * sy * hx * hz / 4. +
                    gradV[nuzl][2] * sz * hx * hy / 4.);
                if (Integral != Integral) { printf("Nan %f %f %f : %f %f %f\n", x1, y1, z1,x2,y2,z2); }
                pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
            }
        }
    }
    printf("non zero elements = %f\n", counter);
}

void FEM::add_b_segregation_3(int NOMEL) {
    const double Xi[3] = { -sqrt(0.6), 0,sqrt(0.6) };
    const double Ci[3] = { 5. / 9., 8. / 9.,5. / 9. };
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
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
                y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
                z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
            double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
            double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

            // Считаем производные в точках Гаусса
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        grad_rz_in_point(cx + Xi[i] * hx / 2., cy + Xi[j] * hy / 2., cz + Xi[k] * hz / 2., gradV[i][j][k]);
                        if (nel == NOMEL) { printf("i - %d, j - %d, k - %d  || dx = %e | dy = %e | dz = %e\n", i,j,k,gradV[i][j][k][0], gradV[i][j][k][1], gradV[i][j][k][2]); }
                    }
                }
            }

            double dpsidx, dpsidy, dpsidz, sx, sy, sz, Gaus, f;
            for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
                Gaus = 0;
                // Неконстантная часть частных производных по базисным функциям
                if (nuzl % 2 == 0) { dpsidx = -1. / hx; sx = -1.; }
                else { dpsidx = 1. / hx; sx = 1.; }
                if (nuzl / 2 % 2 == 0) { dpsidy = -1. / hy; sy = -1.; }
                else { dpsidy = 1. / hy; sy = 1.; }
                if (nuzl / 4 % 2 == 0) { dpsidz = -1. / hz; sz = -1.; }
                else { dpsidz = 1. / hz; sz = 1.; }

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
                //pr[ssxyz.nvtr[nel][nuzl] - 1] += Gaus*3;
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
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
                y1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][1], y2 = ssxyz.coord[ssxyz.nvtr[nel][2] - 1][1],
                z1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][2], z2 = ssxyz.coord[ssxyz.nvtr[nel][4] - 1][2];
            double hx = x2 - x1, hy = y2 - y1, hz = z2 - z1;
            double cx = (x2 + x1) / 2., cy = (y2 + y1) / 2., cz = (z2 + z1) / 2.;

            // Считаем производную в центре элемента
            grad_rz_in_point(cx, cy, cz, gradV);
            if (nel == NOMEL) { printf("dx = %e | dy = %e | dz = %e\n", gradV[0], gradV[1], gradV[2]); }

            double sx, sy, sz, Integral;
            for (int nuzl = 0; nuzl < FEM_XYZ_NODES_NUM; nuzl++) {
                (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
                (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
                (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;

                Integral = sigma_diff * (gradV[0] * sx * hy * hz / 4. + gradV[1] * sy * hx * hz / 4. + gradV[2] * sz * hx * hy / 4.);
                if (nel == NOMEL) { printf("local b %d = %e\n", nuzl, Integral); }
                pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
            }
        }
    }
    printf("non zero elements = %f\n", counter);
}
void FEM::add_b_segregation_5(int NOMEL) {
    double sigma_diff, counter = 0;
    vector<vector<double>> gradV(FEM_XYZ_NODES_NUM); // Вектор произвожных по x,y и z в каждом узле элемента
    for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { gradV[i].resize(3); }
    for (int nel = 0; nel < mesh.kel; nel++) {
        sigma_diff = sreda.sigma[ssxyz.nvkat2dr[nel] - 1][SIGMA] - sreda.sigma[ssxyz.nvkat2d[nel] - 1][SIGMA];
        if (sigma_diff) {
            counter++;
            double
                x1 = ssxyz.coord[ssxyz.nvtr[nel][0] - 1][0], x2 = ssxyz.coord[ssxyz.nvtr[nel][1] - 1][0],
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
                if (nel == NOMEL) { printf("%d -  || dx = %e | dy = %e | dz = %e\n",nuzl, gradV[nuzl][0], gradV[nuzl][1], gradV[nuzl][2]); }

                (nuzl % 2 == 0) ? sx = -1. : sx = 1.;
                (nuzl / 2 % 2 == 0) ? sy = -1. : sy = 1.;
                (nuzl / 4 % 2 == 0) ? sz = -1. : sz = 1;
                Integral = sigma_diff * (
                    gradV[nuzl][0] * sx * hy * hz / 4. +
                    gradV[nuzl][1] * sy * hx * hz / 4. +
                    gradV[nuzl][2] * sz * hx * hy / 4.);
                if (Integral != Integral) { printf("Nan %f %f %f : %f %f %f\n", x1, y1, z1, x2, y2, z2); }
                pr[ssxyz.nvtr[nel][nuzl] - 1] += Integral;
            }
        }
    }
    printf("non zero elements = %f\n", counter);
}

// Добавит матрцу жесткости и вектор правой части для двумерной задачи в разреженном строчно-столбцовом формате
void FEM::add_G_b_rz()
{
    // Инициализация локального элемента
    vector<vector<double>> local_G;
    local_G.resize(FEM_RZ_NODES_NUM);
    for (int i = 0; i < FEM_RZ_NODES_NUM; i++) { local_G[i].resize(FEM_RZ_NODES_NUM); }

    // Основной цикл по элементам
    double h_r, h_z, lambda;
    vector<int> buff(FEM_RZ_NODES_NUM);
    for (int i = 0; i < mesh_rz.kel; i++) {
        //Сортируем номера локального элемента
        buff = ssrz.nvtr[i];
        sort(buff.begin(), buff.end());

        //Посчитаем длину шага по x и y
        //Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
        h_r = ssrz.coord[buff[1] - 1][0] - ssrz.coord[buff[0] - 1][0];
        h_z = ssrz.coord[buff[2] - 1][1] - ssrz.coord[buff[0] - 1][1];
        lambda = sreda.sigma[ssrz.nvkat2d[i] - 1][SIGMA];

        //Посчитаем соответствующую локальную матрицу
        local_G_rz(local_G, h_r, h_z, lambda,i);

        //Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
        //Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
        for (int l = 0; l < FEM_RZ_NODES_NUM; l++) {
            for (int m = 0; m < FEM_RZ_NODES_NUM; m++) {
                if (buff[m] < buff[l]) { add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_G[l][m]); }
                else if (buff[m] == buff[l]) {
                    di[buff[m] - 1] += local_G[m][m];
                }
            }
        }

    }

    // Для учета правой части
    /*
    Для учета тока мы ставим ненулевое значение тока в первый узел

        ->  0-------------> r
            |         |
            |         |
            |         |
            |         |
            |----------
            |
            z
    Получим решение для I=1 затем, т.к. решение линейно по току, мы получим решение для необходимой силы тока
    умножением решения на константу
    */
    pr[0] = I / (2. * 3.141592653589793l);
}


// Добавит правую часть в вектор FEM.pr тест на полином
void FEM::add_G_b_p()
{
    // Инициализация локального элемента
    vector<vector<double>> local_g;
    local_g.resize(FEM_XYZ_NODES_NUM);
    for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { local_g[i].resize(FEM_XYZ_NODES_NUM); }

    //Инициализация локальной правой части
    vector<double> local_b;
    local_b.resize(FEM_XYZ_NODES_NUM);

    // Основной цикл по элементам
    double h_x, h_y, h_z, lambda;
    vector<int> buff(FEM_XYZ_NODES_NUM);
    for (int i = 0; i < mesh.kel; i++) {
        //Сортируем номера локального элемента
        buff = ssxyz.nvtr[i];
        sort(buff.begin(), buff.end());

        //Посчитаем длину шага по x, y и z
        //Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
        h_x = ssxyz.coord[buff[1] - 1][0] - ssxyz.coord[buff[0] - 1][0];
        h_y = ssxyz.coord[buff[2] - 1][1] - ssxyz.coord[buff[0] - 1][1];
        h_z = ssxyz.coord[buff[4] - 1][2] - ssxyz.coord[buff[0] - 1][2];
        lambda = sreda.sigma[ssxyz.nvkat2d[i] - 1][SIGMA];

        //Посчитаем соответствующую локальную матрицу
        local_G(local_g, h_x, h_y, h_z, lambda);

        //Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
        //Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
        for (int l = 0; l < FEM_XYZ_NODES_NUM; l++) {
            for (int m = 0; m < FEM_XYZ_NODES_NUM; m++) {
                if (buff[m] < buff[l]) { add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_g[l][m]); }
                else if (buff[m] == buff[l]) {
                    di[buff[m] - 1] += local_g[m][m];
                }
            }
        }

        local_b_p(local_b, h_x, h_y, h_z, i);
        for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) { pr[buff[i] - 1] += local_b[i]; }
    }
}
void FEM::add_G_b_rz_p()
{
    // Инициализация локального элемента
    vector<vector<double>> local_G;
    local_G.resize(FEM_RZ_NODES_NUM);
    for (int i = 0; i < FEM_RZ_NODES_NUM; i++) { local_G[i].resize(4); }

    //Инициализация локальной правой части
    vector<double> local_b;
    local_b.resize(FEM_RZ_NODES_NUM);

    // Основной цикл по элементам
    double h_x, h_y, lambda;
    vector<int> buff(FEM_RZ_NODES_NUM);
    for (int i = 0; i < mesh_rz.kel; i++) {
        //Сортируем номера локального элемента
        buff = ssrz.nvtr[i];
        sort(buff.begin(), buff.end());

        //Посчитаем длину шага по x и y
        //Для этого получим номера нижнего левого, нижнего правого и верхнего левого узлов и посчитаем разности их координат
        h_x = ssrz.coord[buff[1] - 1][0] - ssrz.coord[buff[0] - 1][0];
        h_y = ssrz.coord[buff[2] - 1][1] - ssrz.coord[buff[0] - 1][1];
        lambda = sreda.sigma[ssrz.nvkat2d[i] - 1][SIGMA];

        //Посчитаем соответствующую локальную матрицу
        local_G_rz(local_G, h_x, h_y, lambda,i);

        //Если номер столбца элемента меньше номера строки (лежит ниже диагонали), то вносим его в разреженню матрицу
        //Если номер столбца элемента равен номеру строки (лежит на диагонали), то складываем значение в диагональ
        for (int l = 0; l < FEM_RZ_NODES_NUM; l++) {
            for (int m = 0; m < FEM_RZ_NODES_NUM; m++) {
                if (buff[m] < buff[l]) { add_to_sparse(ia, ja, ggl, buff[l], buff[m], local_G[l][m]); }
                else if (buff[m] == buff[l]) {
                    di[buff[m] - 1] += local_G[m][m];
                }
            }
        }

        // Учет правой части
        local_b_rz_p(local_b,h_x,h_y,i);
        for (int i = 0; i < FEM_RZ_NODES_NUM; i++) { pr[buff[i] - 1] += local_b[i]; }
    }
}

void FEM::edge_cond_1() {
    // Для сохранения симметричности матрицы будем также занулять столбцы вместе со строками
    // при этом pr = pr - a*u, но т.к. u = 0 ( условие на удаленной границе), то правая часть не изменится при занулении столбцов
    for (int i = 0; i < ssxyz.l1.size(); i++) {
        int node = ssxyz.l1[i] - 1;
        di[ssxyz.l1[i] - 1] = 1.;
        pr[ssxyz.l1[i] - 1] = 0;

        // Зануляем строку нижнего треугольника ( и столбец верхнего треугольника т.к. матрица симметричная)
        long start = ia[ssxyz.l1[i] - 1] - 1, end = ia[ssxyz.l1[i]] - 1;
        for (long j = start; j < end; j++) {
            ggl[j] = 0;
        }

        // Зануляем столбец нижнего треугольника ( и строку верхнего треугольника т.к. матрицы симметричная)
        /*for (int j = ssxyz.l1[i] + 1; j <= di.size(); j++) {
            replace_in_sparse(ia, ja, ggl, j, ssxyz.l1[i], 0);
        }*/
        for (std::set<long>::iterator it = A[node].begin(); it != A[node].end(); it++) {
            if (node + 1 < *it) { // Если узел лежит в верхнем треугольнике
                replace_in_sparse(ia, ja, ggl, *it, node + 1, 0);
            }
        }

    }
}
void FEM::edge_cond_1_rz() {
    // Для сохранения симметричности матрицы будем также занулять столбцы вместе со строками
    // при этом pr = pr - a*u, но т.к. u = 0 ( условие на удаленной границе), то правая часть не изменится при занулении столбцов
    for (int i = 0; i < ssrz.l1.size(); i++) {
        int node = ssrz.l1[i] - 1;
        di[ssrz.l1[i] - 1] = 1.;
        pr[ssrz.l1[i] - 1] = 0;

        // Зануляем строку нижнего треугольника ( и столбец верхнего треугольника т.к. матрица симметричная)
        long start = ia[ssrz.l1[i] - 1] - 1, end = ia[ssrz.l1[i]] - 1;
        for (long j = start; j < end; j++) {
            ggl[j] = 0;
        }

        // Зануляем столбец нижнего треугольника ( и строку верхнего треугольника т.к. матрицы симметричная)
        /*for (int j = ssrz.l1[i] + 1; j <= di.size(); j++) {
            replace_in_sparse(ia, ja, ggl, j, ssrz.l1[i], 0);
        }*/
        for (std::set<long>::iterator it = A[node].begin(); it != A[node].end(); it++) {
            if (node + 1 < *it) { // Если узел лежит в верхнем треугольнике
                replace_in_sparse(ia, ja, ggl, *it, node + 1, 0);
            }
        }
    }
}

void FEM::print_matrix_plenum() {
    int nonzero = 0, dis = di.size(), ggls=ggl.size();
    for (int i = 0; i < ggl.size(); i++) { if (ggl[i]) { nonzero++; } }
    double ggl_matrix, nonzero_matrix, nonzero_ggl;
    ggl_matrix = (double)ggls/(double)(dis*dis);
    nonzero_matrix = (double)nonzero / (double)(dis * dis);
    nonzero_ggl = (double)nonzero /(double)ggls;
    printf("GGL/MATRIX = %f%% | NONZERO/MATRIX = %f%% | NONZERO/GGL = %.2f%% | ggl size = %d\n", 
        ggl_matrix * 100., nonzero_matrix * 100., nonzero_ggl * 100.,ggls);
}

void FEM::edge_cond_1_p() {
    vector<double> coords;
    double reminder, u_f;
    int node,col;

    // Для сохранения симметричности матрицы будем также занулять столбцы вместе со строками
    // при этом pr = pr - a*u
    for (int i = 0; i < ssxyz.l1.size(); i++) {
        node = ssxyz.l1[i] - 1;
        di[node] = 1.;
        coords = ssxyz.coord[node];
        u_f = u(coords[0], coords[1],coords[2]);

        // Зануляем строку нижнего треугольника ( и столбец верхнего треугольника т.к. матрица симметричная)
        // учитываем замененное значение в правой части
        long start = ia[node] - 1, end = ia[node+1] - 1;
        for (long j = start; j < end; j++) {
            pr[ja[j] - 1] -= ggl[j] * u_f;
            ggl[j] = 0;
        }

        // Зануляем столбец нижнего треугольника ( и строку верхнего треугольника т.к. матрицы симметричная)
        // учитываем замененное значение в правой части

        /*for (int j = ssxyz.l1[i] + 1; j <= di.size(); j++) {
            reminder = replace_in_sparse_r(ia, ja, ggl, j, ssxyz.l1[i], 0);
            pr[j - 1] -= reminder * u_f;
        }*/

        for (std::set<long>::iterator it = A[node].begin(); it!=A[node].end();it++){
            if (node + 1 < *it) { // Если узел лежит в верхнем треугольнике
                reminder = replace_in_sparse_r(ia, ja, ggl, *it, node + 1, 0);
                pr[*it-1] -= reminder * u_f;
            }
        }
    }

    // После того как мы учли зануленные столбцы в правой части
    // заменим правую часть первым краевым
    for (int i = 0; i < ssxyz.l1.size(); i++) {
        coords = ssxyz.coord[ssxyz.l1[i] - 1];
        pr[ssxyz.l1[i] - 1] = u(coords[0], coords[1],coords[2]);
    }
}
void FEM::edge_cond_1_rz_p() {
    vector<double> coords;
    double reminder, u;
    int node, col;

    // Для сохранения симметричности матрицы будем также занулять столбцы вместе со строками
    // при этом pr = pr - a*u
    for (int i = 0; i < ssrz.l1.size(); i++) {
        node = ssrz.l1[i] - 1;
        di[node] = 1.;
        coords = ssrz.coord[node];
        u = u_rz(coords[0], coords[1]);

        // Зануляем строку нижнего треугольника ( и столбец верхнего треугольника т.к. матрица симметричная)
        // учитываем замененное значение в правой части
        long start = ia[node] - 1, end = ia[node+1] - 1;
        for (long j = start; j < end; j++) {
            pr[ja[j]-1] -= ggl[j] * u;
            ggl[j] = 0;
        }

        // Зануляем столбец нижнего треугольника ( и строку верхнего треугольника т.к. матрицы симметричная)
        // учитываем замененное значение в правой части
        /*for (int j = ssrz.l1[i] + 1; j <= di.size(); j++) {
            reminder = replace_in_sparse_r(ia, ja, ggl, j, ssrz.l1[i], 0);
            pr[j - 1] -= reminder * u;
        }*/
        for (std::set<long>::iterator it = A[node].begin(); it != A[node].end(); it++) {
            if (node + 1 < *it) { // Если узел лежит в верхнем треугольнике
                reminder = replace_in_sparse_r(ia, ja, ggl, *it, node + 1, 0);
                pr[*it - 1] -= reminder * u;
            }
        }
    }

    // После того как мы учли зануленные столбцы в правой части
    // заменим правую часть первым краевым
    for (int i = 0; i < ssrz.l1.size(); i++) {
        coords = ssrz.coord[ssrz.l1[i] - 1];
        pr[ssrz.l1[i] - 1] = u_rz(coords[0], coords[1]);
    }
}

FEM::FEM(string filename) {
	sreda.read_problem(filename);
	mesh.gen_mesh(sreda);

    sreda_n.read_problem_n(filename);
    mesh_n.gen_mesh(sreda_n);
    mesh_n.add_mesh(mesh); // Прибавим узлы 

	sreda_rz.read_sreda(filename);
	mesh_rz.gen_mesh(sreda_rz, sreda);
    printf("Edge conditions: X -%f,%f | Y - %f,%f | Z0 - %d, Z1 - %d\n\n",
        sreda.edge_conditions[0], sreda.edge_conditions[1], sreda.edge_conditions[2],
        sreda.edge_conditions[3], sreda.edge_conditions[4], sreda.edge_conditions[5]);
}

void FEM::solve(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ--------------------------\n");
    ssxyz.gen_structures(sreda, mesh);

    unsigned int start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    add_G();
    add_b();
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    print_matrix_plenum();

    start = clock();
    edge_cond_1();
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    solve_cgm(q, max_iter, eps, relax);
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    solution_export(mesh,solution_filename,q);
    printf("END XYZ--------------------\n");
}

int FEM::find_nel(double x, double y, double z, class MESH& mesh, class STRCTRS& ss, std::vector<double>& q) {
    int res = 0;
    //static int xi = 0, yi = 0, zi = 0; // Индексы нижних границ в сетке
    if (mesh.x[0] <= x && x <= mesh.x[mesh.x.size() - 1] &&
        mesh.y[0] <= y && y <= mesh.y[mesh.y.size() - 1] &&
        mesh.z[0] <= z && z <= mesh.z[mesh.z.size() - 1]) {// Если точка лежит внутри
        while (!(mesh.x[xi] <= x && x <= mesh.x[xi + 1])) {
            if (x > mesh.x[xi + 1]) { xi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (x < mesh.x[xi]) { xi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        while (!(mesh.y[yi] <= y && y <= mesh.y[yi + 1])) {
            if (y > mesh.y[yi + 1]) { yi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (y < mesh.y[yi]) { yi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        while (!(mesh.z[zi] <= z && z <= mesh.z[zi + 1])) {
            if (z > mesh.z[zi + 1]) { zi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (z < mesh.z[zi]) { zi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        int n_el = (mesh.x.size() - 1) * (mesh.y.size() - 1) * zi + (mesh.x.size() - 1) * yi + xi; // Номер элемента в котором содержится точка
       /* double hx = mesh.x[xi + 1] - mesh.x[xi], hy = mesh.y[yi + 1] - mesh.y[yi], hz = mesh.z[zi + 1] - mesh.z[zi];
        double X1 = (mesh.x[xi + 1] - x) / hx, X2 = (x - mesh.x[xi]) / hx;
        double Y1 = (mesh.y[yi + 1] - y) / hy, Y2 = (y - mesh.y[yi]) / hy;
        double Z1 = (mesh.z[zi + 1] - z) / hz, Z2 = (z - mesh.z[zi]) / hz;
        vector<double> psi = {
            X1 * Y1 * Z1,X2 * Y1 * Z1,
            X1 * Y2 * Z1,X2 * Y2 * Z1,
            X1 * Y1 * Z2,X2 * Y1 * Z2,
            X1 * Y2 * Z2,X2 * Y2 * Z2 };*/
        for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
            printf("Qi = %e\n", q[ss.nvtr[n_el][i] - 1]);
        }
        res = n_el;
    }
    else {
        res = -1;
    }
    return res;
}
void FEM::solve_segregation1(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ: SEGREGATION 1--------------------------\n");
    unsigned int solve_segregation_start = clock();
    unsigned int start = clock();
    ssxyz.gen_structures(sreda, mesh);
    printf("gen_structures sreda -       %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    ssxyzn.gen_structures(sreda_n, mesh_n);
    printf("gen_structures sreda_n -     %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    // Сперва построим матрицу для sigma = sigma(normal) - sigma, ее мы умножим на вектор 
    add_G_sigma_diff();
    q_rz_to_xyz1(Vnb,mesh,ssxyz); // Получим решение двумерной задачи для трехмерной сетки из файла sreda
    matrix_vector_mul(Vnb, pr); // Умножим глобальную матрицу на вектор Vnb, это будет правая часть для задачи с выделением
    /*int bel = find_nel(P_X, P_Y, P_Z, mesh, ssxyz, Vnb);
    for (int i = 0; i < 8; i++) { printf("bi = %e\n", pr[ssxyz.nvtr[bel][i] - 1]); }
    printf("b in point: %e\n", solution_xyz_in_point(P_X, P_Y, P_Z, pr, mesh, ssxyz));
    xi = 0; yi = 0; zi = 0;*/
    //for (int i = 0; i < pr.size(); i++) { printf("pr%d = %e\n", i, pr[i]); }

    add_G(); // Матрица жесткости
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    edge_cond_1(); // Применение первых нулевых краевых условий
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    solve_cgm(Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    q_rz_to_xyz1(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    q_xyz_to_xyz(Va,mesh, ssxyz,Vas, mesh_n,ssxyzn); // Решение задачи с аномальным полем расширить до большей сетки 
    vec_vec_sum(Vas, Vns, Vs); // Сложим решения q = Va + Vн
    
    solution_export(mesh_n, solution_filename + "vn", Vns); // Выведем Vn
    solution_export(mesh_n, solution_filename + "va", Vas); // Выведем Va
    solution_export(mesh_n,solution_filename, Vs); // Выведем суммарное решение для общей сетки
    printf("SOLVE SEGREGATION -          %.2f sec\n", (double)(clock() - solve_segregation_start) / 1000.);
    printf("END XYZ: SEGREGATION 1--------------------\n");
}

void FEM::print_matrix() {
    int jaind = 0;
    for (int i = 0; i < 8; i++) {
        int elms = ia[i + 1] - ia[i];
        for (int j = 0; j < 8; j++) {
            if (elms) {
                if (j == ja[jaind] - 1) {
                    printf("%.3e ", ggl[jaind]);
                    jaind++;
                }
                else {
                    printf("%.3e ", 0);
                }
                elms--;
            }
            else {
                printf("%.3e ", 0);
            }
        }
        printf("\n");
    }
}

void FEM::solve_segregation2(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ: SEGREGATION 2--------------------------\n");
    unsigned int solve_segregation_start = clock();
    unsigned int start = clock();
    ssxyz.gen_structures(sreda, mesh);
    printf("gen_structures sreda -       %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    ssxyzn.gen_structures(sreda_n, mesh_n);
    printf("gen_structures sreda_n -     %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    // Сперва построим матрицу для sigma = sigma(normal) - sigma, ее мы умножим на вектор 
    add_G_sigma_diff();
    q_rz_to_xyz2(Vnb, mesh, ssxyz); // Получим решение двумерной задачи для трехмерной сетки из файла sreda
    //Vnb = { 1,1.001,1.001,1.002001,1.001,1.002001,1.002001,1.003003001};
    matrix_vector_mul(Vnb, pr); // Умножим глобальную матрицу на вектор Vnb, это будет правая часть для задачи с выделением
   /* xi = 0; yi = 0; zi = 0;
    int bel = find_nel(P_X, P_Y, P_Z, mesh, ssxyz, Vnb);
    for (int i = 0; i < 8; i++) { printf("bi = %e\n", pr[ssxyz.nvtr[bel][i] - 1]); }
    printf("b in point: %e\n", solution_xyz_in_point(P_X, P_Y, P_Z, pr, mesh, ssxyz));*/
    //for (int i = 0; i < pr.size(); i++) { printf("pr%d = %e\n", i, pr[i]); }

    add_G(); // Матрица жесткости
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    edge_cond_1(); // Применение первых нулевых краевых условий
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    solve_cgm(Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    q_rz_to_xyz2(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    q_xyz_to_xyz(Va, mesh, ssxyz, Vas, mesh_n, ssxyzn); // Решение задачи с аномальным полем расширить до большей сетки 
    vec_vec_sum(Vas, Vns, Vs); // Сложим решения q = Va + Vн

    solution_export(mesh_n, solution_filename + "vn", Vns); // Выведем Vn
    solution_export(mesh_n, solution_filename + "va", Vas); // Выведем Va
    solution_export(mesh_n, solution_filename, Vs); // Выведем суммарное решение для общей сетки
    printf("SOLVE SEGREGATION -          %.2f sec\n", (double)(clock() - solve_segregation_start) / 1000.);
    printf("END XYZ: SEGREGATION 2--------------------\n");
}
void FEM::solve_segregation3(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ: SEGREGATION 3--------------------------\n");
    unsigned int solve_segregation_start = clock();
    unsigned int start = clock();
    ssxyz.gen_structures(sreda, mesh);
    printf("gen_structures sreda -       %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    ssxyzn.gen_structures(sreda_n, mesh_n);
    printf("gen_structures sreda_n -     %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    add_G(); // Матрица жесткости
    add_b_segregation_3();
    //for (int i = 0; i < pr.size(); i++) { printf("pr%d = %e\n", i, pr[i]); }

    /*xi = 0; yi = 0; zi = 0;
    int NOMEL = find_nel(P_X, P_Y, P_Z, mesh, ssxyz, q);
    add_b_segregation_3(NOMEL);
    for (int i = 0; i < 8; i++) { printf("bi = %e\n", pr[ssxyz.nvtr[NOMEL][i] - 1]); }
    printf("b in point: %e\n", solution_xyz_in_point(P_X, P_Y, P_Z, pr, mesh, ssxyz));*/
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    edge_cond_1(); // Применение первых нулевых краевых условий
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    solve_cgm(Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    q_rz_to_xyz2(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    q_xyz_to_xyz(Va, mesh, ssxyz, Vas, mesh_n, ssxyzn); // Решение задачи с аномальным полем расширить до большей сетки 
    vec_vec_sum(Vas, Vns, Vs); // Сложим решения q = Va + Vн

    solution_export(mesh_n, solution_filename + "vn", Vns); // Выведем Vn
    //solution_export(mesh_n, solution_filename + "va", Vas); // Выведем Va
    solution_export(mesh, solution_filename + "va", Va); // Выведем Va
    solution_export(mesh_n, solution_filename, Vs); // Выведем суммарное решение для общей сетки
    printf("SOLVE SEGREGATION -          %.2f sec\n", (double)(clock() - solve_segregation_start) / 1000.);
    printf("END XYZ: SEGREGATION 3--------------------\n");
}
void FEM::solve_segregation4(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ: SEGREGATION 4--------------------------\n");
    unsigned int solve_segregation_start = clock();
    unsigned int start = clock();
    ssxyz.gen_structures(sreda, mesh);
    printf("gen_structures sreda -       %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    ssxyzn.gen_structures(sreda_n, mesh_n);
    printf("gen_structures sreda_n -     %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    add_G(); // Матрица жесткости
    add_b_segregation_4();
    //for (int i = 0; i < pr.size(); i++) { printf("pr%d = %e\n", i, pr[i]); }
    /*xi = 0; yi = 0; zi = 0;
    int NOMEL = find_nel(P_X, P_Y, P_Z, mesh, ssxyz, q);
    add_b_segregation_4(NOMEL);
    printf("x1 = %f, x2 = %f\n y1 = %f, y2 = %f\n, z1 = %f, z2 = %f\n",
        ssxyz.coord[ssxyz.nvtr[NOMEL][0] - 1][0], ssxyz.coord[ssxyz.nvtr[NOMEL][1] - 1][0],
        ssxyz.coord[ssxyz.nvtr[NOMEL][0] - 1][1], ssxyz.coord[ssxyz.nvtr[NOMEL][2] - 1][1],
        ssxyz.coord[ssxyz.nvtr[NOMEL][0] - 1][2], ssxyz.coord[ssxyz.nvtr[NOMEL][4] - 1][2]);
    xi = 0; yi = 0; zi = 0;
    for (int i = 0; i < 8; i++) { printf("bi = %e\n", pr[ssxyz.nvtr[NOMEL][i] - 1]); }
    printf("b in point: %e\n", solution_xyz_in_point(P_X, P_Y, P_Z, pr, mesh, ssxyz));*/
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    edge_cond_1(); // Применение первых нулевых краевых условий
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    //CGM solver;
    //solver.solve_msg(ggl,ia,ja,di,pr,Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    solve_cgm(Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    q_rz_to_xyz2(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    q_xyz_to_xyz(Va, mesh, ssxyz, Vas, mesh_n, ssxyzn); // Решение задачи с аномальным полем расширить до большей сетки 
    vec_vec_sum(Vas, Vns, Vs); // Сложим решения q = Va + Vн

    solution_export(mesh_n, solution_filename + "vn", Vns); // Выведем Vn
    //solution_export(mesh_n, solution_filename + "va", Vas); // Выведем Va
    solution_export(mesh, solution_filename + "va", Va); // Выведем Va
    solution_export(mesh_n, solution_filename, Vs); // Выведем суммарное решение для общей сетки
    printf("SOLVE SEGREGATION -          %.2f sec\n", (double)(clock() - solve_segregation_start) / 1000.);
    printf("END XYZ: SEGREGATION 4--------------------\n");
}
void FEM::solve_segregation5(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nSOLVE XYZ: SEGREGATION 5--------------------------\n");
    unsigned int solve_segregation_start = clock();
    unsigned int start = clock();
    ssxyz.gen_structures(sreda, mesh);
    printf("gen_structures sreda -       %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    ssxyzn.gen_structures(sreda_n, mesh_n);
    printf("gen_structures sreda_n -     %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    build_matrix_profile();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();

    start = clock();
    add_G(); // Матрица жесткости
    add_b_segregation_5();
    //for (int i = 0; i < pr.size(); i++) { printf("pr%d = %e\n", i, pr[i]); }
    /*xi = 0; yi = 0; zi = 0;
    int NOMEL = find_nel(P_X, P_Y, P_Z, mesh, ssxyz, q);
    add_b_segregation_5(NOMEL);
    xi = 0; yi = 0; zi = 0;
    for (int i = 0; i < 8; i++) { printf("bi = %e\n", pr[ssxyz.nvtr[NOMEL][i] - 1]); }
    printf("b in point: %e\n", solution_xyz_in_point(P_X, P_Y, P_Z, pr, mesh, ssxyz));*/
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    edge_cond_1(); // Применение первых нулевых краевых условий
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);

    start = clock();
    solve_cgm(Va, max_iter, eps, relax); // Решение задачи на добавочное поле в узлах sreda
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    q_rz_to_xyz2(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    q_xyz_to_xyz(Va, mesh, ssxyz, Vas, mesh_n, ssxyzn); // Решение задачи с аномальным полем расширить до большей сетки 
    vec_vec_sum(Vas, Vns, Vs); // Сложим решения q = Va + Vн

    solution_export(mesh_n, solution_filename + "vn", Vns); // Выведем Vn
    //solution_export(mesh_n, solution_filename + "va", Vas); // Выведем Va
    solution_export(mesh, solution_filename + "va", Va); // Выведем Va
    solution_export(mesh_n, solution_filename, Vs); // Выведем суммарное решение для общей сетки
    printf("SOLVE SEGREGATION -          %.2f sec\n", (double)(clock() - solve_segregation_start) / 1000.);
    printf("END XYZ: SEGREGATION 5--------------------\n");
}
 
void FEM::solve_rz(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("SOLVE RZ--------------------------\n");
    ssrz.gen_structures_rz(sreda, mesh_rz);

    unsigned int start = clock();
    build_matrix_profile_rz();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);

    slae_init();
    
    start = clock();
    add_G_b_rz();
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);

    //print_matrix_plenum();

    start = clock();
    edge_cond_1_rz();
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);
    
    start = clock();
    solve_cgm(q, max_iter, eps, relax);
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);

    qrz = q;
    solution_rz_export(solution_filename);
    printf("END RZ--------------------\n\n");
}


// Тесты на полином
void FEM::test_polinome(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("\nTEST POLINOME XYZ--------------------------\n");
    ssxyz.gen_structures_p(sreda,mesh);

    unsigned int start = clock();
    build_matrix_profile();
    printf("build matrix profile - %.2f sec\n", (double)(clock() - start)/1000.);
    
    slae_init();
    
    start = clock();
    add_G_b_p();
    printf("global G and b - %.2f sec\n", (double)(clock() - start) / 1000.);
    
    //print_matrix_plenum();
    
    start = clock();
    edge_cond_1_p();
    printf("first edge condition apply - %.2f sec\n", (double)(clock()-start) / 1000.);

    start = clock();
    solve_cgm(q, max_iter, eps, relax);
    printf("SLAE solution - %.2f sec\n", (double)(clock() - start) / 1000.);

    print_solution_p_miss();

    solution_export(mesh,solution_filename, q);
    printf("END TEST POLINOME XYZ--------------------\n");
}
void FEM::test_polinome_rz(int max_iter, double eps, double relax, std::string solution_filename) {
    printf("TEST POLINOME RZ--------------------------\n");
	ssrz.gen_structures_rz_p(sreda,mesh_rz);

    unsigned int start = clock();
    build_matrix_profile_rz();
    printf("build matrix profile -       %.2f sec\n", (double)(clock() - start) / 1000.);
    slae_init();
    start = clock();
    add_G_b_rz_p();
    printf("global G and b -             %.2f sec\n", (double)(clock() - start) / 1000.);
    print_matrix_plenum();
    start = clock();
    edge_cond_1_rz_p();
    printf("first edge condition apply - %.2f sec\n", (double)(clock() - start) / 1000.);
    start = clock();
    solve_cgm(q, max_iter, eps, relax);
    //solve_msg(q, max_iter, eps, relax);
    printf("SLAE solution -              %.2f sec\n", (double)(clock() - start) / 1000.);
    print_solution_rz_p_miss();
    qrz = q;
    solution_rz_export(solution_filename);
    printf("END TEST POLINOME RZ--------------------------\n\n");
}

void FEM::solve_cgm(vector<double>& q, int max_iter, double eps, double relax) {
    CGM cgm;
    cgm.init(ia, ja, di, ggl, pr);
    cgm.solve(q, max_iter, eps, relax);
}

void FEM::q_rz_to_xyz1(vector<double>& ans, class MESH& mesh, class STRCTRS& ss) {
    unsigned int start = clock();
    double x, y, z, r, sol;
    ans.resize(0);
    ans.resize(mesh.kuzlov);
    for (int i = 0; i < mesh.kel; i++) {
        x = (ss.coord[ss.nvtr[i][0] - 1][0]+ ss.coord[ss.nvtr[i][1] - 1][0])/2.;
        y = (ss.coord[ss.nvtr[i][0] - 1][1] + ss.coord[ss.nvtr[i][2] - 1][1]) / 2.;
        z = (ss.coord[ss.nvtr[i][0] - 1][2] + ss.coord[ss.nvtr[i][4] - 1][2]) / 2.;
        r = sqrtl(x * x + y * y);
        sol = solution_rz_in_point(r, z);
        for (int j = 0; j < FEM_XYZ_NODES_NUM; j++) {
            ans[ss.nvtr[i][j] - 1] = sol;
        }
    }
    printf("q_rz_to_xyz1   -             %.2f sec\n", (double)(clock() - start) / 1000.);
}

// ans - Вектор в который запишется решение
// mesh - Трехмерная сетка
// ss - Струткуры трехмерной сетки
void FEM::q_rz_to_xyz2(vector<double>& ans, class MESH& mesh, class STRCTRS& ss) {
    unsigned int start = clock();
    double x, y, z, r; // Координаты узла
    ans.resize(0);
    ans.resize(mesh.kuzlov);
    for (int i = 0; i < mesh.kuzlov; i++) {
        x = ss.coord[i][0];
        y = ss.coord[i][1];
        z = ss.coord[i][2];
        r = sqrtl(x * x + y * y);
        
        ans[i] = solution_rz_in_point(r, z);
    }
    printf("q_rz_to_xyz2  -              %.2f sec\n", (double)(clock() - start) / 1000.);
}

void FEM::q_xyz_to_xyz(vector<double>& q1, class MESH& mesh1, class STRCTRS& ss1, vector<double>& q2, class MESH& mesh2, class STRCTRS& ss2) {
    unsigned int start = clock();
    double x, y, z; // Координаты узла
    q2.resize(0);
    q2.resize(mesh2.kuzlov);
    xi = 0; yi = 0; zi = 0;
    for (int i = 0; i < mesh2.kuzlov; i++) {
        x = ss2.coord[i][0];
        y = ss2.coord[i][1];
        z = ss2.coord[i][2];
        //if (i % 10000 == 0)printf("%d %f %f %f ", i, x, y, z);
        q2[i] = solution_xyz_in_point(x,y,z,q1,mesh1,ss1);
        
    }
    printf("q_xyz_to_xyz  -              %.2f sec\n", (double)(clock() - start) / 1000.);
}

double norm(vector<double>& x) {
    double res = 0;
    for (int i = 0; i < x.size(); i++) { res += x[i] * x[i]; }
    return sqrtl(res);
}
void FEM::print_solution_p_miss() {
    /*vector<double> diff,solution;
    diff.resize(q.size());
    solution.resize(q.size());
    for (int i = 0; i < q.size(); i++) { solution[i] = u(ssxyz.coord[i][0], ssxyz.coord[i][1], ssxyz.coord[i][2]); }
    for (int i = 0; i < q.size(); i++) { diff[i] =solution[i]-q[i]; }
    printf("Norm of error vector: %e\n", norm(diff) / norm(solution));*/
    vector<double> diff, solution, solution_a;
    double x,y, z;
    diff.resize(mesh.kel);
    solution.resize(mesh.kel);
    solution_a.resize(mesh.kel);
    for (int i = 0; i < mesh.kel; i++) {
        // Для каждого элемента получаем решение в центре
        for (int j = 0; j < FEM_XYZ_NODES_NUM; j++) {
            solution[i] += q[ssxyz.nvtr[i][j] - 1] / 8.;
        }
        x = (ssxyz.coord[ssxyz.nvtr[i][0] - 1][0] + ssxyz.coord[ssxyz.nvtr[i][1] - 1][0]) / 2.;
        y = (ssxyz.coord[ssxyz.nvtr[i][0] - 1][1] + ssxyz.coord[ssxyz.nvtr[i][2] - 1][1]) / 2.;
        z = (ssxyz.coord[ssxyz.nvtr[i][0] - 1][2] + ssxyz.coord[ssxyz.nvtr[i][4] - 1][2]) / 2.;
        solution_a[i] = u(x,y, z);
    }
    for (int i = 0; i < solution.size(); i++) { diff[i] = solution[i] - solution_a[i]; }
    printf("Norm of error vector: %e\n", norm(diff) / norm(solution));
}
void FEM::print_solution_rz_p_miss() {
    /*vector<double> diff, solution;
    diff.resize(q.size());
    solution.resize(q.size());
    for (int i = 0; i < q.size(); i++) { solution[i] = u_rz(ssrz.coord[i][0], ssrz.coord[i][1]); }
    for (int i = 0; i < q.size(); i++) { diff[i] = solution[i] - q[i]; }
    printf("Norm of error vector: %e\n", norm(diff) / norm(solution));*/

    vector<double> diff, solution, solution_a;
    double r, z;
    diff.resize(mesh_rz.kel);
    solution.resize(mesh_rz.kel);
    solution_a.resize(mesh_rz.kel);
    for (int i = 0; i < mesh_rz.kel; i++) {
        // Для каждого элемента получаем решение в центре
        for (int j = 0; j < FEM_RZ_NODES_NUM; j++) {
            solution[i] += q[ssrz.nvtr[i][j] - 1] / 4.;
        }
        r = (ssrz.coord[ssrz.nvtr[i][0] - 1][0] + ssrz.coord[ssrz.nvtr[i][1] - 1][0]) / 2.;
        z = (ssrz.coord[ssrz.nvtr[i][0] - 1][1] + ssrz.coord[ssrz.nvtr[i][2] - 1][1]) / 2.;
        solution_a[i] = u_rz(r, z);
    }
    for (int i = 0; i < solution.size(); i++) { diff[i] = solution[i] - solution_a[i]; }
    printf("Norm of error vector: %e\n", norm(diff) / norm(solution));
}

double FEM::solution_rz_in_point(double r, double z) {
    double res = 0;
    static int ri = 0, zi = 0; // Индексы нижних границ отрезков в одномерных сетках
    if (mesh_rz.r[0] <= r && r <= mesh_rz.r[mesh_rz.r.size() - 1] && // Если точка лежит внутри
        mesh_rz.z[0] <= z && z <= mesh_rz.z[mesh_rz.z.size() - 1]) {

        while (!(mesh_rz.r[ri] <= r && r <= mesh_rz.r[ri + 1])) {
            if (r > mesh_rz.r[ri + 1]) { ri++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (r < mesh_rz.r[ri]) { ri--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        while (!(mesh_rz.z[zi] <= z && z <= mesh_rz.z[zi + 1])) {
            if (z > mesh_rz.z[zi + 1]) { zi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (z < mesh_rz.z[zi]) { zi--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        int n_el = (mesh_rz.r.size() - 1) * zi + ri;
        double hr = mesh_rz.r[ri + 1] - mesh_rz.r[ri], hh = mesh_rz.z[zi + 1] - mesh_rz.z[zi];

        double R1 = (mesh_rz.r[ri + 1] - r) / hr, R2 = (r - mesh_rz.r[ri]) / hr;
        double H1 = (mesh_rz.z[zi + 1] - z) / hh, H2 = (z - mesh_rz.z[zi]) / hh;
        vector<double> psi = { R1 * H1,R2 * H1,R1 * H2,R2 * H2 };
        for (int i = 0; i < FEM_RZ_NODES_NUM; i++) {
            res += psi[i] * qrz[ssrz.nvtr[n_el][i] - 1];
        }
    } else {
        res = NAN;
        printf("Error for point: r = %f, z = %f\n", r, z);
    }
    return res;
}

void FEM::grad_rz_in_point(double x, double y, double z, std::vector<double>& grad) {
    static int ri = 0, zi = 0; // Индексы нижних границ отрезков в одномерных сетках
    double r = sqrt(x*x+y*y);
    //if (r == 0) { r = 1e-30; }
    if (mesh_rz.r[0] <= r && r <= mesh_rz.r[mesh_rz.r.size() - 1] && // Если точка лежит внутри
        mesh_rz.z[0] <= z && z <= mesh_rz.z[mesh_rz.z.size() - 1]) {
        while (!(mesh_rz.r[ri] <= r && r <= mesh_rz.r[ri + 1])) {
            if (r > mesh_rz.r[ri + 1]) { ri++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            else if (r < mesh_rz.r[ri]) { ri--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        while (!(mesh_rz.z[zi] <= z && z <= mesh_rz.z[zi + 1])) {
            if (z > mesh_rz.z[zi + 1]) { zi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            else if (z < mesh_rz.z[zi]) { zi--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        int n_el = (mesh_rz.r.size() - 1) * zi + ri;
        double hr = mesh_rz.r[ri + 1] - mesh_rz.r[ri], hh = mesh_rz.z[zi + 1] - mesh_rz.z[zi];
        double R1 = (mesh_rz.r[ri + 1] - r) / hr, R2 = (r - mesh_rz.r[ri]) / hr;
        double H1 = (mesh_rz.z[zi + 1] - z) / hh, H2 = (z - mesh_rz.z[zi]) / hh;
        double dVdri[4] = { -H1 / hr ,H1 / hr,-H2 / hr,H2 / hr };
        double dVdzi[4] = { -R1 / hh ,-R2 / hh,R1 / hh,R2 / hh };
        //vector<double> psi = { R1 * H1,R2 * H1,R1 * H2,R2 * H2 };
        grad[0] = 0;
        grad[1] = 0;
        grad[2] = 0;
        for (int i = 0; i < FEM_RZ_NODES_NUM; i++) {
            grad[0] += dVdri[i] * qrz[ssrz.nvtr[n_el][i] - 1];
            grad[1] += dVdri[i] * qrz[ssrz.nvtr[n_el][i] - 1];
            grad[2] += dVdzi[i] * qrz[ssrz.nvtr[n_el][i] - 1];
        }
        if (r) {
            grad[0] *= x / r;
            grad[1] *= y / r;
        }
    }
    else {
        printf("Error for point: r = %f, z = %f\n", r, z);
    }
}
void FEM::grad_rz_in_point2(double x, double y, double z, std::vector<double>& grad) {
    static int ri = 0, zi = 0; // Индексы нижних границ отрезков в одномерных сетках
    double r = sqrt(x * x + y * y);
    //if (r == 0) { r = 1e-30; }
    if (mesh_rz.r[0] <= r && r <= mesh_rz.r[mesh_rz.r.size() - 1] && // Если точка лежит внутри
        mesh_rz.z[0] <= z && z <= mesh_rz.z[mesh_rz.z.size() - 1]) {
        while (!(mesh_rz.r[ri] <= r && r <= mesh_rz.r[ri + 1])) {
            if (r > mesh_rz.r[ri + 1]) { ri++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            else if (r < mesh_rz.r[ri]) { ri--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        while (!(mesh_rz.z[zi] <= z && z <= mesh_rz.z[zi + 1])) {
            if (z > mesh_rz.z[zi + 1]) { zi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            else if (z < mesh_rz.z[zi]) { zi--; } // Если точка левее текущего отрезка, то сместить отрезок влево
        }
        int n_el = (mesh_rz.r.size() - 1) * zi + ri;
        vector<double> V = { qrz[ssrz.nvtr[n_el][0] - 1],qrz[ssrz.nvtr[n_el][1] - 1],
            qrz[ssrz.nvtr[n_el][2] - 1],qrz[ssrz.nvtr[n_el][3] - 1] };
        double hr = mesh_rz.r[ri + 1] - mesh_rz.r[ri], hh = mesh_rz.z[zi + 1] - mesh_rz.z[zi];
        double R1 = (mesh_rz.r[ri + 1] - r) / hr, R2 = (r - mesh_rz.r[ri]) / hr;
        double H1 = (mesh_rz.z[zi + 1] - z) / hh, H2 = (z - mesh_rz.z[zi]) / hh;
        grad[0] = (H1 * (V[1] - V[0]) + H2 * (V[3] - V[2])) / hr;
        grad[1] = grad[0];
        grad[2] = (R1*(V[2]-V[0])+R2*(V[3]-V[1]))/hh;
        if (r) {
            grad[0] *= x / r;
            grad[1] *= y / r;
        }
    }
    else {
        printf("Error for point: r = %f, z = %f\n", r, z);
    }
}


int FEM::n_el_xyz(double x, double y, double z) {
    int n_el = -1;
    if (mesh.x[0] <= x && x <= mesh.x[mesh.x.size() - 1] && // Если точка лежит внутри
        mesh.y[0] <= y && y <= mesh.y[mesh.y.size() - 1] &&
        mesh.z[0] <= z && z <= mesh.z[mesh.z.size() - 1]) {
        int xi = 0, yi = 0, zi = 0; // Индексы нижних границ в сетке
        while (!(mesh.x[xi] <= x && x <= mesh.x[xi + 1])) { xi++; }
        while (!(mesh.y[yi] <= y && y <= mesh.y[yi + 1])) { yi++; }
        while (!(mesh.z[zi] <= z && z <= mesh.z[zi + 1])) { zi++; }
        n_el = (mesh.x.size() - 1) * (mesh.y.size() - 1) * zi + (mesh.x.size() - 1) * yi + xi;
    }
    return n_el;
}

double FEM::solution_xyz_in_point(double x, double y, double z) {
    double res = 0;
    if (mesh.x[0] <= x && x <= mesh.x[mesh.x.size() - 1] && // Если точка лежит внутри
        mesh.y[0] <= y && y <= mesh.y[mesh.y.size() - 1] &&
        mesh.z[0] <= z && z <= mesh.z[mesh.z.size() - 1]) {
        int xi = 0, yi = 0, zi = 0; // Индексы нижних границ в сетке
        while (!(mesh.x[xi] <= x && x <= mesh.x[xi + 1])) { xi++; }
        while (!(mesh.y[yi] <= y && y <= mesh.y[yi + 1])) { yi++; }
        while (!(mesh.z[zi] <= z && z <= mesh.z[zi + 1])) { zi++; }
        int n_el = (mesh.x.size() - 1) * (mesh.y.size() - 1) * zi + (mesh.x.size()-1) * yi + xi; // Номер элемента в котором содержится точка
        double hx = mesh.x[xi + 1] - mesh.x[xi], hy = mesh.y[yi + 1] - mesh.y[yi], hz = mesh.z[zi + 1] - mesh.z[zi];
        double X1 = (mesh.x[xi + 1] - x) / hx, X2 = (x - mesh.x[xi]) / hx;
        double Y1 = (mesh.y[yi + 1] - y) / hy, Y2 = (y - mesh.y[yi]) / hy;
        double Z1 = (mesh.z[zi + 1] - z) / hz, Z2 = (z - mesh.z[zi]) / hz;
        vector<double> psi = {
            X1 * Y1 * Z1,X2 * Y1 * Z1,
            X1 * Y2 * Z1,X2 * Y2 * Z1,
            X1 * Y1 * Z2,X2 * Y1 * Z2,
            X1 * Y2 * Z2,X2 * Y2 * Z2 };
        for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
            res += psi[i] * q[ssxyz.nvtr[n_el][i] - 1];
        }
    } else {
        res = NAN;
    }
    return res;
}

double FEM::solution_xyz_in_point(double x, double y, double z,
    std::vector<double>& q, class MESH& mesh, class STRCTRS& ss) {
    double res = 0;
    //static int xi = 0, yi = 0, zi = 0; // Индексы нижних границ в сетке
    if (mesh.x[0] <= x && x <= mesh.x[mesh.x.size() - 1] && 
        mesh.y[0] <= y && y <= mesh.y[mesh.y.size() - 1] &&
        mesh.z[0] <= z && z <= mesh.z[mesh.z.size() - 1]) {// Если точка лежит внутри
        while (!(mesh.x[xi] <= x && x <= mesh.x[xi + 1])) { 
            if (x > mesh.x[xi + 1]) { xi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (x < mesh.x[xi]) { xi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        while (!(mesh.y[yi] <= y && y <= mesh.y[yi + 1])) { 
            if (y > mesh.y[yi + 1]) { yi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (y < mesh.y[yi]) { yi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        while (!(mesh.z[zi] <= z && z <= mesh.z[zi + 1])) {
            if (z > mesh.z[zi + 1]) { zi++; } // Если точка правее текущего отрезка, то сместить отрезок вправо
            if (z < mesh.z[zi]) { zi--; } // Если точка левее текущего отрезка, то сместить отрезок влево 
        }
        int n_el = (mesh.x.size() - 1) * (mesh.y.size() - 1) * zi + (mesh.x.size() - 1) * yi + xi; // Номер элемента в котором содержится точка
        double hx = mesh.x[xi + 1] - mesh.x[xi], hy = mesh.y[yi + 1] - mesh.y[yi], hz = mesh.z[zi + 1] - mesh.z[zi];
        double X1 = (mesh.x[xi + 1] - x) / hx, X2 = (x - mesh.x[xi]) / hx;
        double Y1 = (mesh.y[yi + 1] - y) / hy, Y2 = (y - mesh.y[yi]) / hy;
        double Z1 = (mesh.z[zi + 1] - z) / hz, Z2 = (z - mesh.z[zi]) / hz;
        vector<double> psi = {
            X1 * Y1 * Z1,X2 * Y1 * Z1,
            X1 * Y2 * Z1,X2 * Y2 * Z1,
            X1 * Y1 * Z2,X2 * Y1 * Z2,
            X1 * Y2 * Z2,X2 * Y2 * Z2 };
        for (int i = 0; i < FEM_XYZ_NODES_NUM; i++) {
            res += psi[i] * q[ss.nvtr[n_el][i] - 1];
        }
    }
    else {
        res = NAN;
    }
    return res;
}


void FEM::solution_xyz_in_points(std::string filename) {
    ofstream ifs;
    ifs.open(filename, ios::app);
    double x, y, z;
    ifs.setf(ios::scientific | ios::showpos);
    for (int i = 0; i < sreda.points.size(); i++) {
        x = sreda.points[i][0];
        y = sreda.points[i][1];
        z = sreda.points[i][2];
        ifs <<  x << " "<<  y<< " " << z  << " " << solution_xyz_in_point(x, y, z,Vs,mesh_n,ssxyzn) << endl;
    }
    ifs << endl;
    ifs.close();
}
void FEM::solution_rz_in_points(std::string filename) {
    ofstream ifs;
    ifs.open(filename, ios::app);
    double r, z;
    ifs.setf(ios::scientific | ios::showpos);
    for (int i = 0; i < sreda_rz.points.size(); i++) {
        r = sreda_rz.points[i][0];
        z = sreda_rz.points[i][1];
        ifs << r << " " << z  << " " << solution_rz_in_point(r, z) << endl;
    }
    ifs << endl;
    ifs.close();
}

void FEM::save_rz(FEM& fem, string solution_filename) {
    q_rz_to_xyz2(Vns, fem.mesh_n, fem.ssxyzn);
    solution_export(fem.mesh_n, solution_filename, Vns); // Выведем Vn
}

void FEM::save_xyz(FEM& fem, string solution_filename) {
    q_xyz_to_xyz(Vs, mesh_n,ssxyzn,fem.Vs,fem.mesh_n,fem.ssxyzn);
    solution_export(fem.mesh_n, solution_filename, fem.Vs); // Выведем Vn
}


void FEM::solution_export(class MESH& mesh, std::string solution_filename, std::vector<double>& solution) {
    ofstream ifs;
    ifs.open(solution_filename, ios::binary);
    __int32 size;
    // Запишем в файл размеры сеток по x, y и z
    size = mesh.x.size();
    ifs.write((char*)&size, sizeof(size));
    size = mesh.y.size();
    ifs.write((char*)&size, sizeof(size));
    size = mesh.z.size();
    ifs.write((char*)&size, sizeof(size));

    // Запишем координаты узлов
    // Для x
    for (int i = 0; i < mesh.x.size(); i++) {
        ifs.write((char*)&mesh.x[i], sizeof(mesh.x[0]));
    }
    // Для y
    for (int i = 0; i < mesh.y.size(); i++) {
        ifs.write((char*)&mesh.y[i], sizeof(mesh.x[0]));
    }
    // Для z
    for (int i = 0; i < mesh.z.size(); i++) {
        ifs.write((char*)&mesh.z[i], sizeof(mesh.x[0]));
    }

    //Вывод точек решения q
    for (int i = 0; i < solution.size(); i++) {
        ifs.write((char*)&solution[i], sizeof(solution[0]));
    }

    ifs.close();
}
void FEM::solution_rz_export(std::string filename) {
    ofstream ifs;
    ifs.open(filename, ios::binary);
    __int32 size;
    // Запишем в файл размеры сеток по r и z
    size = mesh_rz.r.size();
    ifs.write((char*)&size, sizeof(size));
    size = mesh_rz.z.size();
    ifs.write((char*)&size, sizeof(size));

    // Запишем координаты узлов
    // Для r
    for (int i = 0; i < mesh_rz.r.size(); i++) {
        ifs.write((char*)&mesh_rz.r[i], sizeof(mesh_rz.r[0]));
    }
    // Для z
    for (int i = 0; i < mesh_rz.z.size(); i++) {
        ifs.write((char*)&mesh_rz.z[i], sizeof(mesh_rz.r[0]));
    }

    //Вывод точек решения q; R*Z точек
    for (int i = 0; i < q.size(); i++) {
        ifs.write((char*)&q[i], sizeof(q[0]));
    }

    ifs.close();
}

void FEM::solution_rz_to_xyz2(std::string filename)
{
    ssxyzn.gen_structures(sreda_n, mesh_n);
    q_rz_to_xyz2(Vns, mesh_n, ssxyzn); // Получим решение двумерной задачи для трехмерной сетки sreda+sreda_n
    solution_export(mesh_n,filename,Vns);
}