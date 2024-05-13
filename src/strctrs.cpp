#include "STRCTRS.h"

// Генерирует номера узлов элементов для трехмерной сетки
void STRCTRS::gen_nvtr(class MESH& mesh) {
  double kel = mesh.kel;
  nvtr.resize(mesh.kel);
  for (int i = 0; i < kel; i++) {
    nvtr[i].resize(8);
  }

  int x, y, z;
  // Количество элементов в слое x
  int kx = (mesh.x.size() - 1);
  // Количество элементов в слое xy
  int kxy = (mesh.x.size() - 1) * (mesh.y.size() - 1);
  for (int i = 0; i < kel; i++) {
    x = i % kx + 1;
    y = (i % kxy) / (kx) + 1;
    z = i / kxy + 1;

    nvtr[i][0] = x + (y - 1) * (kx + 1) + (z - 1) * mesh.x.size() * mesh.y.size();
    nvtr[i][1] = nvtr[i][0] + 1;
    nvtr[i][2] = nvtr[i][0] + mesh.x.size();
    nvtr[i][3] = nvtr[i][2] + 1;
    nvtr[i][4] = nvtr[i][0] + mesh.x.size() * mesh.y.size();
    nvtr[i][5] = nvtr[i][4] + 1;
    nvtr[i][6] = nvtr[i][4] + mesh.x.size();
    nvtr[i][7] = nvtr[i][6] + 1;
  }
}

// Генерирует номера узлов элементов для двумерной сетки
void STRCTRS::gen_nvtr_rz(class MESH_RZ& mesh) {
  double kel = mesh.kel;
  nvtr.resize(mesh.kel);
  for (int i = 0; i < kel; i++) {
    nvtr[i].resize(4);
  }

  int n;
  for (int i = 0; i < kel; i++) {
    n = i / (mesh.r.size() - 1) * (mesh.r.size()) + i % (mesh.r.size() - 1) + 1;

    nvtr[i][0] = n;
    nvtr[i][1] = n + 1;
    nvtr[i][2] = n + mesh.r.size();
    nvtr[i][3] = n + mesh.r.size() + 1;
  }
}

// Генерирует координаты узлов трехмерной сетки
void STRCTRS::gen_coord(class MESH& mesh) {
  coord.resize(mesh.kuzlov);
  int xs = mesh.x.size();
  int ys = mesh.y.size();
  int zs = mesh.z.size();
  int x, y, z;
  for (int i = 0; i < mesh.kuzlov; i++) {
    x = i % xs;
    y = (i % (xs * ys)) / xs;
    z = i / (xs * ys);
    coord[i].resize(3);
    coord[i] = {mesh.x[x], mesh.y[y], mesh.z[z]};
  }
}

// Генерирует коориднаты узлов двумерной сетки
void STRCTRS::gen_coord_rz(class MESH_RZ& mesh) {
  coord.resize(mesh.kuzlov);
  int rs = mesh.r.size();
  int zs = mesh.z.size();
  int r, h;
  for (int i = 0; i < mesh.kuzlov; i++) {
    r = i % rs;
    h = i / rs;
    coord[i].resize(2);
    coord[i] = {mesh.r[r], mesh.z[h]};
  }
}

// По кооринатам определит номера материала
// Обход обратный: будет получен материал аномального поля, если он не существует, то нормального
int STRCTRS::get_mat_from_sreda_reverse(class SREDA& sreda, double x, double y, double z) {
  // Будем обходить подобласти в обратном порядке (учитывая что в файле среда подоблости заданны вложением по порядку)
  for (int i = sreda.elms.size() - 1; i >= 0; i--) {
    /*printf("%lf %lf %lf\n", sreda.elms[i][X0_COORD], x, sreda.elms[i][X1_COORD]);
    printf("%lf %lf %lf\n", sreda.elms[i][Y0_COORD], y, sreda.elms[i][Y1_COORD]);
    printf("%lf %lf %lf\n", sreda.elms[i][Z0_COORD], z, sreda.elms[i][Z1_COORD]);*/
    if (sreda.elms[i][X0_COORD] < x && x < sreda.elms[i][X1_COORD] &&
        sreda.elms[i][Y0_COORD] < y && y < sreda.elms[i][Y1_COORD] &&
        sreda.elms[i][Z0_COORD] < z && z < sreda.elms[i][Z1_COORD]) {
      return sreda.elms[i][MAT_N];
    }
  }
  return 0;
}

// По кооринатам определит номера материала
// Обход прямой: в первую очередь будет получен материал нормального поля
int STRCTRS::get_mat_from_sreda_direct(class SREDA& sreda, double x, double y, double z) {
  // Будем обходить подобласти в обратном порядке (учитывая что в файле среда подоблости заданны вложением по порядку)
  for (int i = 0; i < sreda.elms.size(); i++) {
    /*printf("%lf %lf %lf\n", sreda.elms[i][X0_COORD], x, sreda.elms[i][X1_COORD]);
    printf("%lf %lf %lf\n", sreda.elms[i][Y0_COORD], y, sreda.elms[i][Y1_COORD]);
    printf("%lf %lf %lf\n", sreda.elms[i][Z0_COORD], z, sreda.elms[i][Z1_COORD]);*/
    if (sreda.elms[i][X0_COORD] < x && x < sreda.elms[i][X1_COORD] &&
        sreda.elms[i][Y0_COORD] < y && y < sreda.elms[i][Y1_COORD] &&
        sreda.elms[i][Z0_COORD] < z && z < sreda.elms[i][Z1_COORD]) {
      return sreda.elms[i][MAT_N];
    }
  }
  return 0;
}

// По кооринатам определит номер материала для двумерной задачи
int STRCTRS::get_mat_from_sreda_rz(class SREDA& sreda, double z) {
  // Будем обходить подобласти в прямом порядке (учитывая что в файле среда подоблости заданны вложением по порядку)
  // Считается, что для двумерной задачи элементы, которые не располагаются от начала до конца границы являются аномальными
  for (int i = 0; i < sreda.elms.size(); i++) {
    if (sreda.elms[i][Z0_COORD] < z && z < sreda.elms[i][Z1_COORD]) {
      return sreda.elms[i][MAT_N];
    }
  }
  return 0;
}

// Генерирует номера материалов параллелепипедов для трехмерной задачи
void STRCTRS::gen_nvkat2d(class SREDA& sreda, class MESH& mesh) {
  nvkat2d.resize(mesh.kel);
  double x, y, z;
  for (int i = 0; i < mesh.kel; i++) {
    // Координаты центров элементов
    x = (coord[nvtr[i][0] - 1][0] + coord[nvtr[i][1] - 1][0]) / 2.;
    y = (coord[nvtr[i][0] - 1][1] + coord[nvtr[i][2] - 1][1]) / 2.;
    z = (coord[nvtr[i][0] - 1][2] + coord[nvtr[i][4] - 1][2]) / 2.;
    nvkat2d[i] = get_mat_from_sreda_reverse(sreda, x, y, z);
  }
}

void STRCTRS::gen_nvkat2dr(class SREDA& sreda, class MESH& mesh) {
  nvkat2dr.resize(mesh.kel);
  double x, y, z;
  for (int i = 0; i < mesh.kel; i++) {
    // Координаты центров элементов
    x = (coord[nvtr[i][0] - 1][0] + coord[nvtr[i][1] - 1][0]) / 2.;
    y = (coord[nvtr[i][0] - 1][1] + coord[nvtr[i][2] - 1][1]) / 2.;
    z = (coord[nvtr[i][0] - 1][2] + coord[nvtr[i][4] - 1][2]) / 2.;
    nvkat2dr[i] = get_mat_from_sreda_direct(sreda, x, y, z);
  }
}

// Генерирует номера материалов прямоугольников для двумерной задачи
void STRCTRS::gen_nvkat2d_rz(class MESH_RZ& mesh, class SREDA& sreda) {
  nvkat2d.resize(mesh.kel);
  double r, z;
  for (int i = 0; i < mesh.kel; i++) {
    // Координаты центров элементов
    // r = (coord[nvtr[i][0] - 1][0] + coord[nvtr[i][1] - 1][0]) / 2.;
    z = (coord[nvtr[i][0] - 1][1] + coord[nvtr[i][2] - 1][1]) / 2.;
    nvkat2d[i] = get_mat_from_sreda_rz(sreda, z);
  }
}

// Генерирует номера узлов с первыми нулевыми краевыми для трехмерной задачи
void STRCTRS::gen_l1(class MESH& mesh, vector<double>& edge_conditions) {
  l1.resize(0);
  set<int> buff;
  for (int i = 0; i < mesh.kuzlov; i++) {
    // Если указано, что на грани x = x0 заданы первые краевые,
    // то если координата текущей точки совпадает с x0, тогда добавить ее номер.
    // С остальными гранями аналогично
    if (edge_conditions[L1_X0]) {
      if (coord[i][0] == mesh.x[0]) {
        buff.insert(i + 1);
      }
    }
    if (edge_conditions[L1_X1]) {
      if (coord[i][0] == mesh.x[mesh.x.size() - 1]) {
        buff.insert(i + 1);
      }
    }

    if (edge_conditions[L1_Y0]) {
      if (coord[i][1] == mesh.y[0]) {
        buff.insert(i + 1);
      }
    }
    if (edge_conditions[L1_Y1]) {
      if (coord[i][1] == mesh.y[mesh.y.size() - 1]) {
        buff.insert(i + 1);
      }
    }

    if (edge_conditions[L1_Z0]) {
      if (coord[i][2] == mesh.z[0]) {
        buff.insert(i + 1);
      }
    }
    if (edge_conditions[L1_Z1]) {
      if (coord[i][2] == mesh.z[mesh.z.size() - 1]) {
        buff.insert(i + 1);
      }
    }
  }
  l1.insert(l1.begin(), buff.begin(), buff.end());
}
// Генерирует номера узлов с первыми нулевыми краевыми для двумерной задачи
void STRCTRS::gen_l1_rz(class MESH_RZ& mesh, vector<double>& edge_conditions) {
  l1.resize(0);
  set<int> buff;
  /* Краевые условия в двумерной задаче заданы следующим образом

     S2=0
   ----------
   |        |
s2=0 |        | S1 = 0
   |        |
   ----------
    S1 = 0
   */
  for (int i = 0; i < mesh.kuzlov; i++) {
    if (coord[i][0] == mesh.r[mesh.r.size() - 1]) {
      buff.insert(i + 1);
    } else if (coord[i][1] == mesh.z[mesh.z.size() - 1]) {
      buff.insert(i + 1);
    }
  }
  l1.insert(l1.begin(), buff.begin(), buff.end());
}

void STRCTRS::gen_l1_p(class MESH& mesh, vector<double>& edge_conditions) {
  l1.resize(0);
  set<int> buff;
  for (int i = 0; i < mesh.kuzlov; i++) {
    // Для теста на полином зададим первые краевые на всех границах
    if (coord[i][0] == mesh.x[0]) {
      buff.insert(i + 1);
    } else if (coord[i][0] == mesh.x[mesh.x.size() - 1]) {
      buff.insert(i + 1);
    } else if (coord[i][1] == mesh.y[0]) {
      buff.insert(i + 1);
    } else if (coord[i][1] == mesh.y[mesh.y.size() - 1]) {
      buff.insert(i + 1);
    } else if (coord[i][2] == mesh.z[0]) {
      buff.insert(i + 1);
    } else if (coord[i][2] == mesh.z[mesh.z.size() - 1]) {
      buff.insert(i + 1);
    }
  }
  l1.insert(l1.begin(), buff.begin(), buff.end());
}
void STRCTRS::gen_l1_rz_p(class MESH_RZ& mesh, vector<double>& edge_conditions) {
  l1.resize(0);
  set<int> buff;
  /* Для теста на полином зададим первые краевые на всех границах

       S1
   ----------
   |        |
  S1 |        | S1
   |        |
   ----------
     S1
   */
  for (int i = 0; i < mesh.kuzlov; i++) {
    if (coord[i][0] == mesh.r[mesh.r.size() - 1]) {
      buff.insert(i + 1);
    } else if (coord[i][0] == mesh.r[0]) {
      buff.insert(i + 1);
    } else if (coord[i][1] == mesh.z[mesh.z.size() - 1]) {
      buff.insert(i + 1);
    } else if (coord[i][1] == mesh.z[0]) {
      buff.insert(i + 1);
    }
  }
  l1.insert(l1.begin(), buff.begin(), buff.end());
}

// nvtr : 8 * int номера вершин прямоугольников
// nvkat2d : 1 * int номера матреиала прямоугольников
// coord : 3 * double координаты вершин
// l1 : 1 * int (номера вершин с первым нулевым краевым условием)
void STRCTRS::gen_structures(SREDA& sreda, MESH& mesh) {
  gen_nvtr(mesh);
  gen_coord(mesh);
  gen_nvkat2d(sreda, mesh);
  gen_nvkat2dr(sreda, mesh);
  gen_l1(mesh, sreda.edge_conditions);
}
// nvtr : 4 * int номера вершин прямоугольников
// nvkat2d : 1 * int номера матреиала прямоугольников
// coord : 2 * double координаты вершин
// l1 : 1 * int (номера вершин с первым нулевым краевым условием)
void STRCTRS::gen_structures_rz(SREDA& sreda, MESH_RZ& mesh) {
  gen_nvtr_rz(mesh);
  gen_coord_rz(mesh);
  gen_nvkat2d_rz(mesh, sreda);
  gen_l1_rz(mesh, sreda.edge_conditions);
}

// Тест на полиномы
void STRCTRS::gen_structures_p(SREDA& sreda, MESH& mesh) {
  gen_nvtr(mesh);
  gen_coord(mesh);
  gen_nvkat2d(sreda, mesh);
  gen_l1_p(mesh, sreda.edge_conditions);
}
void STRCTRS::gen_structures_rz_p(SREDA& sreda, MESH_RZ& mesh) {
  gen_nvtr_rz(mesh);
  gen_coord_rz(mesh);
  gen_nvkat2d_rz(mesh, sreda);
  gen_l1_rz_p(mesh, sreda.edge_conditions);
}