#ifndef MESH_H_
#define MESH_H_

#include <list>
#include <set>
#include <vector>

// MAT_N - номер материала
//
enum { X0_COORD,
       X1_COORD,
       Y0_COORD,
       Y1_COORD,
       Z0_COORD,
       Z1_COORD,
       MAT_N,
       ANOMAL };

// Описание массива источников
enum { SOURCE_X,
       SOURCE_Y,
       SOURCE_Z,
       SOURCE_POW,
       SOURCE_N_UZLA };

// Описание массива с первыми краевыми
enum { L1_X0,
       L1_X1,
       L1_Y0,
       L1_Y1,
       L1_Z0,
       L1_Z1 };

using namespace std;

int get_mat_from_sreda_direct(class SREDA& sreda, double x, double y, double z);

class MESH_RZ {
 public:
  int kuzlov;
  int kel;
  std::vector<double> r;
  std::vector<double> z;
  void gen_mesh(class SREDA_RZ& sreda_rz, class SREDA& sreda);
};

struct MESH {
 public:
  int kuzlov;
  int kel;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> r;
  std::vector<double> h;

  /*vector<vector<int>> nvtr;
  vector<int> nvkat2d;
  vector<vector<double>> rz;
  vector<int> l1;
  vector<double> sigma;*/

  void gen_mesh(class SREDA& sreda);
  void gen_structures(class FEM& fem, class SREDA& sreda);

  void gen_mesh_rh(class SREDA& sreda);
  void gen_structures_rh(class FEM& fem, class SREDA& sreda);

  void add_mesh(class MESH& mesh);
};

#endif MESH_H_