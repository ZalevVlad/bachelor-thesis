#include <string.h>

#include <iostream>

#include "fem.h"
// #include "fem_n.h"

#define MESH "./meshs/simple_mesh"
#define MESH "./meshs/test_segregation"
#define MESH "./meshs/mesh1"
#define MESH "./meshs/polinom/polinom_lambda_test_6"
#define MESH "./meshs/mesh_2x2"
#define MESH "./meshs/polinom/test_splitting_rz_z"
#define MESH "./meshs/rz_optim/h0"
#define MESH "./meshs/rz_optim/optim"
#define MESH "./meshs/polinom/polinom_lambda_1"

#define MAXITER 100
#define EPS 1e-30
#define RELAX 1.

#define SEG1 1
#define SEG2 1
#define SEG3 1
#define SEG4 1
#define SEG5 1

void fem_new() {
  string filename = "/polinome/polinome_lambda_1";
  FEM fem("./meshs/" + filename);
  // fem.solve(2000, EPS, 1.);
  fem.test_polinome_rz(2000, EPS, 1., "./solutions/tests_polinome_rz/delete");
  // fem.solution_xyz_in_points("solution_in_point.txt");
  // fem.test_polinome_rz(2000, EPS, 1., "./solutions/tests_polinome_rz/delete");
  // fem.solution_rz_in_points("solution_in_point.txt");
  // fem.solve_segregation1(2000, EPS, 1.); // rz->xyz :Р РµС€РµРЅРёРµ РєРѕРЅСЃС‚Р°РЅС‚Р° РЅР° СЌР»РµРјРµРЅС‚Рµ
  // fem.solve_segregation2(2000, EPS, 1.); // rz->xyz: Р РµС€РµРЅРёРµ РІ СѓР·Р»Р°С… СЌР»РµРјРµРЅС‚Р°
  // fem.solve_segregation3(2000, EPS, 1.); // rz->xyz: РўРѕС‡РєРё Р“Р°СѓСЃР°?
  // fem.solve_rz(2000, EPS, 1., "./solutions/rz_optim/optim/q_rz");
}

void test_polinome() {
  string filename = "/polinom/polinom_lambda_1";
  FEM fem("./meshs/" + filename);
  // fem.solve(2000, EPS, 1.);
  fem.test_polinome_rz(2000, EPS, 1., "./solutions/tests_polinome_rz/delete");
}

void vec_find() {
  std::vector<double> u = {1, 2, 3, 10, 20};
  std::vector<double>::iterator it;
  it = std::upper_bound(u.begin(), u.end(), u[0]);
  printf(" bound for %f on position: %d\n", u[0], (int)(it - u.begin()));
  for (int i = 0; i < u.size() - 1; i++) {
    it = std::upper_bound(u.begin(), u.end(), (u[i] + u[i + 1]) / 2.);
    printf(" bound for %f on position: %d\n", (u[i] + u[i + 1]) / 2., (int)(it - u.begin()));
    it = std::upper_bound(u.begin(), u.end(), u[i + 1]);
    printf(" bound for %f on position: %d\n", u[i + 1], (int)(it - u.begin()));
  }

  printf("\nLOWER\n");
  it = std::lower_bound(u.begin(), u.end(), u[0]);
  printf(" bound for %f on position: %d\n", u[0], (int)(it - u.begin()));
  for (int i = 0; i < u.size() - 1; i++) {
    it = std::lower_bound(u.begin(), u.end(), (u[i] + u[i + 1]) / 2.);
    printf(" bound for %f on position: %d\n", (u[i] + u[i + 1]) / 2., (int)(it - u.begin()));
    it = std::lower_bound(u.begin(), u.end(), u[i + 1]);
    printf(" bound for %f on position: %d\n", u[i + 1], (int)(it - u.begin()));
  }
}

void matrix_mul() {
  vector<double> ia, ja, ggl, di, pr, q, x;
  ia = {1, 1, 2, 4};
  ja = {1, 1, 2};
  ggl = {2, 3, 5};
  di = {1, 4, 6};
  q = {1, 2, 3};
  pr = {14, 25, 31};

  int n = di.size();
  x.resize(n);
  for (unsigned int i = 0; i < n; i++) {
    x[i] = di[i] * q[i];
    for (unsigned int k = ia[i] - 1, k1 = ia[i + 1] - 1; k < k1; k++) {
      unsigned int j = ja[k] - 1;
      x[i] += ggl[k] * q[j];
      x[j] += ggl[k] * q[i];
    }
  }
}

template <typename T>
std::string toString(T val) {
  std::stringstream s;
  s << val;  // РїСЂРµРѕР±СЂР°Р·РѕРІР°РЅРёРµ РІ СЃС‚СЂРѕРєСѓ
  return s.str();
}

void test_polinome_rz_z_splitting() {
  for (int i = 0; i < 9; i++) {
    std::cout << "TEST SPLITTING #" << i << endl;
    char buf[20];
    sprintf_s(buf, "/%d", i);
    string filename = MESH + (string)buf;
    FEM fem(filename);
    fem.test_polinome_rz(2000, EPS, 1., "./solutions/tests_polinome_rz/delete");
    fem.solution_rz_in_points("solution_in_point.txt");
  }
}

void test_rz_optim() {
  std::cout << "\nTEST RZ_OPTIM" << endl;
  string filename = "./meshs/rz_optim/optim";
  FEM fem(filename);
  fem.solve_rz(2000, EPS, 1., (string) "./solutions/rz_optim/optim/q_rz");
  // fem.solution_rz_in_points("solution_in_point.txt");
}

void test_xyz_optim() {
  std::cout << "\nTEST XYZ_OPTIM" << endl;
  string filename = "./meshs/xyz_optim";
  FEM fem(filename);
  fem.solve(2000, EPS, 1., (string) "./solutions/xyz_optim/q_xyz");
  // fem.solution_rz_in_points("solution_in_point.txt");
}

void test_rz_optim_h0() {
  for (int i = 1; i < 6; i++) {
    std::cout << "\nTEST RZ_OPTIM H0 #" << i << endl;
    char buf[20];
    sprintf_s(buf, "/1e-%d", i);
    string filename = MESH + (string)buf;
    FEM fem(filename);
    fem.solve_rz(2000, EPS, 1., (string) "./solutions/rz_optim/h0" + (string)buf);
    // fem.solution_rz_in_points("solution_in_point.txt");
  }
}

void test_rz_optim_k() {
  for (int i = 1; i < 11; i++) {
    std::cout << "\nTEST RZ_OPTIM K #" << i << endl;
    char buf[20];
    sprintf_s(buf, "/1_%02d", i);
    string filename = "./meshs/rz_optim/k" + (string)buf;
    FEM fem(filename);
    fem.solve_rz(2000, EPS, 1., (string) "./solutions/rz_optim/k" + (string)buf);
  }
}

void segregation_test() {
  string folder = "test";
  std::cout << "\nSEGREGATION TEST 1\n";
  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  fem.solution_rz_to_xyz2("./solutions/segregation/" + folder + "/rz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
}

void segregation_analitycal() {
  string folder = "analytical";
  std::cout << "\nSEGREGATION ANALYTICAL 1\n";
  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  fem.solution_rz_to_xyz2("./solutions/segregation/" + folder + "/rz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
}

void i2() {
  string folder = "i2";
  std::cout << "\nTEST 2\n";
  FEM femrz("./meshs/segregation/" + folder + "/rz");
  femrz.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz");
  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  femrz.save_rz(fem, "./solutions/segregation/" + folder + "/xyz_accur");
}

void i3() {
  string folder = "i3";
  std::cout << "\nTEST 3\n";
  FEM femrz("./meshs/segregation/" + folder + "/rz");
  femrz.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz");
  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  femrz.save_rz(fem, "./solutions/segregation/" + folder + "/xyz_accur");
}

void i4() {
  string folder = "i4";
  std::cout << "\nTEST 4\n";
  /*FEM_N fem_big("./meshs/segregation/" + folder + "/xyz_big");
  fem_big.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  fem_big.solve_segregation3(2000, EPS, 1., (string)"./solutions/segregation/" + folder + "/xyz_big");
  fem_big.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");*/

  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  // fem_big.save_xyz(fem, "./solutions/segregation/" + folder + "/xyz_accur");
}

void i5() {
  string folder = "i5";
  std::cout << "\nTEST 5\n";

  FEM fem_big("./meshs/segregation/" + folder + "/xyz_big");
  fem_big.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  fem_big.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz_big");
  fem_big.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");

  FEM fem("./meshs/segregation/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation/" + folder + "/xyz/solutions");
  }
  fem_big.save_xyz(fem, "./solutions/segregation/" + folder + "/xyz_accur");
}

void solution_trasnition_test() {
  std::cout << "\nSEGREGATION ANALYTICAL 1\n";
  FEM fem("./meshs/segregation/analytical/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/analytical/rz_for_xyz");
  // fem.solution_rz_to_xyz1("./solutions/segregation/analytical/xyz1");
  fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation/analytical/xyz1");
}

void segregation_3_test() {
  std::cout << "\nSEGREGATION 3 TEST\n";
  FEM fem("./meshs/segregation/seg3test/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation/seg3test/rz_for_xyz");
  fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation/seg3test/xyz3");
  fem.solution_xyz_in_points("./meshs/segregation/seg3test/xyz/solutions");
}

void segregation_analitycal_2() {
  string folder = "analytical";
  std::cout << "\nSEGREGATION ANALYTICAL 1\n";
  FEM fem("./meshs/segregation2/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/segregation2/" + folder + "/rz_for_xyz");
  fem.solution_rz_to_xyz2("./solutions/segregation2/" + folder + "/rz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/segregation2/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/segregation2/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/segregation2/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/segregation2/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/segregation2/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/segregation2/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/segregation2/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/segregation2/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/segregation2/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/segregation2/" + folder + "/xyz/solutions");
  }
}

void i2_2() {
  string folder = "segregation2/i2";
  std::cout << "\nTEST 2\n";
  FEM femrz("./meshs/" + folder + "/rz");
  femrz.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz");
  FEM fem("./meshs/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  femrz.save_rz(fem, "./solutions/" + folder + "/xyz_accur");
}

void i3_2() {
  string folder = "segregation2/i3";
  std::cout << "\nTEST 3\n";
  FEM femrz("./meshs/" + folder + "/rz");
  femrz.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz");
  FEM fem("./meshs/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  femrz.save_rz(fem, "./solutions/" + folder + "/xyz_accur");
}

void i4_2() {
  string folder = "segregation2/i4";
  std::cout << "\nTEST 4\n";
  FEM fem_big("./meshs/" + folder + "/xyz_big");
  fem_big.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  fem_big.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz_big");
  fem_big.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");

  FEM fem("./meshs/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  fem_big.save_xyz(fem, "./solutions/" + folder + "/xyz_accur");
}

void i5_2() {
  string folder = "segregation2/i5";
  std::cout << "\nTEST 5\n";

  FEM fem_big("./meshs/" + folder + "/xyz_big");
  fem_big.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  fem_big.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz_big");
  fem_big.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");

  FEM fem("./meshs/" + folder + "/xyz");
  fem.solve_rz(2000, EPS, 1., "./solutions/" + folder + "/rz_for_xyz");
  if (SEG1) {
    fem.solve_segregation1(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz1");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG2) {
    fem.solve_segregation2(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz2");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG3) {
    fem.solve_segregation3(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz3");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG4) {
    fem.solve_segregation4(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz4");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  if (SEG5) {
    fem.solve_segregation5(2000, EPS, 1., (string) "./solutions/" + folder + "/xyz5");
    fem.solution_xyz_in_points("./meshs/" + folder + "/xyz/solutions");
  }
  fem_big.save_xyz(fem, "./solutions/" + folder + "/xyz_accur");
}

int main() {
  // vec_find();
  // fem_new();
  // test_polinome_rz_z_splitting();
  // test_rz_optim();
  // test_xyz_optim(); // РћРґРЅРѕСЂРѕРґРЅР°СЏ СЃСЂРµРґР°
  // test_rz_optim_h0();
  // test_rz_optim_k();
  test_polinome();

  // segregation_test();
  // segregation_analitycal();
  // i2();
  // i3();
  // i4();
  // i5();
  // solution_trasnition_test();
  // segregation_3_test();
  // segregation_analitycal_2();
  // i2_2();
  // i3_2();
  // i4_2();
  // i5_2();
}