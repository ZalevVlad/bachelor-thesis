#ifndef FEM_RZ_H_
#define FEM_RZ_H_

#include <algorithm>
#include <time.h>

#include "sreda.h"
#include "STRCTRS.h"
#include "cgm.h"


// Количество узлов в элементе
#define FEM_XYZ_NODES_NUM 8
#define FEM_RZ_NODES_NUM 4

class FEM
{
	void print_matrix();
	int find_nel(double x, double y, double z, class MESH& mesh, class STRCTRS& ss, std::vector<double>& q);
	SREDA sreda;
	MESH mesh;

	SREDA sreda_n;
	MESH mesh_n;

	SREDA_RZ sreda_rz;
	MESH_RZ mesh_rz;

	STRCTRS ssxyz, ssxyzn,ssrz;
	// СЛАУ
	// ggl - нижний треугольник
	// di - диагональные элементы
	std::vector<double>ggl, ggu, di, pr, qrz, q;
	std::vector<long> ia, ja;
	std::vector < std::set<long>> A; // Профиль ВСЕЙ матрицы, используется для применения первых краевых
	std::vector<double> Vns; // Решение двумерной задачи в узлах общей сетки
	std::vector<double> Vnb; // Решение двумерной задачи в узлах sreda
	std::vector<double> Va;  // Решение задачи на аномальное поле в узлах sreda 
	std::vector<double> Vas; // Решение задачи на аномальное поле в узлах общей сетки
	std::vector<double> Vs; // Сумма решений Vas и Vns

	int xi, yi, zi; // Переменные поиска элемента для функции xyz_to_xyz
	// Статические переменные не используется для корректной работы последовательного расчета 
public:
	// filename - 
	FEM(string filename);
	void solve(int max_iter, double eps, double relax, std::string solution_filename);
	void solve_segregation1(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz :Решение константа на элементе
	void solve_segregation2(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: Решение в узлах элемента
	void solve_segregation3(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: Точки Гауса?
	void solve_segregation4(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: Точки Гауса?
	void solve_segregation5(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: Точки Гауса?
	void solve_rz(int max_iter, double eps, double relax, std::string solution_filename);

	void solution_xyz_in_points(std::string filename);
	void solution_rz_in_points(std::string filename);

	void save_rz(FEM& fem, string filename);
	void save_xyz(FEM& fem, string filename);


private:
	void build_matrix_profile();
	void build_matrix_profile_rz();
	void slae_init();
	
	void add_G(); // Строит матрицу жесткости для задачи xyz без выделения поля
	void add_b(); // Правая часть для задачи без выделения поля
	void add_G_sigma_diff(); // Строит матрицу с sigma = sigma(normal) - sigma
	void add_b_segregation_1(); // Считает правую часть для первого способа выделения поля
	void add_b_segregation_2(); // Считает правую часть для второго способа выделения поля
	void add_b_segregation_3(); // Считает правую часть для третьего способа выделения поля
	void add_b_segregation_4(); // Считает правую часть для четвертого способа выделения поля
	void add_b_segregation_5(); // Считает правую часть для пятого способа выделения поля
	void add_b_segregation_3(int NOMEL); // Считает правую часть для третьего способа выделения поля
	void add_b_segregation_4(int NOMEL); // Считает правую часть для четвертого способа выделения поля
	void add_b_segregation_5(int NOMEL); // Считает правую часть для пятого способа выделения поля

	void add_G_b_rz(); // Строит глобальные матрицу жесткости и вектор правой части для задачи rz
	void edge_cond_1(); // Первые краевые для задачи xyz
	void edge_cond_1_rz(); // Первые краевые для задачи rz

	// Выведет заполненность матрицы ненулевыми элементами
	void print_matrix_plenum(); 

	// Номер элемента в который попал узел
	int n_el_xyz(double x, double y, double z);

	// Перевести решение двумерной задачи в решение трехмерной
	// 1: Формирует вектор с длиной равной КОЛИЧЕСТВУ УЗЛОВ с решением двумерной задачи в КАЖДОМ ЭЛЕМЕНТЕ
	void q_rz_to_xyz1(vector<double>& ans, class MESH& mesh, class STRCTRS& ss);
	// 2: Формирует вектор с длиной равной КОЛИЧЕСТВУ УЗЛОВ с решением двумерной задачи в КАЖДОМ УЗЛЕ
	void q_rz_to_xyz2(vector<double>& ans, class MESH& mesh, class STRCTRS& ss);


	// Перевести решение трехмерной задачи с одной сетки в другую
	// q1 - исходный вектор решений
	// mesh1 - сетка исходного решения
	// 
	// q2 - получаемый вектор решений
	// mesh2 - сетка для получаемого решения
	// ss2 - струкутры получаемого решения
	void q_xyz_to_xyz(vector<double>& q1, class MESH& mesh1,class STRCTRS& ss1, vector<double>& q2, class MESH& mesh2, class STRCTRS& ss2);

private:
	void add_to_sparse(std::vector<long>& ia, std::vector<long>& ja, std::vector<double>& ggl, long str, long col, double x);
	void replace_in_sparse(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x);
	double replace_in_sparse_r(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x);
	// Сложение двух векторов x = a + b
	void vec_vec_sum(std::vector<double>& a, std::vector<double>& b, std::vector<double>& x);
	// Умножит матрицу на вектор q, результат будет помещен в вектор x
	void matrix_vector_mul(std::vector<double>& q, std::vector<double>& x);

	void solve_cgm(vector<double>& q, int max_iter, double eps, double relax);


	void local_G(vector<vector<double>>& G, double h_x, double h_y, double h_z, double lambda);
	void local_G_rz(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i);
	void local_G_rzn(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i); //Удалить
	// -------Функции тестирования на полиномах---------
public:
	void test_polinome(int max_iter, double eps, double relax, std::string solution_filename);
	void test_polinome_rz(int max_iter, double eps, double relax, std::string solution_filename);
	
	double solution_xyz_in_point(double x, double y, double z); // Решение трехмерной задачи в точке
	double solution_xyz_in_point(double x, double y, double z,
		std::vector<double>& q, class MESH& mesh, class STRCTRS& ss); // Решение трехмерной задачи для сетки mesh с решением q
	double solution_rz_in_point(double r, double z); // Выдает решение двумерной задачи в точке
	void grad_rz_in_point(double x, double y, double z, std::vector<double>& grad); // Значение градиента по декартовым координатам для цилиндрической задачи
	void grad_rz_in_point2(double x, double y, double z, std::vector<double>& grad); // Значение градиента по декартовым координатам для цилиндрической задачи
	void solution_export(class MESH& mesh, std::string solution_filename, std::vector<double>& solution);
	void solution_rz_export(std::string filename);

	void solution_rz_to_xyz2(std::string filename); // Переведет точки решения двумерной задачи в трехмерную
private:
	void add_G_b_p(); // Считает глобальную матрицу жесткости и правую часть
	void add_G_b_rz_p();
	void edge_cond_1_p();
	void edge_cond_1_rz_p();
	void print_solution_p_miss();
	void print_solution_rz_p_miss(); // Норма вектора относительной погрешности решения

	// u - искомая функция дифференциального уравнения
	// f - функция правой части дифференциального уравнения
	double u(double x, double y, double z);
	double f(double x, double y, double z);
	double u_rz(double r, double z);
	double f_rz(double r, double z);

	// Cчитает локальную правую часть
	void local_b_p(vector<double>& b, double h_x, double h_y, double h_z, int i);
	void local_b_rz_p(vector<double>& b, double h_r, double h_z, int i); // тест на полином
	void local_b_rz_pn(vector<double>& b, double h_r, double h_z, int i); // удалить
	//---------------------------------------------------
};

#endif // FEM_RZ_H_
