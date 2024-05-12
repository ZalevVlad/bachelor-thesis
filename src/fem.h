#ifndef FEM_RZ_H_
#define FEM_RZ_H_

#include <algorithm>
#include <time.h>

#include "sreda.h"
#include "STRCTRS.h"
#include "cgm.h"


// ���������� ����� � ��������
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
	// ����
	// ggl - ������ �����������
	// di - ������������ ��������
	std::vector<double>ggl, ggu, di, pr, qrz, q;
	std::vector<long> ia, ja;
	std::vector < std::set<long>> A; // ������� ���� �������, ������������ ��� ���������� ������ �������
	std::vector<double> Vns; // ������� ��������� ������ � ����� ����� �����
	std::vector<double> Vnb; // ������� ��������� ������ � ����� sreda
	std::vector<double> Va;  // ������� ������ �� ���������� ���� � ����� sreda 
	std::vector<double> Vas; // ������� ������ �� ���������� ���� � ����� ����� �����
	std::vector<double> Vs; // ����� ������� Vas � Vns

	int xi, yi, zi; // ���������� ������ �������� ��� ������� xyz_to_xyz
	// ����������� ���������� �� ������������ ��� ���������� ������ ����������������� ������� 
public:
	// filename - 
	FEM(string filename);
	void solve(int max_iter, double eps, double relax, std::string solution_filename);
	void solve_segregation1(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz :������� ��������� �� ��������
	void solve_segregation2(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: ������� � ����� ��������
	void solve_segregation3(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: ����� �����?
	void solve_segregation4(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: ����� �����?
	void solve_segregation5(int max_iter, double eps, double relax, std::string solution_filename); // rz->xyz: ����� �����?
	void solve_rz(int max_iter, double eps, double relax, std::string solution_filename);

	void solution_xyz_in_points(std::string filename);
	void solution_rz_in_points(std::string filename);

	void save_rz(FEM& fem, string filename);
	void save_xyz(FEM& fem, string filename);


private:
	void build_matrix_profile();
	void build_matrix_profile_rz();
	void slae_init();
	
	void add_G(); // ������ ������� ��������� ��� ������ xyz ��� ��������� ����
	void add_b(); // ������ ����� ��� ������ ��� ��������� ����
	void add_G_sigma_diff(); // ������ ������� � sigma = sigma(normal) - sigma
	void add_b_segregation_1(); // ������� ������ ����� ��� ������� ������� ��������� ����
	void add_b_segregation_2(); // ������� ������ ����� ��� ������� ������� ��������� ����
	void add_b_segregation_3(); // ������� ������ ����� ��� �������� ������� ��������� ����
	void add_b_segregation_4(); // ������� ������ ����� ��� ���������� ������� ��������� ����
	void add_b_segregation_5(); // ������� ������ ����� ��� ������ ������� ��������� ����
	void add_b_segregation_3(int NOMEL); // ������� ������ ����� ��� �������� ������� ��������� ����
	void add_b_segregation_4(int NOMEL); // ������� ������ ����� ��� ���������� ������� ��������� ����
	void add_b_segregation_5(int NOMEL); // ������� ������ ����� ��� ������ ������� ��������� ����

	void add_G_b_rz(); // ������ ���������� ������� ��������� � ������ ������ ����� ��� ������ rz
	void edge_cond_1(); // ������ ������� ��� ������ xyz
	void edge_cond_1_rz(); // ������ ������� ��� ������ rz

	// ������� ������������� ������� ���������� ����������
	void print_matrix_plenum(); 

	// ����� �������� � ������� ����� ����
	int n_el_xyz(double x, double y, double z);

	// ��������� ������� ��������� ������ � ������� ����������
	// 1: ��������� ������ � ������ ������ ���������� ����� � �������� ��������� ������ � ������ ��������
	void q_rz_to_xyz1(vector<double>& ans, class MESH& mesh, class STRCTRS& ss);
	// 2: ��������� ������ � ������ ������ ���������� ����� � �������� ��������� ������ � ������ ����
	void q_rz_to_xyz2(vector<double>& ans, class MESH& mesh, class STRCTRS& ss);


	// ��������� ������� ���������� ������ � ����� ����� � ������
	// q1 - �������� ������ �������
	// mesh1 - ����� ��������� �������
	// 
	// q2 - ���������� ������ �������
	// mesh2 - ����� ��� ����������� �������
	// ss2 - ��������� ����������� �������
	void q_xyz_to_xyz(vector<double>& q1, class MESH& mesh1,class STRCTRS& ss1, vector<double>& q2, class MESH& mesh2, class STRCTRS& ss2);

private:
	void add_to_sparse(std::vector<long>& ia, std::vector<long>& ja, std::vector<double>& ggl, long str, long col, double x);
	void replace_in_sparse(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x);
	double replace_in_sparse_r(vector<long>& ia, vector<long>& ja, vector<double>& ggl, long str, long col, double x);
	// �������� ���� �������� x = a + b
	void vec_vec_sum(std::vector<double>& a, std::vector<double>& b, std::vector<double>& x);
	// ������� ������� �� ������ q, ��������� ����� ������� � ������ x
	void matrix_vector_mul(std::vector<double>& q, std::vector<double>& x);

	void solve_cgm(vector<double>& q, int max_iter, double eps, double relax);


	void local_G(vector<vector<double>>& G, double h_x, double h_y, double h_z, double lambda);
	void local_G_rz(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i);
	void local_G_rzn(vector<vector<double>>& G, double h_r, double h_z, double lambda, int i); //�������
	// -------������� ������������ �� ���������---------
public:
	void test_polinome(int max_iter, double eps, double relax, std::string solution_filename);
	void test_polinome_rz(int max_iter, double eps, double relax, std::string solution_filename);
	
	double solution_xyz_in_point(double x, double y, double z); // ������� ���������� ������ � �����
	double solution_xyz_in_point(double x, double y, double z,
		std::vector<double>& q, class MESH& mesh, class STRCTRS& ss); // ������� ���������� ������ ��� ����� mesh � �������� q
	double solution_rz_in_point(double r, double z); // ������ ������� ��������� ������ � �����
	void grad_rz_in_point(double x, double y, double z, std::vector<double>& grad); // �������� ��������� �� ���������� ����������� ��� �������������� ������
	void grad_rz_in_point2(double x, double y, double z, std::vector<double>& grad); // �������� ��������� �� ���������� ����������� ��� �������������� ������
	void solution_export(class MESH& mesh, std::string solution_filename, std::vector<double>& solution);
	void solution_rz_export(std::string filename);

	void solution_rz_to_xyz2(std::string filename); // ��������� ����� ������� ��������� ������ � ����������
private:
	void add_G_b_p(); // ������� ���������� ������� ��������� � ������ �����
	void add_G_b_rz_p();
	void edge_cond_1_p();
	void edge_cond_1_rz_p();
	void print_solution_p_miss();
	void print_solution_rz_p_miss(); // ����� ������� ������������� ����������� �������

	// u - ������� ������� ����������������� ���������
	// f - ������� ������ ����� ����������������� ���������
	double u(double x, double y, double z);
	double f(double x, double y, double z);
	double u_rz(double r, double z);
	double f_rz(double r, double z);

	// C������ ��������� ������ �����
	void local_b_p(vector<double>& b, double h_x, double h_y, double h_z, int i);
	void local_b_rz_p(vector<double>& b, double h_r, double h_z, int i); // ���� �� �������
	void local_b_rz_pn(vector<double>& b, double h_r, double h_z, int i); // �������
	//---------------------------------------------------
};

#endif // FEM_RZ_H_
