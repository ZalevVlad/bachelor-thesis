#ifndef SREDA_H_
#define SREDA_H_

#include <vector>
#include <iostream>
#include <fstream>

#define SREDA_FILENAME "/sreda"
#define SREDA_N_FILENAME "/sreda_n"
#define EDGE_CONDITIONS_FILENAME "/edge_conditions"
#define CURRENT_SOURCES_FILENAME "/current_sources"
#define MATERIALS_FILENAME "/materials"
#define POINTS_FILENAME "/points"
#define POINTS_RZ_FILENAME "/points_rz"


using namespace std;
// ������, ��������� ���, ������������ ���, ����������� ��������
enum {END, H0, HMAX, K};

enum sigma_strucure { SIGMA, SIG_ANOMAL };

class SREDA_RZ {
public:
    // ������, ��������� ���, ������������ ���, ����������� ��������
    std::vector<double> r;
    std::vector<double> z;
    // ��������� ��������� ����� �� r � �� z
    std::vector<double> splitting; 

    //���� points_rz
    std::vector<std::vector<double>> points;

    void read_points(string filename);

    // ������� ���������� � ��������� ������
    void read_sreda(string filename);
};

class SREDA
{
public:
    std::vector < std::vector<double>> elms;

    std::vector<double> x;
    std::vector<double> hx;
    std::vector<double> kx;

    std::vector<double> y;
    std::vector<double> hy;
    std::vector<double> ky;

    std::vector <double> z;
    std::vector<double> hz;
    std::vector<double> kz;
    std::vector<std::vector<int>> left_right;

    std::vector<int> splitting;

    //���� edge_conditions
    vector<double> edge_conditions;

    //���� current_sources
    vector<vector<double>> current_sources;

    //���� materials
    vector<vector<double>> sigma;

    // ���� points
    vector<vector<double>> points;

    //��������� ���� sreda
    void read_sreda(const char* filename);
    void read_sreda(string filename);

    //��������� ���� edge_conditions
    void read_edge_conditions(const char* filename);
    void read_edge_conditions(string filename);

    //��������� ���� current_sources
    void read_current_sources(const char* filename);
    void read_current_sources(string filename);

    //��������� ���� materials
    void read_materials(const char* filename, class FEM_old& fem);
    void read_materials(string filename);

    // ��������� ���� points
    void read_points(string filename);

    // ������� ��� ������� ����� ��� ���������� ������
    void read_problem(string meshname);
    void read_problem_n(string meshname); // ���� sreda_n
};
#endif // SREDA_H_
