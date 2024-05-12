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
// радиус, начальный шаг, максимальный шаг, коэффициент разрядки
enum {END, H0, HMAX, K};

enum sigma_strucure { SIGMA, SIG_ANOMAL };

class SREDA_RZ {
public:
    // радиус, начальный шаг, максимальный шаг, коэффициент разрядки
    std::vector<double> r;
    std::vector<double> z;
    // Параметры дробления сетки по r и по z
    std::vector<double> splitting; 

    //Файл points_rz
    std::vector<std::vector<double>> points;

    void read_points(string filename);

    // Прочтет информацию о двумерной задаче
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

    //Файл edge_conditions
    vector<double> edge_conditions;

    //Файл current_sources
    vector<vector<double>> current_sources;

    //Файл materials
    vector<vector<double>> sigma;

    // Файл points
    vector<vector<double>> points;

    //Прочитает файл sreda
    void read_sreda(const char* filename);
    void read_sreda(string filename);

    //Прочитает файл edge_conditions
    void read_edge_conditions(const char* filename);
    void read_edge_conditions(string filename);

    //Прочитает файл current_sources
    void read_current_sources(const char* filename);
    void read_current_sources(string filename);

    //Прочитает файл materials
    void read_materials(const char* filename, class FEM_old& fem);
    void read_materials(string filename);

    // Прочитает файл points
    void read_points(string filename);

    // Прочтет все входные файлы для трехмерной задачи
    void read_problem(string meshname);
    void read_problem_n(string meshname); // файл sreda_n
};
#endif // SREDA_H_
