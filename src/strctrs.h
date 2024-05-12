#ifndef STRCTRS_H_
#define STRCTRS_H_

#include "mesh.h"
#include "sreda.h"

class STRCTRS
{
public:
    vector<vector<int>> nvtr; // mesh.kel записей номеров узлов элементов
    vector<int> nvkat2d;  // mesh.kel записей номеров материала элементов
    vector<int> nvkat2dr;  // mesh.kel записей номеров материалов элементов для нормального поля
    vector<vector<double>> coord; // mesh.kuzlov координат точек
    vector<int> l1; // номера узлов с первыми краевыми условиями

public: // Вызов функций генерации структур
    void gen_structures(SREDA& sreda, MESH& mesh);
    void gen_structures_rz(SREDA& sreda, MESH_RZ& mesh);

    void gen_structures_p(SREDA& sreda, MESH& mesh);
    void gen_structures_rz_p(SREDA& sreda, MESH_RZ& mesh); // Тест на полиномах


private: // Функции генерации структур для трехмерной задачи
    void gen_nvtr(class MESH& mesh);
    void gen_coord(class MESH& mesh);
    void gen_nvkat2d(class SREDA& sreda, class MESH& mesh);
    void gen_nvkat2dr(class SREDA& sreda, class MESH& mesh); // Номера материалов для нормальной задачи
    void gen_l1(class MESH& mesh, vector<double>& edge_conditions); 
    void gen_l1_p(class MESH& mesh, vector<double>& edge_conditions);

private: // Функции генерации структур для двумерной задачи
    void gen_nvtr_rz(class MESH_RZ& mesh);
    void gen_coord_rz(class MESH_RZ& mesh);
    void gen_nvkat2d_rz(class MESH_RZ& mesh, class SREDA& sreda);
    void gen_l1_rz(class MESH_RZ& mesh, vector<double>& edge_conditions);
    void gen_l1_rz_p(class MESH_RZ& mesh, vector<double>& edge_conditions); // Тест на полиномах


private:
    int get_mat_from_sreda_reverse(class SREDA& sreda, double x, double y, double z);
    int get_mat_from_sreda_direct(class SREDA& sreda, double x, double y, double z);
    int get_mat_from_sreda_rz(class SREDA& sreda, double z);
};

#endif // STRCTRS_H_