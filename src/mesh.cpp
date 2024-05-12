#include "fem.h"

#define MESH_R_K 1.1
#define MESH_R_H0 1e-6
#define MESH_H_K 1.1
#define MESH_H_H0 1e-6

void gen_1d_mesh(vector<double>& sreda_x, vector<double>& hx, vector<double>& kx, vector<int>& left_right, vector<double>& x) {
	list<double> buff;
	list<double> buff_inv;

	double next;
	double k, h0;
	for (int i = 0; i < sreda_x.size() - 1; i++) {
		buff.resize(0);
		//Если требуется разрядка от правой границы
		if (left_right[i]) {
			//Значение первого шага
			h0 = hx[i];
			//Вычислим положение следущего узла
			next = sreda_x[i + 1] - h0;
			//Коэффициент разрядки
			k = kx[i];
			//Пока мы не перешли за левую границу будем добавлять узлы
			while (next > sreda_x[i]) {
				//Добавляем узел
				buff.push_back(next);
				//Вычисляем следующий
				h0 *= k;
				next -= h0;
			}
			//В конец кладем левую границу
			buff.push_back(sreda_x[i]);

			//Обходим buff в обратном порядке и записываем значения в x
			x.insert(x.end(), buff.rbegin(), buff.rend());
		}
		//Если разрядка от левой границы
		else {
			//В начало кладем левую границу
			buff.push_back(sreda_x[i]);
			//Значение первого шага
			h0 = hx[i];
			//Вычислим положение следущего узла
			next = sreda_x[i] + h0;
			//Коэффициент разрядки
			k = kx[i];
			//Пока мы не перешли за правую границу, будем добавлять узлы
			while (next < sreda_x[i + 1]) {
				//Добавляем узел
				buff.push_back(next);
				//Вычисляем следующий
				h0 *= k;
				next += h0;
			}

			//Обходим buff и записываем значения в x
			x.insert(x.end(), buff.begin(), buff.end());
		}
	}
	x.push_back(sreda_x[sreda_x.size() - 1]);

}

// Генерация одномерной сетки r для двумерной задачи
void gen_1d_mesh_r(class SREDA& sreda, vector<double>& r) {
	list<double> buff;
	list<double> buff_inv;

	double next;
	double k = MESH_R_K, h0 = MESH_R_H0,hx,hy;
	double r0, r1;
	// Для выбора радуиса на котором будет решаться двумерная задача выберается максимально возможное расстояние до границы от
	// источника: В случае квадрата его диагональ равна a*sqrt(2)
	hx = sreda.x[sreda.x.size() - 1] - sreda.x[0];
	hy = sreda.y[sreda.y.size() - 1] - sreda.y[0];
	
	r1 = max(hx, hy) * sqrt(2);
	buff.resize(0);
	next = h0;

	r.push_back(0); //Сперва занесем 0 
	//Пока мы не перешли за правую границу, будем добавлять узлы
	while (next < r1) {
		//Добавляем узел
		buff.push_back(next);
		//Вычисляем следующий
		h0 *= k;
		next += h0;
	}

	//Обходим buff и записываем значения в x
	r.insert(r.end(), buff.begin(), buff.end());

	// Заносим последний узел
	r.push_back(r1);
}

// Генерация одномерной сетки h для двумерной задачи
void gen_1d_mesh_h(class SREDA& sreda, vector<double>& h) {
	list<double> buff;
	list<double> buff_inv;

	double next;
	double k = MESH_H_K, h0 = MESH_H_H0, hx, hy;
	double h_end;

	h_end = sreda.z[sreda.z.size()-1];
	buff.resize(0);
	next = h0;

	h.push_back(0); //Сперва занесем 0 
	//Пока мы не перешли за правую границу, будем добавлять узлы
	while (next < h_end) {
		//Добавляем узел
		buff.push_back(next);
		//Вычисляем следующий
		h0 *= k;
		next += h0;
	}

	//Обходим buff и записываем значения в x
	h.insert(h.end(), buff.begin(), buff.end());

	// Заносим последний узел
	h.push_back(h_end);
}


//Подробит сетку
void splitting(vector<double>& x, int k) {
	if (k) {
		list<double> buff;
		double h;
		for (int i = 0; i < x.size() - 1; i++) {
			buff.push_back(x[i]);

			h = (x[i + 1] - x[i]) / pow(2, k);
			for (int j = 1; j < pow(2, k); j++) {
				buff.push_back(x[i] + (j * h));
			}
		}
		buff.push_back(x[x.size() - 1]);

		x.resize(0);
		x.insert(x.begin(), buff.begin(), buff.end());
	}
}

// Функция для объединения двух отсортированных векторов X[] и Y[]
std::vector<double> vector_merge(std::vector<double> const& X, std::vector < double > const& Y)
{
	int k = 0, i = 0, j = 0;
	std::vector<double> aux(X.size() + Y.size());

	// пока есть элементы в первом и втором массивах
	while (i < X.size() && j < Y.size())
	{
		if (X[i] <= Y[j]) {
			aux[k++] = X[i++];
		}
		else {
			aux[k++] = Y[j++];
		}
	}

	// копируем оставшиеся элементы
	while (i < X.size()) {
		aux[k++] = X[i++];
	}

	while (j < Y.size()) {
		aux[k++] = Y[j++];
	}

	aux.resize(std::distance(aux.begin(), std::unique(aux.begin(), aux.end())));
	
	return aux;
}

// Добавит узлы аномальных элементов в сетку по x
void add_elements_x_edges(class SREDA& sreda, std::vector<double>& x) {
	std::set<double> buff;
	for (int i = 0; i < sreda.elms.size(); i++) {
		if (sreda.elms[i][ANOMAL]) { // Если указано что элемент относится к аномальному полю то включаем его границы
			if (sreda.x[0] <= sreda.elms[i][X0_COORD] && sreda.elms[i][X0_COORD] <= sreda.x[sreda.x.size() - 1]) {
				buff.insert(sreda.elms[i][X0_COORD]);
			}
			if (sreda.x[0] <= sreda.elms[i][X1_COORD] && sreda.elms[i][X1_COORD] <= sreda.x[sreda.x.size() - 1]) {
				buff.insert(sreda.elms[i][X1_COORD]);
			}
		}
	}
	std::vector<double> nodes;
	nodes.insert(nodes.end(), buff.begin(), buff.end());

	// Добавим в сетку новые узлы и удалим повторяющиеся
	x = vector_merge(x, nodes);
	x.resize(std::distance(x.begin(), std::unique(x.begin(), x.end())));
}

// Добавит узлы аномальных элементов в сетку по y
void add_elements_y_edges(class SREDA& sreda, std::vector<double>& y) {
	std::set<double> buff;
	for (int i = 0; i < sreda.elms.size(); i++) {
		if (sreda.elms[i][ANOMAL]) { // Если указано что элемент относится к аномальному полю то включаем его границы
			if (sreda.y[0] <= sreda.elms[i][Y0_COORD] && sreda.elms[i][Y0_COORD] <= sreda.y[sreda.y.size() - 1]) {
				buff.insert(sreda.elms[i][Y0_COORD]);
			}
			if (sreda.y[0] <= sreda.elms[i][Y1_COORD] && sreda.elms[i][Y1_COORD] <= sreda.y[sreda.y.size() - 1]) {
				buff.insert(sreda.elms[i][Y1_COORD]);
			}
		}
	}
	std::vector<double> nodes;
	nodes.insert(nodes.end(), buff.begin(), buff.end());

	// Добавим в сетку новые узлы и удалим повторяющиеся
	y = vector_merge(y, nodes);
	y.resize(std::distance(y.begin(), std::unique(y.begin(), y.end())));
}

// Добавит узлы аномальных элементов в сетку по z
void add_elements_z_edges(class SREDA& sreda, std::vector<double>& z) {
	std::set<double> buff;
	for (int i = 0; i < sreda.elms.size(); i++) {
		if (sreda.elms[i][ANOMAL]) { // Если указано что элемент относится к аномальному полю то включаем его границы
			if (sreda.z[0] <= sreda.elms[i][Z0_COORD] && sreda.elms[i][Z0_COORD] <= sreda.z[sreda.z.size() - 1]) {
				buff.insert(sreda.elms[i][Z0_COORD]);
			}
			if (sreda.z[0] <= sreda.elms[i][Z1_COORD] && sreda.elms[i][Z1_COORD] <= sreda.z[sreda.z.size() - 1]) {
				buff.insert(sreda.elms[i][Z1_COORD]);
			}
		}
	}
	std::vector<double> nodes;
	nodes.insert(nodes.end(), buff.begin(), buff.end());

	// Добавим в сетку новые узлы и удалим повторяющиеся
	z = vector_merge(z, nodes);
}


void  MESH::gen_mesh(class SREDA& sreda) {
	//Генерация сетки по x
	gen_1d_mesh(sreda.x, sreda.hx, sreda.kx, sreda.left_right[0], x);
	splitting(x, sreda.splitting[0]);
	add_elements_x_edges(sreda, x);

	//Генерация сетки по y
	gen_1d_mesh(sreda.y, sreda.hy, sreda.ky, sreda.left_right[1], y);
	splitting(y, sreda.splitting[1]);
	add_elements_y_edges(sreda, y);

	//Генерация сетки по z	
	gen_1d_mesh(sreda.z, sreda.hz, sreda.kz, sreda.left_right[2], z);
	splitting(z, sreda.splitting[2]);
	add_elements_z_edges(sreda, z);

	kuzlov = x.size() * y.size() * z.size();
	kel = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);

	printf("3D: kuzlov - %d ; kel - %d; x - %d ; y - %d, z - %d\n",
		kuzlov, kel,x.size(),y.size(),z.size());
}

// Генерация двумерной сетки
void  MESH::gen_mesh_rh(class SREDA& sreda) {
	//Генерация сетки по r
	gen_1d_mesh_r(sreda, r);

	//Генерация сетки по h
	gen_1d_mesh_h(sreda, h);

	kuzlov = r.size() * h.size();
	kel = (r.size() - 1) * (h.size() - 1);

	printf("2D: kuzlov - %d ; kel - %d; r - %d, z - %d\n", kuzlov, kel, r.size(), h.size());
}

// Генерирует номера узлов элементов для трехмерной сетки
void gen_nvtr(class MESH* mesh,
	vector<vector<int>>& nvtr)
{
	double kel = mesh->kel;
	nvtr.resize(mesh->kel);
	for (int i = 0; i < kel; i++) { nvtr[i].resize(8); }

	int x, y, z;
	// Количество элементов в слое x
	int kx = (mesh->x.size() - 1);
	// Количество элементов в слое xy
	int kxy = (mesh->x.size() - 1) * (mesh->y.size() - 1);
	for (int i = 0; i < kel; i++) {
		x = i % kx + 1;
		y = (i % kxy) / (kx)+1;
		z = i / kxy + 1;

		nvtr[i][0] = x + (y - 1) * (kx + 1) + (z - 1) * mesh->x.size() * mesh->y.size();
		nvtr[i][1] = nvtr[i][0] + 1;
		nvtr[i][2] = nvtr[i][0] + mesh->x.size();
		nvtr[i][3] = nvtr[i][2] + 1;
		nvtr[i][4] = nvtr[i][0] + mesh->x.size() * mesh->y.size();
		nvtr[i][5] = nvtr[i][4] + 1;
		nvtr[i][6] = nvtr[i][4] + mesh->x.size();
		nvtr[i][7] = nvtr[i][6] + 1;
	}
}

// Генерирует номера узлов элементов для двумерной сетки
void gen_nvtr_rh(class MESH* mesh,
	vector<vector<int>>& nvtr)
{
	double kel = mesh->kel;
	nvtr.resize(mesh->kel);
	for (int i = 0; i < kel; i++) { nvtr[i].resize(4); }

	int n;
	for (int i = 0; i < kel; i++) {
		n = i / (mesh->r.size()-1)*(mesh->r.size()) + i%(mesh->r.size()-1) + 1;

		nvtr[i][0] = n;
		nvtr[i][1] = n + 1;
		nvtr[i][2] = n + mesh->r.size();
		nvtr[i][3] = n + mesh->r.size() + 1;

	}
}


void gen_coord(class MESH* mesh,
	vector<vector<double>>& rz) {
	rz.resize(mesh->kuzlov);
	int xs = mesh->x.size();
	int ys = mesh->y.size();
	int zs = mesh->z.size();
	int x, y, z;
	for (int i = 0; i < mesh->kuzlov; i++) {
		x = i % xs;
		y = (i % (xs * ys)) / xs;
		z = i / (xs * ys);
		rz[i].resize(3);
		rz[i] = { mesh->x[x], mesh->y[y], mesh->z[z] };
	}

}

void gen_coord_rh(class MESH* mesh,
	vector<vector<double>>& coord) {
	coord.resize(mesh->kuzlov);
	int rs = mesh->r.size();
	int hs = mesh->h.size();
	int r, h;
	for (int i = 0; i < mesh->kuzlov; i++) {
		r = i % rs;
		h = i / rs;
		coord[i].resize(2);
		coord[i] = { mesh->r[r], mesh->h[h]};
	}

}

// По кооринатам определит номера материала
// Обход обратный: будет получен материал аномального поля, если он не существует, то нормального
int get_mat_from_sreda_reverse(class SREDA& sreda, double x, double y, double z) {
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
int get_mat_from_sreda_direct(class SREDA& sreda, double x, double y, double z) {
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

// По кооринатам определит номера материала для двумерной задачи
int get_mat_from_sreda_rz(class SREDA& sreda, double r, double h) {
	// Будем обходить подобласти в прямом порядке (учитывая что в файле среда подоблости заданны вложением по порядку)
	// Считается, что для двумерной задачи элементы, которые не располагаются от начала до конца границы являются аномальными
	for (int i = sreda.elms.size() - 1; i >= 0; i--) {
		if (sreda.elms[i][Z0_COORD] < h && h < sreda.elms[i][Z1_COORD]) {
			return sreda.elms[i][MAT_N];
		}
	}
	return 0;
}

void gen_nvkat2d(class SREDA& sreda, class MESH* mesh, vector<int>& nvkat2d, vector<vector<int>>& nvtr,
	vector<vector<double>>& rz) {
	nvkat2d.resize(mesh->kel);
	double x, y, z;
	for (int i = 0; i < mesh->kel; i++) {
		// Координаты центров элементов
		x = (rz[nvtr[i][0] - 1][0] + rz[nvtr[i][1] - 1][0]) / 2.;
		y = (rz[nvtr[i][0] - 1][1] + rz[nvtr[i][2] - 1][1]) / 2.;
		z = (rz[nvtr[i][0] - 1][2] + rz[nvtr[i][4] - 1][2]) / 2.;
		nvkat2d[i] = get_mat_from_sreda_reverse(sreda, x, y, z);
	}
}

void gen_nvkat2d_rh(class SREDA& sreda, class MESH* mesh, vector<int>& nvkat2d, vector<vector<int>>& nvtr,
	vector<vector<double>>& rz) {
	nvkat2d.resize(mesh->kel);
	double r, z;
	for (int i = 0; i < mesh->kel; i++) {
		// Координаты центров элементов
		r = (rz[nvtr[i][0] - 1][0] + rz[nvtr[i][1] - 1][0]) / 2.;
		z = (rz[nvtr[i][0] - 1][1] + rz[nvtr[i][2] - 1][1]) / 2.;
		nvkat2d[i] = get_mat_from_sreda_rz(sreda, r, z);
	}
}

void gen_l1(class MESH* mesh, vector<double>& edge_conditions, vector<int>& l1, vector<vector<double>>& rz) {
	l1.resize(0);
	set<int> buff;
	for (int i = 0; i < mesh->kuzlov; i++) {
		// Если указано, что на грани x = x0 заданы первые краевые,
		// то если координата текущей точки совпадает с x0, тогда добавить ее номер.
		// С остальными гранями аналогично
		if (edge_conditions[L1_X0]) { if (rz[i][0] == mesh->x[0]) { buff.insert(i + 1); } }
		if (edge_conditions[L1_X1]) { if (rz[i][0] == mesh->x[mesh->x.size() - 1]) { buff.insert(i + 1); } }

		if (edge_conditions[L1_Y0]) { if (rz[i][1] == mesh->y[0]) { buff.insert(i + 1); } }
		if (edge_conditions[L1_Y1]) { if (rz[i][1] == mesh->y[mesh->y.size() - 1]) { buff.insert(i + 1); } }

		if (edge_conditions[L1_Z0]) { if (rz[i][2] == mesh->z[0]) { buff.insert(i + 1); } }
		if (edge_conditions[L1_Z1]) { if (rz[i][2] == mesh->z[mesh->z.size() - 1]) { buff.insert(i + 1); } }
	}
	l1.insert(l1.begin(), buff.begin(), buff.end());
}

void gen_l1_rh(class MESH* mesh, vector<double>& edge_conditions, vector<int>& l1, vector<vector<double>>& coord) {
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
	for (int i = 0; i < mesh->kuzlov; i++) {
		if (coord[i][0] == mesh->r[mesh->r.size()-1]) { buff.insert(i + 1); }
		else if (coord[i][1] == mesh->h[mesh->h.size()-1]) { buff.insert(i + 1); }
	}
	l1.insert(l1.begin(), buff.begin(), buff.end());
}


// nvtr : 4 * int номера вершин прямоугольников
// nvkat2d : 1 * int номера матреиала прямоугольников
// coord : 3 * double (r, z) координаты вершин
// l1 : 1 * int (номера вершин с первым нулевым краевым условием)
// sigma : 1 * double значение параметра сигма по номеру материала
//void MESH::gen_structures(class FEM& fem, class SREDA& sreda)
//{
//	gen_nvtr(this, fem.nvtr);
//	gen_coord(this, fem.coord);
//	gen_nvkat2d(sreda, this, fem.nvkat2d, fem.nvtr, fem.coord);
//	gen_l1(this, sreda.edge_conditions, fem.l1, fem.coord);
//}

// nvtr : 8 * int номера вершин прямоугольников
// nvkat2d : int nvkat2d номера матреиала прямоугольников
// rz : 3 * double (x, y, z) координаты вершин
// l1 : 1 * int (номера вершин с первым нулевым краевым условием)
// sigma : 1 * double значение параметра сигма по номеру материала
//void MESH::gen_structures_rh(class FEM& fem, class SREDA& sreda)
//{
//	gen_nvtr_rh(this, fem.nvtr);
//	gen_coord_rh(this, fem.coord);
//	gen_nvkat2d_rh(sreda, this, fem.nvkat2d, fem.nvtr, fem.coord);
//	fem.nvtr_rh = fem.nvtr; // Сохраним структуру для дальнейшего пользования
//	gen_l1_rh(this, sreda.edge_conditions, fem.l1, fem.coord);
//}

// Генерация двумерной сетки
void  MESH_RZ::gen_mesh(class SREDA_RZ & sreda_rz, class SREDA & sreda) {
	////Генерация сетки по r
	//gen_1d_mesh_r(sreda, r);

	////Генерация сетки по h
	//gen_1d_mesh_h(sreda, z);

	//kuzlov = r.size() * z.size();
	//kel = (r.size() - 1) * (z.size() - 1);

	vector<double> hx = { sreda_rz.r[H0] };
	vector<double> kx = { sreda_rz.r[K] };
	vector<int> left_right = { 0 };
	vector < double> x = { 0,sreda_rz.r[END] };

		//Генерация сетки по r
	gen_1d_mesh(x, hx, kx, left_right, r);
	splitting(r, sreda_rz.splitting[0]);


	//Генерация сетки по h
	hx[0] = sreda_rz.z[H0];
	kx[0] = sreda_rz.z[K];
	x[1] = sreda_rz.z[END];
	gen_1d_mesh(x, hx, kx, left_right, z);
	splitting(z, sreda_rz.splitting[1]);
	add_elements_z_edges(sreda, z);


	kuzlov = r.size() * z.size();
	kel = (r.size() - 1) * (z.size() - 1);

	printf("2D: kuzlov - %d ; kel - %d; r - %d, z - %d\n", kuzlov, kel, r.size(), z.size());
}

void MESH::add_mesh(class MESH& mesh) {
	x = vector_merge(x, mesh.x);
	y = vector_merge(y, mesh.y);
	z = vector_merge(z, mesh.z);

	kuzlov = x.size() * y.size() * z.size();
	kel = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);

	printf("3D SUM NORMAL: kuzlov - %d ; kel - %d; x - %d ; y - %d, z - %d\n",
		kuzlov, kel, x.size(), y.size(), z.size());
	
}