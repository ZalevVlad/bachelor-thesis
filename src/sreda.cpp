#include "fem.h"

//—читает вектор из потока
void read_vector(ifstream* ifs, vector<double>& x, int n) {
	x.resize(n);
	for (int i = 0; i < n; i++) { *ifs >> x[i]; }
}

void read_vector(ifstream* ifs, vector<int>& x, int n) {
	x.resize(n);
	for (int i = 0; i < n; i++) { *ifs >> x[i]; }
}

void SREDA::read_sreda(const char* filename) {
	cout << "Reading file sreda\n";

	ifstream ifs;
	ifs.open(filename, ios::in);

	int num_areas;
	ifs >> num_areas;
	elms.resize(num_areas);
	for (int i = 0; i < num_areas; i++) { read_vector(&ifs, elms[i], 8); }

	int x_num, y_num, z_num;
	left_right.resize(3);
	ifs >> x_num;
	read_vector(&ifs, x, x_num);
	read_vector(&ifs, hx, x_num - 1);
	read_vector(&ifs, kx, x_num - 1);
	read_vector(&ifs, left_right[0], x_num - 1);

	ifs >> y_num;
	read_vector(&ifs, y, y_num);
	read_vector(&ifs, hy, y_num - 1);
	read_vector(&ifs, ky, y_num - 1);
	read_vector(&ifs, left_right[1], y_num - 1);

	ifs >> z_num;
	read_vector(&ifs, z, z_num);
	read_vector(&ifs, hz, z_num - 1);
	read_vector(&ifs, kz, z_num - 1);
	read_vector(&ifs, left_right[2], z_num - 1);

	vector<int> split(3);
	read_vector(&ifs, split, 3);
	splitting = split;
}
void SREDA::read_sreda(string filename) {

	ifstream ifs;
	ifs.open(filename, ios::in);
	if (ifs.is_open()) {
		cout << "Reading file sreda: " << filename << endl;
		int num_areas;
		ifs >> num_areas;
		elms.resize(num_areas);
		for (int i = 0; i < num_areas; i++) { read_vector(&ifs, elms[i], 8); }

		int x_num, y_num, z_num;
		left_right.resize(3);
		ifs >> x_num;
		read_vector(&ifs, x, x_num);
		read_vector(&ifs, hx, x_num - 1);
		read_vector(&ifs, kx, x_num - 1);
		read_vector(&ifs, left_right[0], x_num - 1);

		ifs >> y_num;
		read_vector(&ifs, y, y_num);
		read_vector(&ifs, hy, y_num - 1);
		read_vector(&ifs, ky, y_num - 1);
		read_vector(&ifs, left_right[1], y_num - 1);

		ifs >> z_num;
		read_vector(&ifs, z, z_num);
		read_vector(&ifs, hz, z_num - 1);
		read_vector(&ifs, kz, z_num - 1);
		read_vector(&ifs, left_right[2], z_num - 1);

		vector<int> split(3);
		read_vector(&ifs, split, 3);
		splitting = split;
	}
	else {
		cout << "Error opening file sreda: " << filename << endl;
	}

}

void SREDA::read_edge_conditions(const char* filename) {
	ifstream ifs;
	ifs.open(filename, ios::in);

	edge_conditions.resize(6);
	for (int i = 0; i < 6; i++) { ifs >> edge_conditions[i]; }
}
void SREDA::read_edge_conditions(string filename) {
	ifstream ifs;
	ifs.open(filename, ios::in);

	edge_conditions.resize(6);
	for (int i = 0; i < 6; i++) { ifs >> edge_conditions[i]; }
}

void SREDA::read_current_sources(const char* filename) {
	ifstream ifs;
	ifs.exceptions(std::ios::failbit);   // throw if failbit get set
	try {
		ifs.open(filename, ios::in);
	}
	catch (const std::exception& ex) {
		std::cerr << "Could not open current sourcess: "
			<< ex.what();
	}
	

	int n;
	ifs >> n;

	current_sources.resize(n);
	for (int i = 0; i < n; i++) {
		current_sources[i].resize(5);
		ifs >> current_sources[i][SOURCE_X] >> current_sources[i][SOURCE_Y] >> current_sources[i][SOURCE_Z] >> current_sources[i][SOURCE_POW];
	}
}
void SREDA::read_current_sources(string filename) {
	ifstream ifs;
	ifs.exceptions(std::ios::failbit);   // throw if failbit get set
	try {
		ifs.open(filename, ios::in);
	}
	catch (const std::exception& ex) {
		std::cerr << "Could not open current sourcess: "
			<< ex.what();
	}


	int n;
	ifs >> n;

	current_sources.resize(n);
	for (int i = 0; i < n; i++) {
		current_sources[i].resize(5);
		ifs >> current_sources[i][SOURCE_X] >> current_sources[i][SOURCE_Y] >> current_sources[i][SOURCE_Z] >> current_sources[i][SOURCE_POW];
	}
}

//void SREDA::read_materials(const char* filename, class FEM_old& fem) {
//	ifstream ifs;
//	ifs.exceptions(std::ios::failbit);   // throw if failbit get set
//	try {
//		ifs.open(filename, ios::in);
//	}
//	catch (const std::exception& ex) {
//		std::cerr << "Could not open materials: "
//			<< ex.what();
//	}
//	int n;
//	ifs >> n;
//
//	fem.sigma.resize(n);
//	for (int i = 0; i < n; i++) {
//		fem.sigma[i].resize(2);
//		ifs >> fem.sigma[i][0];
//		ifs >> fem.sigma[i][1];
//	}
//}
void SREDA::read_materials(string filename) {
	ifstream ifs;
	ifs.exceptions(std::ios::failbit);   // throw if failbit get set
	try {
		ifs.open(filename, ios::in);
	}
	catch (const std::exception& ex) {
		std::cerr << "Could not open materials: "
			<< ex.what();
	}
	int n;
	ifs >> n;

	sigma.resize(n);
	for (int i = 0; i < n; i++) {
		sigma[i].resize(2);
		ifs >> sigma[i][0];
		ifs >> sigma[i][1];
	}
}

void SREDA::read_points(string filename) {
	ifstream ifs;
	ifs.open(filename, ios::in);
	std::vector<double> buf;
	while (!ifs.eof()) {
		read_vector(&ifs, buf, 3);
		points.push_back(buf);
	}
}


void SREDA::read_problem(string meshname) {
	read_sreda(meshname + SREDA_FILENAME);
	read_edge_conditions(meshname + EDGE_CONDITIONS_FILENAME);
	read_current_sources(meshname + CURRENT_SOURCES_FILENAME);
	read_materials(meshname + MATERIALS_FILENAME);
	read_points(meshname + POINTS_FILENAME);

}

void SREDA::read_problem_n(string meshname) {
	read_sreda(meshname + SREDA_N_FILENAME);
	read_edge_conditions(meshname + EDGE_CONDITIONS_FILENAME);
	read_current_sources(meshname + CURRENT_SOURCES_FILENAME);
	read_materials(meshname + MATERIALS_FILENAME);
	read_points(meshname + POINTS_FILENAME);
}

void SREDA_RZ::read_points(string filename) {
	ifstream ifs;
	ifs.open(filename, ios::in);
	std::vector<double> buf;
	while (!ifs.eof()) {
		read_vector(&ifs, buf, 2);
		points.push_back(buf);
	}
}

// —читает файл с описанием среды дл€ двумерной задачи
void SREDA_RZ::read_sreda(string filename) {
	cout << "Reading file sreda_rz\n";

	ifstream ifs;
	ifs.open(filename + "/sreda_rz", ios::in);

	read_vector(&ifs, r, 4);
	read_vector(&ifs, z, 4);
	read_vector(&ifs, splitting, 2);
	read_points(filename + POINTS_RZ_FILENAME);
}