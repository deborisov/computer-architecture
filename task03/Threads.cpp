//Borisov Dmitry csbse192 
//Variant 4
#include <iostream>
#include <thread>
#include <string>
#include <fstream>
#include <iomanip>

//-----------------Reads matrix and its dimension from file-------------- 
void read_matrix(double**& matrix, std::string path_to_file, size_t& n) {
	std::fstream in(path_to_file, std::ios::in);
	if (in.is_open()) {
		in >> n;
		matrix = new double* [n];
		for (size_t i = 0; i < n; ++i) {
			matrix[i] = new double[n];
			for (size_t j = 0; j < n; ++j) {
				in >> matrix[i][j];
			}
		}
	}
	in.close();
}

//--------------------Finds minor of matrix with exlusion of exclusive_row row and exclusive_col colomn.---------
//--------------------Writes it to tmp matrix. cur_dim - dimension of matrix-------------------------------------
void cofactor(double** const& matrix, double**& tmp, size_t exclusive_row, size_t exclusive_col, size_t cur_dim) {
	size_t i = 0, j = 0;
	for (size_t row = 0; row < cur_dim; ++row) {
		for (size_t col = 0; col < cur_dim; ++col) {
			if (row != exclusive_row && col != exclusive_col) {
				tmp[i][j] = matrix[row][col];
				++j;
				if (j == cur_dim - 1) {
					j = 0;
					++i;
				}
			}
		}
	}
}

//--------------------------Find determinant of matrix with dimension n--------------------
double determinant(double** const& matrix, size_t n) {
	double det = 0;
	int sign = 1;
	if (n == 1) {
		return matrix[0][0];
	}
	double** tmp = new double* [n - 1];
	for (size_t i = 0; i < n - 1; ++i) {
		tmp[i] = new double[n - 1];
	}
	for (size_t i = 0; i < n; ++i) {
		cofactor(matrix, tmp, 0, i, n);
		det += sign * determinant(tmp, n - 1) * matrix[0][i];
		sign = -sign;
	}
	for (size_t i = 0; i < n - 1; ++i) {
		delete[] tmp[i];
	}
	delete[] tmp;
	return det;
}

//-----------------------Finds and writes minors to inv_matrix for elements-----------------------------
//-----------------------Elements numbers depend on number_of_threads and thread number-----------------
void thread_function(double** const& matrix, double**& inv_matrix, size_t number_of_threads, size_t n, size_t thread_number) {
	//std::cout << std::this_thread::get_id() << std::endl;
	double** tmp = new double* [n - 1];
	for (size_t i = thread_number; i < n * n; i += number_of_threads) {
		for (size_t j = 0; j < n - 1; ++j) {
			tmp[j] = new double[n - 1];
		}
		size_t row = i / n, col = i % n;
		cofactor(matrix, tmp, row, col, n);
		inv_matrix[col][row] = (row + col) % 2 == 0 ? determinant(tmp, n - 1) : - determinant(tmp, n -1);
	}
	for (size_t i = 0; i < n - 1; ++i) {
		delete[] tmp[i];
	}
	delete[] tmp;
}

//-------------------Multithread function writes inverted matrix to inv_matrix------------------
//-------------------If matrix is non-invertable returns false----------------------------------
bool invert_matrix(double** const& matrix, double**& inv_matrix, size_t number_of_threads, size_t n) {
	std::thread* threads = new std::thread[number_of_threads];
	for (size_t i = 0; i < number_of_threads; ++i) {
		threads[i] = std::thread(thread_function, std::ref(matrix), std::ref(inv_matrix), number_of_threads, n, i);
	}
	double det = determinant(matrix, n);
	for (size_t i = 0; i < number_of_threads; ++i) {
		threads[i].join();  
	}
	if (det == 0) {
		delete[] threads;
		return false;
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			inv_matrix[i][j] /= det;
		}
	}
	delete[] threads;
	return true;
}

//--------------------Function writes inverse matrix to a file or ------------
//--------------------a message that it was non-invertable--------------------
void write_output(double** const& inv_matrix, size_t n, std::string path_to_write, bool invertable) {
	std::fstream out(path_to_write, std::ios::out);
	if (out.is_open()) {
		if (invertable) {
			out << std::fixed << std::setprecision(2);
			for (size_t i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					out << inv_matrix[i][j] << " ";
				}
				out << std::endl;
			}
		}
		else {
			out << "Determinant of matrix is equal to zero. Matrix is non-invertable.";
		}
	}
	out.close();
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cout << "Please provide three argument: path to file with matrix, path to output file and number of threads.";
		return 1;
	}
	std::string path = argv[1];
	std::string path_to_write = argv[2];
	if ((int)(argv[3]) < 1) {
		std::cout << "Please input positive number of threads.";
		return 1;
	}
	size_t number_of_threads = std::stoi(argv[3]);
	size_t n;
	double** matrix = nullptr;
	read_matrix(matrix, path, n);
	double** inv_matrix = new double* [n];
	for (size_t i = 0; i < n; ++i) {
		inv_matrix[i] = new double[n];
	}
	bool invertable = invert_matrix(matrix, inv_matrix, number_of_threads, n);
	write_output(inv_matrix, n, path_to_write, invertable);
	for (size_t i = 0; i < n; ++i) {
		delete[] matrix[i];
		delete[] inv_matrix[i];
	}
	delete[] matrix;
	delete[] inv_matrix;
}

