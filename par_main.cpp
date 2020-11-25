#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <fstream>

#include <omp.h>

double dotKernel(const std::vector<double>& vec_1, const std::vector<double>& vec_2) {
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (size_t i = 0; i < vec_1.size(); ++i) {
		sum += vec_1[i] * vec_2[i];
	}
	return sum;
}

void axpbyKernel(const double a, std::vector<double>& x,
	             const double b, const std::vector<double>& y) {
#pragma omp parallel for
	for (size_t i = 0; i < x.size(); ++i) {
		x[i] = a * x[i] + b * y[i];
	}
	return;
}

void SpMVKernel(const std::vector<double>& matrix, const std::vector<int>& row_ptr,
	            const std::vector<int>& col_ptr, const std::vector<double>& vec,
	            std::vector<double>& res) {
#pragma omp parallel for
	for (size_t i = 0; i < res.size(); ++i) {
		double sum = 0;
		for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			sum += matrix[j] * vec[col_ptr[j]];
		}
		res[i] = sum;
	}
	return;
}

int obliquesBefore(const int K1, const int K2, const int cell_idx) {
	return cell_idx / (K1 + K2) * K2 + std::max(cell_idx % (K1 + K2) - K1, 0);
}

bool hasOblique(const int K1, const int K2, const int cell_idx) {
	return cell_idx / (K1 + K2) * (K1 + K2) + K1 <= cell_idx;
}

void generate(const int Nx, const int Ny, const int K1, const int K2,
	          std::vector<int>& row_ptr, std::vector<int>& col_ptr) {
	row_ptr[0] = 0;
	int steps = 100;
#pragma omp parallel for schedule(dynamic)
	for (int idx = 0; idx < (Nx + 1) * (Ny + 1); idx += steps) {
		const int cur_thread_steps = std::min((Nx + 1) * (Ny + 1) - idx, steps);
		int i = idx / (Nx + 1);
		int j = idx % (Nx + 1);
		int col_idx = 0;
		if (i > 1) {
				col_idx += 2 * (Nx + 1) * (i - 1);	
			}
		if (i > 0) {
			col_idx += 2 * j + (Nx + 1 - j);
		}
		if (i != Ny) {
			col_idx += j;	
		}
		col_idx += 2 * Nx * i;
		if (j > 0) {
			col_idx += 1;
			col_idx += 2 * (j - 1);
		}
		if (j != 0 and i != Ny) {
			int tmp = Nx * i + j;
			col_idx += obliquesBefore(K1, K2, tmp);
			if (i != 0) {
				tmp = Nx * i;
				col_idx -= obliquesBefore(K1, K2, tmp);
			}
		}
		if (i != 0) {
			int tmp = Nx * i;
			tmp = obliquesBefore(K1, K2, tmp);
			col_idx += 2 * tmp;
			int ptrL, ptrB = Nx * i;
			if (j == 0) {
				ptrL = ptrB - Nx;
			} else {
				ptrL = Nx * (i - 1) + (j - 1);
			}
			col_idx -= obliquesBefore(K1, K2, ptrB) - obliquesBefore(K1, K2, ptrL);
		}
		col_idx += (Nx + 1) * i + j;
		for (int k = 0; k < cur_thread_steps; ++k) {
			i = (idx + k) / (Nx + 1);
			j = (idx + k) % (Nx + 1);
			if (i != 0 and j != 0 and hasOblique(K1, K2, Nx * (i - 1) + (j - 1))) {	
				col_ptr[col_idx] = (Nx + 1) * (i - 1) + (j - 1);
				++col_idx;
			}
			if (i != 0) {
				col_ptr[col_idx] = (Nx + 1) * (i - 1) + j;
				++col_idx;
			}
			if (j != 0) {
				col_ptr[col_idx] = (Nx + 1) * i + (j - 1);
				++col_idx;
			}
			col_ptr[col_idx] = (Nx + 1) * i + j;
			++col_idx;
			if (j != Nx) {
				col_ptr[col_idx] = (Nx + 1) * i + (j + 1);
				++col_idx;
			}
			if (i != Ny) {
				col_ptr[col_idx] = (Nx + 1) * (i + 1) + j;
				++col_idx;
			}
			if (i != Ny and j != Nx and hasOblique(K1, K2, Nx * i + j)) {
				col_ptr[col_idx] = (Nx + 1) * (i + 1) + (j + 1);
				++col_idx;
			}
			row_ptr[(Nx + 1) * i + j + 1] = col_idx;
		}
		
	}
	return;
}

void fill(const int Nx, const int Ny,
	      const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	      std::vector<double>& A_arr, std::vector<double>& b_vec) {
#pragma omp parallel for schedule(dynamic, 100)
	for (int cur_idx = 0; cur_idx < (Nx + 1) * (Ny + 1); ++cur_idx) {
		double sum = 0;
		int diag_idx;
		for (int k = row_ptr[cur_idx]; k < row_ptr[cur_idx + 1]; ++k) {
			const int neid_idx = col_ptr[k];
			if (neid_idx == cur_idx) {
				diag_idx = k;
			} else {
				A_arr[k] = std::cos(static_cast<double>(cur_idx) * neid_idx + cur_idx + neid_idx);
				sum += std::abs(A_arr[k]);
			}
		}
		A_arr[diag_idx] = sum * 1.234;
		b_vec[cur_idx] = std::sin(cur_idx);
	}
	return;
}

void Solve(const std::vector<double>& A_arr, const std::vector<double>& b_vec,
	       const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	       std::vector<double>& x_vec, const double TOL = 1e-7) {
	const size_t vec_size = b_vec.size();
	std::vector<double> p_vec(vec_size), z_vec(vec_size),
	                    q_vec(vec_size), r_vec(b_vec),
	                    inverse_M(vec_size);
	std::vector<int> M_row_ptr(vec_size + 1), M_col_ptr(vec_size);
#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(vec_size); ++i) {
		int diag_idx;
		for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
			if (col_ptr[k] == i) {
				diag_idx = k;
				break;
			}
		}
		inverse_M[i] = 1. / A_arr[diag_idx];
		M_row_ptr[i] = i;
		M_col_ptr[i] = i;
	}
	M_row_ptr[vec_size] = static_cast<int>(vec_size);
	bool do_cycle = true;
	int it_num = 1;
	const int MAX_IT_NUM = 100;
	double dot_time = 0., axpby_time = 0., spmv_time = 0.;
	dot_time -= omp_get_wtime();
	const double b_norm = dotKernel(b_vec, b_vec);
	dot_time += omp_get_wtime();
	const double EPS = TOL * b_norm;
	double alpha, rho_prev, rho_cur, beta, residual_norm;
	do {
		spmv_time -= omp_get_wtime();
		SpMVKernel(inverse_M, M_row_ptr, M_col_ptr, r_vec, z_vec);
		spmv_time += omp_get_wtime();
		dot_time -= omp_get_wtime();
		rho_cur = dotKernel(r_vec, z_vec);
		dot_time += omp_get_wtime();
		if (it_num == 1) {
#pragma omp parallel for
			for (size_t i = 0; i < z_vec.size(); ++i) {
				p_vec[i] = z_vec[i];
			}
		} else {
			beta = rho_cur / rho_prev;
			axpby_time -= omp_get_wtime();
			axpbyKernel(beta, p_vec, 1., z_vec);
			axpby_time += omp_get_wtime();
		}
		spmv_time -= omp_get_wtime();
		SpMVKernel(A_arr, row_ptr, col_ptr, p_vec, q_vec);
		spmv_time += omp_get_wtime();
		dot_time -= omp_get_wtime();
		alpha = rho_cur / dotKernel(p_vec, q_vec);
		dot_time += omp_get_wtime();
		axpby_time -= omp_get_wtime();
		axpbyKernel(1., x_vec,  alpha, p_vec);
		axpby_time += omp_get_wtime();
		axpby_time -= omp_get_wtime();
		axpbyKernel(1., r_vec, -alpha, q_vec);
		axpby_time += omp_get_wtime();
		dot_time -= omp_get_wtime();
		residual_norm = dotKernel(r_vec, r_vec);
		dot_time += omp_get_wtime();
		std::cout << "It - " << it_num << ", RES_NORM - " << residual_norm << ", tol - " << residual_norm / b_norm << '\n';
		if (residual_norm < EPS or it_num >= MAX_IT_NUM) {
			do_cycle = false;
		} else {	
			++it_num;
			rho_prev = rho_cur;
		}
	} while (do_cycle);
	std::cout << "Dot func time - " << dot_time <<
	             "\nAxpby func time - " << axpby_time <<
	             "\nSpmv func time - " << spmv_time << '\n';
	return;
}

void printHelp() {
	std::cout << "Incorrect input.\n"
	          << "Program mode:\n" 
	          << "- filename\n"
	          << "- filename debug_flag(+)\n"
	          << "- Nx Ny K1 K2 T\n"
	          << "- Nx Ny K1 K2 T debug_flag(+)\n"
	          << "File with name filename must contain parameters Nx Ny K1 K2 T debug_flag.\n"
	          << "Nx, Ny, K1, K2, T have integer type. 1 <= Nx * Ny <= 10^7.\n";
}

int main(int argc, char* argv[]) {
	int Nx, Ny, K1, K2, threads_num;
	if (argc == 2 or argc == 3) {
		std::ifstream in_file(argv[1]);
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
		in_file >> Nx >> Ny >> K1 >> K2 >> threads_num;
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
	} else if (argc == 6 or argc == 7) {
		try {
			Nx = std::stoi(argv[1]);
			Ny = std::stoi(argv[2]);
			K1 = std::stoi(argv[3]);
			K2 = std::stoi(argv[4]);
			threads_num = std::stoi(argv[5]);
		} catch (const std::invalid_argument&) {
			std::cout << "Invalid number" << '\n';
			printHelp();
			return 1;
		} catch (const std::out_of_range&) {
			std::cout << "Out of range" << '\n';
			printHelp();
			return 1;
		}
	} else {
		printHelp();
		return 1;
	}
	if (Nx <= 0 or Ny <= 0 or K1 <= 0 or K2 <= 0 or Nx * Ny >= 1e9) {
		std::cout << "Bad parameter values" << '\n';
		printHelp();
		return 1;
	}
	bool debug_flag = false;
	if (argc == 3 or argc == 7) {
		debug_flag = true;
	}

	omp_set_num_threads(threads_num);

	const size_t sz_row_ptr = (Nx + 1) * (Ny + 1) + 1;
	const size_t sz_col_ptr =
	    2 * (Nx * (Ny + 1) + Ny * (Nx + 1) +                                    // horizontal, vertical
	    (Nx * Ny) / (K1 + K2) * K2 + std::max((Nx * Ny) % (K1 + K2) - K1, 0)) + // oblique
        (Nx + 1) * (Ny + 1);                                                    // with itself
    // std::cout << "Graph vertex - " << sz_row_ptr - 1 << '\n';
    // std::cout << "Graph edges - " << (sz_col_ptr - (Nx + 1) * (Ny + 1)) / 2 << '\n';
    // std::cout << "Non-zero elements - " << sz_col_ptr << '\n';
    // size_t X = sz_row_ptr, Y = sz_col_ptr;
    // std::cout << (2 * X) << ' ' << 3 * X << ' ' << 2 * Y << '\n';
    // return 0;
        
	std::vector<int> row_ptr(sz_row_ptr);
	std::vector<int> col_ptr(sz_col_ptr);
	double time = -omp_get_wtime();
	generate(Nx, Ny, K1, K2, row_ptr, col_ptr);
	time += omp_get_wtime();
	std::cout << "Gen(parallel). Time - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("parallel_IA_JA.txt");
		for (int i = 0; i < static_cast<int>(row_ptr.size() - 1); ++i) {
			oFile << i << " - ";
			for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
				oFile << col_ptr[j] << ' ';
			}
			oFile << '\n';
		}
	}
	std::vector<double> A_arr(sz_col_ptr);
	std::vector<double> b_vec(sz_row_ptr - 1);
	time = -omp_get_wtime();
	fill(Nx, Ny, row_ptr, col_ptr, A_arr, b_vec);
	time += omp_get_wtime();
	std::cout << "Filling(parallel). Time - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("parallel_A_b.txt");
		for (int i = 0; i < static_cast<int>(row_ptr.size() - 1); ++i) {
			oFile << i << " - ";
			for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
				oFile << A_arr[j] << ' ';
			}
			oFile << " | " << b_vec[i];
			oFile << '\n';
		}
	}
	time = -omp_get_wtime();
	std::vector<double> x_vec(b_vec.size(), 0);
	Solve(A_arr, b_vec, row_ptr, col_ptr, x_vec);
	time += omp_get_wtime();
	std::cout << "Solve(parallel). Time - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("parallel_X.txt");
		for (int i = 0; i < static_cast<int>(x_vec.size()); ++i) {
			oFile << x_vec[i] << ' ';
		}
		oFile << '\n';
	}

	return 0;
}
