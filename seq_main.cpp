#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <omp.h>

double dotKernel(const std::vector<double>& vec_1, const std::vector<double>& vec_2) {
	return std::inner_product(vec_1.cbegin(), vec_1.cend(), vec_2.cbegin(), 0.);
}

void axpbyKernel(const double a, std::vector<double>& x,
	             const double b, const std::vector<double>& y) {
	std::transform(x.cbegin(), x.cend(),
		           y.cbegin(), x.begin(),
	   			   [a, b](const double x_i, const double y_i) {
	   				   return a * x_i + b * y_i;
	   			   });
	return;
}

void SpMVKernel(const std::vector<double>& matrix, const std::vector<int>& row_ptr,
	            const std::vector<int>& col_ptr, const std::vector<double>& vec,
	            std::vector<double>& res) {
	for (int i = 0; i < static_cast<int>(res.size()); ++i) {
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
	size_t idx_col_ptr = 0;
	for (int i = 0; i < (Ny + 1); ++i) {
		for (int j = 0; j < (Nx + 1); ++j) {
			int neib_num = 0;
			if (i != 0 and j != 0 and hasOblique(K1, K2, Nx * (i - 1) + (j - 1))) {	
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * (i - 1) + (j - 1);
				++idx_col_ptr;
			}
			if (i != 0) {
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * (i - 1) + j;
				++idx_col_ptr;
			}
			if (j != 0) {
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * i + (j - 1);
				++idx_col_ptr;
			}
			++neib_num;
			col_ptr[idx_col_ptr] = (Nx + 1) * i + j;
			++idx_col_ptr;
			if (j != Nx) {
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * i + (j + 1);
				++idx_col_ptr;
			}
			if (i != Ny) {
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * (i + 1) + j;
				++idx_col_ptr;
			}
			if (i != Ny and j != Nx and hasOblique(K1, K2, Nx * i + j)) {
				++neib_num;
				col_ptr[idx_col_ptr] = (Nx + 1) * (i + 1) + (j + 1);
				++idx_col_ptr;
			}
			row_ptr[(Nx + 1) * i + j + 1] = row_ptr[(Nx + 1) * i + j] + neib_num;
		}
	}
	return;
}

void fill(const int Nx, const int Ny,
	      const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	      std::vector<double>& A_arr, std::vector<double>& b_vec) {
	for (int i = 0; i < (Ny + 1); ++i) {
		for (int j = 0; j < (Nx + 1); ++j) {
			const int cur_idx = (Nx + 1) * i + j;
			double sum = 0;
			int diag_idx;
			for (int k = row_ptr[cur_idx]; k < row_ptr[cur_idx + 1]; ++k) {
				const int neib_idx = col_ptr[k];
				if (neib_idx == cur_idx) {
					diag_idx = k;
				} else {
					A_arr[k] = std::cos(static_cast<double>(cur_idx) * neib_idx + cur_idx + neib_idx);
					sum += std::abs(A_arr[k]);
				}
			}
			A_arr[diag_idx] = sum * 1.234;
			b_vec[cur_idx] = std::sin(cur_idx);
		}
	}
	return;
}

void Solve(const std::vector<double>& A_arr, const std::vector<double>& b_vec,
	       const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	       std::vector<double>& x_vec, const double TOL = 1e-4) {
	const size_t vec_size = b_vec.size();
	std::vector<double> p_vec(vec_size), z_vec(vec_size), q_vec(vec_size), r_vec(b_vec);
	bool do_cycle = true;
	int it_num = 1;
	const int MAX_IT_NUM = 100;
	std::vector<double> inverse_M(vec_size);
	std::vector<int> M_row_ptr(vec_size + 1), M_col_ptr(vec_size);
//#pragma omp parallel for
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
	double dot_time = 0., axpby_time = 0., spmv_time = 0.;
	double alpha, ro_prev, ro_cur, beta;
	do {
		spmv_time -= omp_get_wtime();
		SpMVKernel(inverse_M, M_row_ptr, M_col_ptr, r_vec, z_vec);
		spmv_time += omp_get_wtime();
		
		dot_time -= omp_get_wtime();
		ro_cur = dotKernel(r_vec, z_vec);
		dot_time += omp_get_wtime();
		if (it_num == 1) {
			p_vec = z_vec; // мб распараллелить
		} else {
			beta = ro_cur / ro_prev;
			axpby_time -= omp_get_wtime();
			axpbyKernel(beta, p_vec, 1., z_vec);
			axpby_time += omp_get_wtime();
		}
		spmv_time -= omp_get_wtime();
		SpMVKernel(A_arr, row_ptr, col_ptr, p_vec, q_vec);
		spmv_time += omp_get_wtime();
		dot_time -= omp_get_wtime();
		alpha = ro_cur / dotKernel(p_vec, q_vec);
		dot_time += omp_get_wtime();
		axpby_time -= omp_get_wtime();
		axpbyKernel(1., x_vec,  alpha, p_vec);
		axpbyKernel(1., r_vec, -alpha, q_vec);
		axpby_time += omp_get_wtime();
		std::cout << it_num << " , " << ro_cur << '\n';
		if (ro_cur < TOL or it_num >= MAX_IT_NUM) {
			do_cycle = false;
		} else {	
			++it_num;
			ro_prev = ro_cur;
		}
	} while (do_cycle);
	std::cout << dot_time << ' ' << axpby_time << ' ' << spmv_time << '\n';
	return;
}

void printHelp() {
	std::cout << "Incorrect input.\n"
	          << "Program mode:\n" 
	          << "- filename\n"
	          << "- filename debug_flag(+)\n"
	          << "- Nx Ny K1 K2\n"
	          << "- Nx Ny K1 K2 debug_flag(+)\n"
	          << "File with name filename must contain parameters Nx Ny K1 K2.\n"
	          << "Nx, Ny, K1, K2 have integer type. 1 <= Nx * Ny <= 10^7.\n";
}

int main(int argc, char* argv[]) {
	int Nx, Ny, K1, K2;
	if (argc == 2 or argc == 3) {
		std::ifstream in_file(argv[1]);
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
		in_file >> Nx >> Ny >> K1 >> K2;
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
	} else if (argc == 5 or argc == 6) {
		try {
			Nx = std::stoi(argv[1]);
			Ny = std::stoi(argv[2]);
			K1 = std::stoi(argv[3]);
			K2 = std::stoi(argv[4]);
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
	if (Nx <= 0 or Ny <= 0 or K1 <= 0 or K2 <= 0 or (static_cast<long long>(Nx) * Ny >= 1'000'000'000)) {
		std::cout << "Bad parameter values" << '\n';
		printHelp();
		return 1;
	}
	bool debug_flag = false;
	if (argc == 3 or argc == 6) {
		debug_flag = true;
	}

	const size_t sz_row_ptr = (Nx + 1) * (Ny + 1) + 1;
	const size_t sz_col_ptr =
	    2 * (Nx * (Ny + 1) + Ny * (Nx + 1) +                                    // horizontal, vertical
	    (Nx * Ny) / (K1 + K2) * K2 + std::max((Nx * Ny) % (K1 + K2) - K1, 0)) + // oblique
        (Nx + 1) * (Ny + 1);                                                    // with itself

    double time = -omp_get_wtime();
	std::vector<int> row_ptr(sz_row_ptr);
	std::vector<int> col_ptr(sz_col_ptr);
	generate(Nx, Ny, K1, K2, row_ptr, col_ptr);
	time += omp_get_wtime();
	std::cout << "Generation (sequential) completed.\nTime - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("seq_IA_JA.txt");
		for (int i = 0; i < static_cast<int>(row_ptr.size() - 1); ++i) {
			oFile << i << " - ";
			for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
				oFile << col_ptr[j] << ' ';
			}
			oFile << '\n';
		}
	}
	
	time = -omp_get_wtime();
	std::vector<double> A_arr(sz_col_ptr);
	std::vector<double> b_vec(sz_row_ptr - 1);
	fill(Nx, Ny, row_ptr, col_ptr, A_arr, b_vec);
	time += omp_get_wtime();
	std::cout << "Filling (sequential) completed.\nTime - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("seq_A_b.txt");
		for (int i = 0; i < static_cast<int>(row_ptr.size() - 1); ++i) {
			oFile << i << " - ";
			for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
				oFile << A_arr[j] << ' ';
			}
			oFile << " | " << b_vec[i];
			oFile << '\n';
		}
	}
/*
	// kernels test
	size_t vec_size = b_vec.size();
	std::vector<double> a(vec_size, 1.), b(vec_size, 2);

	std::cout << "dot - " << dotKernel(a, b) << '\n';
	axpbyKernel(2, a, 1, b);
	std::cout << "axpby: sum - " << std::accumulate(a.cbegin(), a.cend(), 0.)
	          << " | l2_norm - " << std::sqrt(dotKernel(a, a)) << '\n';
	// for (int i = 0; i < a.size(); ++i) {
	// 	a[i] = i;
	// }
	std::vector<double> res(vec_size, 0);
	SpMVKernel(A_arr, row_ptr, col_ptr, a, res);
	std::cout << "SpMV: sum - " << std::accumulate(res.cbegin(), res.cend(), 0.)
	          << " | l2_norm - " << std::sqrt(dotKernel(res, res)) << '\n';
	// for (int i = 0; i < a.size(); ++i) {
	// 	std::cout << res[i] << ' ';//a[i] = i;
	// }
	// std::cout << '\n';
*/

	time = -omp_get_wtime();
	std::vector<double> x_vec(b_vec.size(), 0);
	Solve(A_arr, b_vec, row_ptr, col_ptr, x_vec);
	time += omp_get_wtime();
	std::cout << "Solve (sequential) completed.\nTime - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("seq_X.txt");
		for (int i = 0; i < static_cast<int>(x_vec.size()); ++i) {
			oFile << x_vec[i] << ' ';
		}
		oFile << '\n';
	}
	
	return 0;
}
