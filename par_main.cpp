#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <omp.h>

int obliquesBefore(const int K1, const int K2, const int cell_idx) {
	return cell_idx / (K1 + K2) * K2 + std::max(cell_idx % (K1 + K2) - K1, 0);
}

bool hasOblique(const int K1, const int K2, const int cell_idx) {
	return cell_idx / (K1 + K2) * (K1 + K2) + K1 <= cell_idx;
}

void generate(const int Nx, const int Ny, const int K1, const int K2,
	          std::vector<int>& row_ptr, std::vector<int>& col_ptr) {
	row_ptr[0] = 0;
#pragma omp parallel for
	for (int i = 0; i < (Ny + 1); ++i) {
		size_t my_idx = 0;
			//vertical
			if (i > 1) {
				my_idx += 2 * (Nx + 1) * (i - 1);	
			}
			if (i > 0) {
				//my_idx += 2 * j + (Nx + 1 - j);
				my_idx += Nx + 1;
			}
			// if (i != Ny) {
			// 	my_idx += j;	
			// }
			// horizontal
			my_idx += 2 * Nx * i;
			// if (j > 0) {
			// 	my_idx += 1;
			// 	my_idx += 2 * (j - 1);
			// }
			// oblique
			// if (j != 0 and i != Ny) {
			// 	int tmp = Nx * i + j;
			// 	my_idx += obliquesBefore(K1, K2, tmp);
			// 	if (i != 0) {
			// 		tmp = Nx * i;
			// 		my_idx -= obliquesBefore(K1, K2, tmp);
			// 	}
			// }
			if (i != 0) {
				int tmp = Nx * i;
				tmp = obliquesBefore(K1, K2, tmp);
				my_idx += 2 * tmp;
				int ptrL, ptrB = Nx * i;
				//if (j == 0) {
					ptrL = ptrB - Nx;
				//} else {
				//	ptrL = Nx * (i - 1) + (j - 1);
				//}
				my_idx -= obliquesBefore(K1, K2, ptrB) - obliquesBefore(K1, K2, ptrL);
			}
			my_idx += (Nx + 1) * i;
		for (int j = 0; j < (Nx + 1); ++j) {
			const int idx_begin = static_cast<int>(my_idx);
			int neib_num = 0;
			if (i != 0 and j != 0 and hasOblique(K1, K2, Nx * (i - 1) + (j - 1))) {	
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * (i - 1) + (j - 1);
				++my_idx;
			}
			if (i != 0) {
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * (i - 1) + j;
				++my_idx;
			}
			if (j != 0) {
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * i + (j - 1);
				++my_idx;
			}
			++neib_num;
			col_ptr[my_idx] = (Nx + 1) * i + j;
			++my_idx;
			if (j != Nx) {
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * i + (j + 1);
				++my_idx;
			}
			if (i != Ny) {
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * (i + 1) + j;
				++my_idx;
			}
			if (i != Ny and j != Nx and hasOblique(K1, K2, Nx * i + j)) {
				++neib_num;
				col_ptr[my_idx] = (Nx + 1) * (i + 1) + (j + 1);
				++my_idx;
			}
			row_ptr[(Nx + 1) * i + j + 1] = idx_begin + neib_num;
		}
	}
	return;
}

void fill(const int Nx, const int Ny,
	     const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	     std::vector<double>& A_arr, std::vector<double>& b_vec) {
#pragma omp parallel for
	for (int i = 0; i < (Ny + 1); ++i) {
		for (int j = 0; j < (Nx + 1); ++j) {
			const int cur_idx = (Nx + 1) * i + j;
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
	}
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
    std::cout << "Graph vertex - " << sz_row_ptr - 1 << '\n';
    std::cout << "Graph edges - " << (sz_col_ptr - (Nx + 1) * (Ny + 1)) / 2 << '\n';
    std::cout << "Non-zero elements - " << sz_col_ptr << '\n';

	std::vector<int> row_ptr(sz_row_ptr);
	std::vector<int> col_ptr(sz_col_ptr);
	double time = -omp_get_wtime();
	generate(Nx, Ny, K1, K2, row_ptr, col_ptr);
	time += omp_get_wtime();
	std::cout << "Generation (parallel) completed.\nTime - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("parallel_IA_JA.txt");
		for (int i = 0; i < row_ptr.size() - 1; ++i) {
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
	std::cout << "Filling (parallel) completed.\nTime - " << time << '\n';
	if (debug_flag) {
		std::ofstream oFile ("parallel_A_b.txt");
		for (int i = 0; i < row_ptr.size() - 1; ++i) {
			oFile << i << " - ";
			for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
				oFile << A_arr[j] << ' ';
			}
			oFile << " | " << b_vec[i];
			oFile << '\n';
		}
	}
	return 0;
}
