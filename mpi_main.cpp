#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <set>

#include <mpi.h>

double dotKernel(const int N_own,
                 const std::vector<double>& vec_1,
                 const std::vector<double>& vec_2) {
	double local_sum = 0;
	for (int i = 0; i < N_own; ++i) {
		local_sum += vec_1[i] * vec_2[i];
	}
	double global_sum;
	MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return global_sum;
}

void axpbyKernel(const double a, std::vector<double>& x,
	             const double b, const std::vector<double>& y) {
	std::transform(x.cbegin(), x.cend(),
	               y.cbegin(), x.begin(),
	               [a, b] (const double x_i, const double y_i) {
	                   return a * x_i + b * y_i;
	                }
	);
}

struct Com_struct {
	std::vector<int> neighbours;
	std::vector<int> send_offset = {0},
	                 recv_offset = {0};
	std::vector<int> send, recv;
};

void update(std::vector<double>& vec, const Com_struct& com, const std::vector<std::vector<MPI_Datatype>>& types) {
	static const int num_of_neib_proc = static_cast<int>(com.neighbours.size());
	static std::vector<MPI_Status> status_vec(2 * num_of_neib_proc);
	static std::vector<MPI_Request> request_vec(2 * num_of_neib_proc, MPI_REQUEST_NULL);

	for (int p = 0; p < num_of_neib_proc; ++p) {
		int neib_idx = com.neighbours[p];
		MPI_Irecv(vec.data(), 1, types[p][1], neib_idx,
		          0, MPI_COMM_WORLD, request_vec.data() + p);
		MPI_Isend(vec.data(), 1, types[p][0], neib_idx,
		          0, MPI_COMM_WORLD, request_vec.data() + p + num_of_neib_proc);
	}

	MPI_Waitall(num_of_neib_proc, request_vec.data(), status_vec.data());
}

void SpMVKernel(const int N_own,
                const std::vector<double>& matrix, const std::vector<int>& row_ptr,
                const std::vector<int>& col_ptr, std::vector<double>& vec,
                const Com_struct& com, const std::vector<std::vector<MPI_Datatype>>& types,
                std::vector<double>& res) {
	update(vec, com, types);
	for (int i = 0; i < N_own; ++i) {
		double sum = 0;
		for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
			sum += matrix[j] * vec[col_ptr[j]];
		}
		res[i] = sum;
	}
	return;
}

bool hasOblique(const int K1, const int K2, const int cell_idx) {
	return cell_idx / (K1 + K2) * (K1 + K2) + K1 <= cell_idx;
}

void get_local_begin_end(const int N, const int p, const int idx, int& b, int& e) {
	b = N / p * idx;
	e = N / p * (idx + 1);
	if (idx < N % p) {
		b += idx;
	} else {
		b += N % p;
	}
	if ((idx + 1) < N % p) {
		e += idx + 1;
	} else {
		e += N % p;
	}
}

void generate(const int Nx, const int Ny,
	          const int K1, const int K2,
	          const int Px, const int Py,
	          const int rank, const int size,
	          int& N, int& N_own, 
	          std::vector<int>& row_ptr,
	          std::vector<int>& col_ptr,
	          std::vector<int>& Part,
	          std::vector<int>& L2G,
	          std::vector<int>& G2L) {

	const int Pi = rank / Px;
	const int Pj = rank % Px;
	int ib, ie, jb, je;
	get_local_begin_end(Ny + 1, Py, Pi, ib, ie);
	get_local_begin_end(Nx + 1, Px, Pj, jb, je);
	int N_halo = 0;

	// own
	for (int i = ib; i < ie; ++i) {
		for (int j = jb; j < je; ++j) {
			L2G.push_back(i * (Nx + 1) + j);
			Part.push_back(rank);
		}
	}

	//halo
	if (Pi != 0) {
		int i = ib - 1;
		for (int j = jb; j < je; ++j) {
			L2G.push_back(i * (Nx + 1) + j);
			++N_halo;
			Part.push_back((Pi - 1) * Px + Pj);
		}
	}
	if (Pi != Py - 1) {
		int i = ie;
		for (int j = jb; j < je; ++j) {
			L2G.push_back(i * (Nx + 1) + j);
			++N_halo;
			Part.push_back((Pi + 1) * Px + Pj);
		}
	}
	if (Pj != 0) {
		int j = jb - 1;
		for (int i = ib; i < ie; ++i) {
			L2G.push_back(i * (Nx + 1) + j);
			++N_halo;
			Part.push_back(Pi * Px + (Pj - 1));
		}
	}
	if (Pj != Px - 1) {
		int j = je;
		for (int i = ib; i < ie; ++i) {
			L2G.push_back(i * (Nx + 1) + j);
			++N_halo;
			Part.push_back(Pi * Px + (Pj + 1));
		}
	}
	if (Pi != 0 and Pj != 0 and hasOblique(K1, K2, (ib - 1) * Nx + (jb - 1))) {
		L2G.push_back((ib - 1) * (Nx + 1) + (jb - 1));
		++N_halo;
		Part.push_back((Pi - 1) * Px + (Pj - 1));
	}
	if (Pi != Py - 1 and Pj != Px - 1 and hasOblique(K1, K2, (ie - 1) * Nx + (je - 1))) {
		L2G.push_back(ie * (Nx + 1) + je);
		++N_halo;
		Part.push_back((Pi + 1) * Px + (Pj + 1));
	}

	N = static_cast<int>(L2G.size());
	N_own = N - N_halo;
	G2L.resize((Nx + 1) * (Ny + 1), -1);
	for (int i = 0; i < N; ++i) {
		G2L[L2G[i]] = i;
	}

	row_ptr.push_back(0);
	for (int i = ib; i < ie; ++i) {
		for (int j = jb; j < je; ++j) {
			int neib_num = 0;
			if (i != 0 and j != 0 and hasOblique(K1, K2, Nx * (i - 1) + (j - 1))) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * (i - 1) + (j - 1)]);
			}
			if (i != 0) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * (i - 1) + j]);
			}
			if (j != 0) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * i + (j - 1)]);
			}
			++neib_num;
			col_ptr.push_back(G2L[(Nx + 1) * i + j]);
			if (j != Nx) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * i + (j + 1)]);
			}
			if (i != Ny) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * (i + 1) + j]);
			}
			if (i != Ny and j != Nx and hasOblique(K1, K2, Nx * i + j)) {
				++neib_num;
				col_ptr.push_back(G2L[(Nx + 1) * (i + 1) + (j + 1)]);
			}
			row_ptr.push_back(row_ptr.back() + neib_num);
		}
	}
}

void fill(const int N, const int N_own,
          const std::vector<int>& row_ptr,
          const std::vector<int>& col_ptr,
          const std::vector<int>& L2G,
          std::vector<double>& A_arr,
          std::vector<double>& b_vec) {
	for (int cur_idx = 0; cur_idx < N_own; ++cur_idx) {
		double sum = 0;
		int diag_idx;
		for (int k = row_ptr[cur_idx]; k < row_ptr[cur_idx + 1]; ++k) {
			const int neib_idx = col_ptr[k];
			if (neib_idx == cur_idx) {
				diag_idx = k;
			} else {
				A_arr[k] = std::cos(static_cast<double>(L2G[cur_idx]) * L2G[neib_idx] + L2G[cur_idx] + L2G[neib_idx]);
				sum += std::abs(A_arr[k]);
			}
		}
		A_arr[diag_idx] = sum * 1.234;
	}
	for (int cur_idx = 0; cur_idx < N; ++cur_idx) {
		b_vec[cur_idx] = std::sin(L2G[cur_idx]);
	}
}

void Com(const int rank, const int size,
         const int N, const int N_own,
         const std::vector<int>& row_ptr,
         const std::vector<int>& col_ptr,
         const std::vector<int>& Part,
         const std::vector<int>& L2G,
         const std::vector<int>& G2L,
         Com_struct& com) {

	const auto comp_for_set =
	    [&L2G] (const int l, const int r) {
	        return L2G[l] < L2G[r];
	    };
	std::vector<std::set<int, decltype(comp_for_set)>>
	    send_to_process(size, std::set<int, decltype(comp_for_set)>(comp_for_set)),
	    recv_from_process(size, std::set<int, decltype(comp_for_set)>(comp_for_set));

	for (int row_idx = 0; row_idx < N_own; ++row_idx) {
		for (int k = row_ptr[row_idx]; k < row_ptr[row_idx + 1]; ++k) {
			const int col_idx = col_ptr[k];
			const int proc_rank = Part[col_idx];
			if (proc_rank != rank) {
				send_to_process[proc_rank].insert(row_idx);
				recv_from_process[proc_rank].insert(col_idx);
			}
		}
	}
	for (int neib_idx = 0; neib_idx < size; ++neib_idx) {
		if (not send_to_process[neib_idx].empty()) {
			com.neighbours.push_back(neib_idx);
			std::copy(send_to_process[neib_idx].cbegin(),
			          send_to_process[neib_idx].cend(),
			          std::back_inserter(com.send));
			std::copy(recv_from_process[neib_idx].cbegin(),
			          recv_from_process[neib_idx].cend(),
			          std::back_inserter(com.recv));
			com.send_offset.push_back(static_cast<int>(com.send.size()));
			com.recv_offset.push_back(static_cast<int>(com.recv.size()));
		}
	}
}

void Solve(const int N, const int N_own, const int rank,
	       const std::vector<double>& A_arr, const std::vector<double>& b_vec,
	       const std::vector<int>& row_ptr, const std::vector<int>& col_ptr,
	       const Com_struct& com, const std::vector<std::vector<MPI_Datatype>>& types,
	       std::vector<double>& x_vec, const double TOL = 1e-7) {
	       
	std::vector<double> p_vec(N), z_vec(N),
	                    q_vec(N), r_vec(b_vec),
	                    inverse_M;
	std::vector<int> M_row_ptr = {0}, M_col_ptr;
	for (int i = 0; i < N_own; ++i) {
		int diag_idx;
		for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
			if (col_ptr[k] == i) {
				diag_idx = k;
				break;
			}
		}
		inverse_M.push_back(1. / A_arr[diag_idx]);
		M_row_ptr.push_back(M_row_ptr.back() + 1);
		M_col_ptr.push_back(i);
	}
	bool do_cycle = true;
	int it_num = 1;
	const int MAX_IT_NUM = 100;
	double dot_time = 0., axpby_time = 0., spmv_time = 0.;
	dot_time -= MPI_Wtime();
	const double b_norm = dotKernel(N_own, b_vec, b_vec);
	dot_time += MPI_Wtime();
	const double EPS = TOL * b_norm;
	double alpha, rho_prev, rho_cur, beta, residual_norm;
	do {
		spmv_time -= MPI_Wtime();
		SpMVKernel(N_own, inverse_M, M_row_ptr, M_col_ptr, r_vec, com, types, z_vec);
		spmv_time += MPI_Wtime();
		dot_time -= MPI_Wtime();
		rho_cur = dotKernel(N_own, r_vec, z_vec);
		dot_time += MPI_Wtime();
		if (it_num == 1) {
			p_vec = z_vec;
		} else {
			beta = rho_cur / rho_prev;
			axpby_time -= MPI_Wtime();
			axpbyKernel(beta, p_vec, 1., z_vec);
			axpby_time += MPI_Wtime();
		}
		spmv_time -= MPI_Wtime();
		SpMVKernel(N_own, A_arr, row_ptr, col_ptr, p_vec, com, types, q_vec);
		spmv_time += MPI_Wtime();
		dot_time -= MPI_Wtime();
		alpha = rho_cur / dotKernel(N_own, p_vec, q_vec);
		dot_time += MPI_Wtime();
		axpby_time -= MPI_Wtime();
		axpbyKernel(1., x_vec,  alpha, p_vec);
		axpbyKernel(1., r_vec, -alpha, q_vec);
		axpby_time += MPI_Wtime();
		dot_time -= MPI_Wtime();
		residual_norm = dotKernel(N_own, r_vec, r_vec);
		dot_time += MPI_Wtime();
/*
		if (rank == 0) {
			std::cout << "It - " << it_num <<
			             ", RES_NORM - " << residual_norm <<
			             ", tol - " << residual_norm / b_norm << '\n';
		}
*/
		if (residual_norm < EPS or it_num >= MAX_IT_NUM) {
			do_cycle = false;
		} else {	
			++it_num;
			rho_prev = rho_cur;
		}
	} while (do_cycle);
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &dot_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &axpby_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &spmv_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		std::cout << "Dot func time - " << dot_time <<
		             "\nAxpby func time - " << axpby_time <<
		             "\nSpmv func time - " << spmv_time << '\n';
	} else {
		MPI_Reduce(&dot_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&axpby_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&spmv_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}
}

void printHelp() {
	std::cout << "Incorrect input.\n"
	          << "Program mode:\n" 
	          << "- filename\n"
	          << "- Nx Ny K1 K2 Px Py\n"
	          << "File with name filename must contain parameters Nx Ny K1 K2 Px Py.\n"
	          << "Nx, Ny, K1, K2, Px, Py have integer type. Nx, Ny, K1, K2, Px, Py > 0, 1 <= Nx * Ny <= 10^7.\n";
}

void create_types(const Com_struct& com, std::vector<std::vector<MPI_Datatype>> &types) {
	const int num_of_neib_proc = static_cast<int>(com.neighbours.size());
	for (int p = 0; p < num_of_neib_proc; ++p) {
		int message_size = com.send_offset[p + 1] - com.send_offset[p];
		MPI_Type_indexed(message_size, std::vector<int>(message_size, 1).data(),
		                 com.send.data() + com.send_offset[p], MPI_DOUBLE, &types[p][0]);

		message_size = com.recv_offset[p + 1] - com.recv_offset[p];
		MPI_Type_indexed(message_size, std::vector<int>(message_size, 1).data(),
		                 com.recv.data() + com.recv_offset[p], MPI_DOUBLE, &types[p][1]);
		MPI_Type_commit(&types[p][0]);
		MPI_Type_commit(&types[p][1]);
	}
}

int main(int argc, char* argv[]) {
	int Nx, Ny, K1, K2, Px, Py;
	if (argc == 2) {
		std::ifstream in_file(argv[1]);
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
		in_file >> Nx >> Ny >> K1 >> K2 >> Px >> Py;
		if (not in_file.good()) {
			printHelp();
			return 1;
		}
	} else if (argc == 7) {
		try {
			Nx = std::stoi(argv[1]);
			Ny = std::stoi(argv[2]);
			K1 = std::stoi(argv[3]);
			K2 = std::stoi(argv[4]);
			Px = std::stoi(argv[5]);
			Py = std::stoi(argv[6]);
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
	if (Nx <= 0 or Ny <= 0 or
	    K1 <= 0 or K2 <= 0 or
	    Px <= 0 or Py <= 0 or
	    (static_cast<long long>(Nx) * Ny >= 1000000000)) {
		std::cout << "Bad parameter values" << '\n';
		printHelp();
		return 1;
	}

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	double stage_time;
// 0
	int N, N_own;
	std::vector<int> row_ptr, col_ptr, Part, L2G, G2L;
	stage_time = -MPI_Wtime();
	generate(Nx, Ny, K1, K2, Px, Py, rank, size,
	        N, N_own, row_ptr, col_ptr, Part, L2G, G2L);
	stage_time += MPI_Wtime();
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &stage_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		std::cout << "Gen(mpi). Time - " << stage_time << '\n';
	} else {
		MPI_Reduce(&stage_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}

// 1
	std::vector<double> A_arr(col_ptr.size()), b_vec(N);
	stage_time = -MPI_Wtime();
	fill(N, N_own, row_ptr, col_ptr,
	     L2G, A_arr, b_vec);
	stage_time += MPI_Wtime();
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &stage_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		std::cout << "Fill(mpi). Time - " << stage_time << '\n';
	} else {
		MPI_Reduce(&stage_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}

// 2
	Com_struct com;
	Com(rank, size, N, N_own,
	    row_ptr, col_ptr,
	    Part, L2G, G2L,
	    com);


// 3
	std::vector<double> x_vec(N, 0);
	std::vector<std::vector<MPI_Datatype>> types(com.neighbours.size(), std::vector<MPI_Datatype>(2));
	create_types(com, types);
	stage_time = -MPI_Wtime();
	Solve(N, N_own, rank, A_arr,  b_vec,
	      row_ptr, col_ptr, com, types,
	      x_vec);
	stage_time += MPI_Wtime();
	if (rank == 0) {
		MPI_Reduce(MPI_IN_PLACE, &stage_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		std::cout << "Solve(mpi). Time - " << stage_time << '\n';
	} else {
		MPI_Reduce(&stage_time, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

