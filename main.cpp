#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

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
			Ny = std::stoi(atgv[2]);
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
	bool debug_flag = false;
	if (argc == 3 or argc == 6) {
		debug_flag = true;
	}




	return 0;
}