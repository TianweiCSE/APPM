#include "Main.h"

#ifdef isWindows
std::string input_dir = "..\\";
#endif
#ifdef isLinux
std::string input_dir = "../";
#endif

int main(int argc, char** argv) {
	std::cout << "***********************" << std::endl;
	std::cout << "*        APPM         *" << std::endl;
	std::cout << "***********************" << std::endl;
	//omp_set_num_threads(4);
	//Eigen::setNbThreads(12);
	std::cout << "threads number = " << Eigen::nbThreads() << std::endl;
	if (argc == 2) {
		input_dir = std::string(argv[1]);
	}
	else if (argc > 2) {
		std::cout << "Too much command line argument!" << std::endl;
		return EXIT_FAILURE;
	}
	Main main;
	main.run();
	std::cout << "TERMINATED" << std::endl;
	return EXIT_SUCCESS;
}


Main::Main()
{
}


Main::~Main()
{
}

void Main::run()
{	
	PrimalMesh::PrimalMeshParams primalMeshParams;
	primalMeshParams = PrimalMesh::PrimalMeshParams(input_dir + "primalMeshParams.txt");

	AppmSolver appm(primalMeshParams);
	appm.run();	
}
