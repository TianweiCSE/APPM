#include "Main.h"

int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*        APPM         *" << std::endl;
	std::cout << "***********************" << std::endl;
	//omp_set_num_threads(4);
	//Eigen::setNbThreads(12);
	std::cout << "threads number = " << Eigen::nbThreads() << std::endl;
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
	primalMeshParams = PrimalMesh::PrimalMeshParams("primalMeshParams.txt");

	AppmSolver appm(primalMeshParams);
	appm.run();	
}
