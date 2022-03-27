#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <chrono>
using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main(int argc, char** argv){

	//read data from file
	if(argc!=2){
		std::cerr << "Error: expected one and only one argument.\n";
		return -1;
	}
	char* filename = argv[1];
	FILE *fs = fopen(filename, "r");
	int n = 0;
	fscanf(fs, "%d", &n);
	cout << "N: " << n << endl;
	int *IA = (int*)malloc(sizeof(int) * (n + 1));
	int i = 0, j = 0;
	for (i = 0; i < n + 1; i++){
		fscanf(fs, "%d", &IA[i]);
		IA[i]--;
	}
	int nnz = IA[n]; //уже уменьшила на 1
	int *JA = (int*)malloc(sizeof(int) * nnz);
	for (i = 0; i < nnz; i++){
		fscanf(fs, "%d", &JA[i]);
		JA[i]--;
	}
	double *A = (double*)malloc(sizeof(double) * nnz);
	for (i = 0; i < nnz; i++){
		fscanf(fs, "%lf", &A[i]);
	}

	// convert data into Eigen SparseMatrix format
	vector<T> tripletList;
	tripletList.reserve(n);
	int num_elems_row = 0, k = 0;
	double v = 0;
	for (i = 0; i < n; i++){
		// i - current row number
		num_elems_row = IA[i + 1] - IA[i];
		for (k = 0; k < num_elems_row; k++){
			j = JA[IA[i] + k];
			v = A[IA[i] + k];
			tripletList.push_back(T(i, j, v));
		}
	}
	SparseMatrix<double> M(n, n);
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	

	// solve the system
	VectorXd b(n), x(n);
	for (i = 0; i < n; i++){
		b[i] = 1;
		x[i] = -1;
	}
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
	//initialization: columns reordering and A factorization
	auto t1 = high_resolution_clock::now();
	solver.analyzePattern(M); 
	// Compute the numerical factorization 
	solver.factorize(M); 
	auto t2 = high_resolution_clock::now();
	auto t = duration_cast<milliseconds>(t2 - t1);
	cout << "Initialization time: " << t.count() << " ms" << endl;
	// Use the factors to solve the linear system 
	t1 = high_resolution_clock::now();
	x = solver.solve(b);
	t2 = high_resolution_clock::now();
	t = duration_cast<milliseconds>(t2 - t1);
	cout << "Execution time: " << t.count() << " ms" << endl;
	// calculate l2-norm of the residual
	VectorXd res(n);
	res = M * x - b;
	cout << "Residual norm: " << res.norm() << endl;
	free(IA);
	free(JA);
	free(A);
	return 0;





}	
