#include "funcs.h"

// columns = x
// rows    = y - and it's inversed (points with greater y will be lower)

// boundary and initial values for mesh of solution

double g(double x, double y) {
	return 0;
}

 // the function in RHS

double f(double x, double y) {
	if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.2*0.2) {
		return -0.5;
	}
	else {
		return 0;
	}
}

void output_mesh(double** mesh, int n, double elapsed)
{
	std::ofstream f(OUTPUT_FILE);

	f << n << std::endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			f << mesh[i][j] << " ";
		}
		f << std::endl;
	}

	f.close();
}

void calculate_parallel(int n, int num_t)
{
	omp_set_num_threads(num_t);

	// n - 2 = number of my rows, 2 rows are constant as a boundary condition
	int delta = (n - 2)/num_t + 1;
	if (n % num_t == 0) delta--;

	double start, end;

	double** mesh = new double*[n];
	double** f_m  = new double*[n];

	// make meshes

	for (int i = 0; i < n; ++i) {
		mesh[i] = new double[n];

		for (int j = 0; j < n; ++j) {
			mesh[i][j] = g((1.0*j)/(n - 1), (i*1.0)/(n - 1));
		}
	}

	for (int i = 0; i < n; ++i) {
		f_m[i] = new double[n];

		for (int j = 0; j < n; ++j) {
			f_m[i][j] = f((1.0*j)/(n - 1), (i*1.0)/(n - 1));
		}
	}

	start = omp_get_wtime();

	#pragma omp parallel 
	{
		int tid = omp_get_thread_num();

		// the begin and the end of the group of rows which will be updated by this thread
		int i_start = std::min(1 + delta*(tid), n - 2); 
		int i_end   = std::min(1 + delta*(tid + 1) - 1, n - 2);

		int chunck_size = i_end - i_start + 1 + 2; // + 2 rows in the top and in the bottom, which will not be updated
		// in this chunck; we need them to calculate the rows of elements which are close to them

		double** tchunck1 = new double*[chunck_size];
		double** tchunck2 = new double*[chunck_size];

		// init only inner rows, which will be updated

		for (int i = 0; i <= chunck_size - 1; ++i) {
			tchunck1[i] = new double[n];
			tchunck2[i] = new double[n];

			for (int j = 0; j < n; ++j) {
				tchunck1[i][j] = mesh[i_start + i - 1][j];
				tchunck2[i][j] = mesh[i_start + i - 1][j];
			}
		}

		// flag of current chunck
		double flag = 2;
		double** tchunck;
		double** tchunck_prev;

		for (int ITER = 1; ITER <= MAX_ITER; ++ITER) {

			// change previous mesh and new one

			if (flag == 1) {
				flag = 2;
			}
			else {
				flag = 1;
			}

			// choose current mesh

			if (flag == 1) {
				tchunck = tchunck1;
				tchunck_prev = tchunck2;
			}
			else {
				tchunck = tchunck2;
				tchunck_prev = tchunck1;
			}

			// actualize top and bottom rows, which was updated by other two threads

			for (int j = 1; j <= n - 2; ++j) {
				tchunck_prev[0][j] = mesh[i_start - 1][j];
				tchunck_prev[chunck_size - 1][j] = mesh[i_end + 1][j];

				tchunck[0][j] = mesh[i_start - 1][j];
				tchunck[chunck_size - 1][j] = mesh[i_end + 1][j];
			}
			
			// update inner rows in chunck

			for (int i = 1; i <= chunck_size - 2; ++i) {
				for (int j = 1; j <= n - 2; ++j) {
					tchunck[i][j] = 0.25*(tchunck_prev[i+1][j] + tchunck_prev[i][j+1] + tchunck_prev[i-1][j] + tchunck_prev[i][j-1] - f_m[i_start + i -1][j]);
				}
			}

			// save new borders with other chuncks

			for (int j = 1; j <= n - 2; ++j) {
				mesh[i_start][j] = tchunck[1][j];
				mesh[i_end][j] = tchunck[chunck_size-2][j];
			}

			#pragma omp barrier
		}

		// throw all calculated chunck of values to to the mesh

		for (int i = 1; i <= chunck_size - 2; ++i) {
			if (flag == 1) {
				for (int j = 1; j <= n - 2; ++j) {
					mesh[i_start + i - 1][j] = tchunck1[i][j];
				}
			}
			else {
				for (int j = 1; j <= n - 2; ++j) {
					mesh[i_start + i - 1][j] = tchunck2[i][j];
				}				
			}
		}
	}

	end = omp_get_wtime();

	output_mesh(mesh, n, end - start);
}