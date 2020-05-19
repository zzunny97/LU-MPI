#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <cmath>
#include <sys/time.h>

using namespace std;
typedef double Data;
struct timeval start_point, end_point;
double optime;

int N;							// final matrix size
int original_N;				    // original matrix size
Data** A_original;				// copied A, which is original
int rank;						// # of process
int process_size;				// total process number

void Print_mat(Data** p, int size) {
    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            printf("%.3f\t", p[i][j]);
        }
        cout << endl;
    }
    cout << endl;
}

void Print_mat2(Data** p, int size1, int size2) {
	for(int i=0; i<size1; i++) {
		for(int j=0; j<size2; j++) {
			printf("%.3f\t\t", p[i][j]);
		}
		cout << endl;
	}
	cout << endl;
}

void Copy_mat(Data** src, Data** dst, int size) {
    for(int i=0; i<size; i++) 
        memcpy(dst[i], src[i], sizeof(Data)* size);
}

void Free_mat(Data** p, int size) {
    for(int i=0; i<size; i++)
        delete[] p[i];
    delete[] p;
}

void Inverse_mat(Data** p, int size) {
    if(size < 2)
        return;
	
	// use ad-bc, if size == 2
    else if(size == 2) {
        Data prefix = 1 / (p[0][0]*p[1][1] - p[0][1]*p[1][0]);
        for(int i=0; i<size; i++) {
            for(int j=0; j<size; j++) {
                p[i][j] *= prefix;
            }
        }
        Data tmp = p[1][1];
        p[1][1] = p[0][0];
        p[0][0] = tmp;
        p[0][1] *= -1;
        p[1][0] *= -1;
        return;
    }

	Data tmp;
    Data** mat = new Data*[size];
    for(int i=0; i<size; i++) {
        mat[i] = new Data[2*size];
		memset(mat[i], 0, sizeof(Data)*2*size);
	}

    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            mat[i][j] = p[i][j];
        }
    }


    for(int i=0; i<size; i++) {
        for(int j=size; j< 2*size; j++) {
            if(j==(i+size)) {
				mat[i][j] = 1;
			}
			// else {
			// 	mat[i][j] = 0;
			// }
        }
    }


    // for(int i=size-1; i>0; i--) {
    //     if(mat[i-1][0] < mat[i][0]) {
	// 		Data* tmp2;
	// 		/*
    //         for(int j=0; j<2*size; j++) {
    //             tmp2 = mat[i];
    //             mat[i] = mat[i-1];
    //             mat[i-1] = tmp2;
    //         }*/
	// 		tmp2 = mat[i];
	// 		mat[i] = mat[i-1];
	// 		mat[i-1] = tmp2;
	// 	}
    // }
    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            if(j!=i) {
                tmp = mat[j][i] / mat[i][i];
                for(int k=0; k<size*2; k++)
                    mat[j][k] -= mat[i][k]*tmp;
            }
        }
    }
	for(int i=0; i<size; i++) {
		tmp = mat[i][i];
		for(int j=0; j<2*size; j++) {
			mat[i][j] /= tmp;

		}
	}
    for(int i=0; i<size; i++) 
        for(int j=size; j<size*2; j++) 
            p[i][j-size] = mat[i][j];
}

void verify(Data** A, Data** L, Data** U, int size) {
    double errsum = 0.0;
	Data norm = 0.0, tmp = 0.0;
	// calculate errsum, Euclidean norm
    for(int i=0; i<size; i++) {
		norm = 0;
        for(int j=0; j<size; j++) {
			tmp = 0;
            for(int k=0; k<size; k++) {
                // C[i][j] += L[i][k] * U[k][j];
				tmp += L[i][k] * U[k][j];
			}
			norm += pow(A[i][j] - tmp, 2);
		}
		printf("Errsum %d: %.9lf\n", i, errsum);
		errsum += sqrt(norm);
	}
    printf("Errsum: %.9lf\n", errsum);
}

void verify2(Data** A, Data** L, Data** U, int size) {
    double errsum = 0.0;
	Data norm = 0.0, tmp = 0.0;
	// calculate errsum, Euclidean norm
    for(int i=0; i<size; i++) {
		norm = 0;
        for(int j=0; j<size; j++) {
			tmp = 0;
            for(int k=0; k<size; k++) {
                // C[i][j] += L[i][k] * U[k][j];
				tmp += L[i][k] * U[k][j];
			}
			norm += (A[i][j] - tmp);
		}
		printf("Errsum: %.9lf\n", errsum);
		errsum += norm;
	}
    printf("Errsum: %.9lf\n", errsum);
}



void Mul_mat(Data** a, Data** b, Data** ret, int size) {
    for(int i=0; i<size; i++) {
        memset(ret[i], 0, sizeof(Data) * size);
    }
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++) 
                ret[i][j] += a[i][k] * b[k][j];
}

void Sub_mat(Data** a, Data** b, int size) {
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			a[i][j] -= b[i][j];
}


void block_LU(Data** a, Data** l, Data** u, int size) {
    for(int k=0; k<size; k++) {
        u[k][k] = a[k][k];
        for(int i=k+1; i<size; i++) {
            l[i][k] = a[i][k] / u[k][k];
            u[k][i] = a[k][i];
        }
        for(int i=k+1; i<size; i++)
            for(int j=k+1; j<size; j++)
                a[i][j] -= l[i][k] * u[k][j];
    }
}

void Copy_to_one(Data** src, Data* dst, int size) {
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			dst[i*size+j] = src[i][j];
}

void Copy_to_two(Data* src, Data** dst, int size) {
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			dst[i][j] = src[i*size+j];
}


void LU(Data** A, Data** L, Data** U) {
	if(rank == 0)
		cout << "Start LU" << endl;

    int block_sq = sqrt(process_size);     
    int row_per_block = (N / block_sq);   
    int x_start = (rank / block_sq) * row_per_block;    // start of x in the process
    int x_end = x_start + row_per_block;                // end of x in the process
    int y_start = (rank % block_sq) * row_per_block;    // start of y in the process
    int y_end = y_start + row_per_block;                // end of y in the process
	long oned_size = pow(row_per_block, 2);             // number of elements of a block
	long idx = 0;
   
    // ====== INITIALIZE SMALL MATRICES ==== 
    Data** small_A = new Data*[row_per_block];          // block A
    Data** small_L = new Data*[row_per_block];          // block L
    Data** small_U = new Data*[row_per_block];          // block U
    for(int i=0; i<row_per_block; i++) {
        small_A[i] = new Data[row_per_block];
        small_L[i] = new Data[row_per_block];
        small_U[i] = new Data[row_per_block];
    }
	Data* oned_mat = new Data[oned_size];               // oned_size == pow(row_per_block ,2)
	Data* oned_mat2 = new Data[oned_size];

    for(int i=0; i<row_per_block; i++) small_L[i][i] = 1;
    
    for(int i=x_start; i<x_end; i++)
        for(int j=y_start; j<y_end; j++)
            small_A[i-x_start][j-y_start] = A[i][j];

    // if the process is master
	if(rank == 0) {
		block_LU(small_A, small_L, small_U, row_per_block);	
		for(int i=0; i<row_per_block; i++)
			for(int j=0; j<row_per_block; j++)
				L[i][j] = small_L[i][j];
		for(int i=0; i<row_per_block; i++)
			for(int j=0; j<row_per_block; j++)
				U[i][j] = small_U[i][j];

        Inverse_mat(small_L, row_per_block);
        Inverse_mat(small_U, row_per_block);
		Copy_to_one(small_L, oned_mat, row_per_block);
		Copy_to_one(small_U, oned_mat2, row_per_block);
		// send to row
        for(int dst = 1; dst<block_sq; dst++) {
			MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD); // send to row the inverse of L
        }
        // send to col
        for(int dst = block_sq; dst<process_size; dst+=block_sq) {
			MPI_Send(&oned_mat2[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD); // send to col the inverse of U
        }

		for(int i=0; i<block_sq; i++){
			int center = (i * block_sq) + i;
			//cout << "center: " << center << endl;
			int x_start, y_start;
			// gather from center
			if(center != 0)	{
				x_start = (center / block_sq) * row_per_block;
				y_start = (center % block_sq) * row_per_block;
				//cout << "center: " << center << " x_start: " << x_start << " y_start: " << y_start << endl;
				MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, center, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, center, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				idx = 0;
				for(int i=x_start; i<x_start+row_per_block; i++) 
					for(int j=y_start; j<y_start+row_per_block; j++)
						L[i][j] = oned_mat[idx++];
				idx = 0;
				for(int i=x_start; i<x_start+row_per_block; i++) 
					for(int j=y_start; j<y_start+row_per_block; j++)
						U[i][j] = oned_mat2[idx++];
				//cout << "gathered " << center << endl;
			}
			else if(center == process_size - 1) {
				break;
			}

			// gather from row
				for(int src=center+1; src<(i+1)*block_sq; src++) {
					x_start = (src / block_sq) * row_per_block;
					y_start = (src % block_sq) * row_per_block;
					//cout << "normal: " << src << " x_start: " << x_start << " y_start: " << y_start << endl;

					MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv L
					idx = 0;
					for(int i=x_start; i<x_start+row_per_block; i++) 
						for(int j=y_start; j<y_start+row_per_block; j++)
							U[i][j] = oned_mat[idx++];
					//cout << "gathered " << src << endl;
				}
				// gather from col
				for(int src=center+block_sq; src<process_size; src+=block_sq) {
					x_start = (src / block_sq) * row_per_block;
					y_start = (src % block_sq) * row_per_block;
					//cout << "normal: " << src << " x_start: " << x_start << " y_start: " << y_start << endl;

					MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv L
					idx = 0;
					for(int i=x_start; i<x_start+row_per_block; i++) 
						for(int j=y_start; j<y_start+row_per_block; j++)
							L[i][j] = oned_mat2[idx++];
					//cout << "gathered " << src << endl;
				}


		}
		// // gather all from other processes
		// for(int i=0; i<block_sq; i++) {
		// 	for(int j=0; j<block_sq; j++) {
		// 		int src = i*block_sq + j;
		// 		x_start = (src / block_sq) * row_per_block;
		// 		y_start = (src % block_sq) * row_per_block;
		// 		// get both
		// 		if(i==j) {
		// 			if(i==0 && j==0)
		// 				continue;
		// 			MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	// recv L
		// 			MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv U	
		// 			idx = 0;
		// 			for(int i=x_start; i<x_start+row_per_block; i++) 
		// 				for(int j=y_start; j<y_start+row_per_block; j++)
		// 					L[i][j] = oned_mat[idx++];
		// 			idx = 0;
		// 			for(int i=x_start; i<x_start+row_per_block; i++) 
		// 				for(int j=y_start; j<y_start+row_per_block; j++)
		// 					U[i][j] = oned_mat2[idx++];
		// 			cout << "gathered " << src << endl;
		// 		}
		// 		// get only L
		// 		else if(i > j) {
		// 			MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv L
		// 			idx = 0;
		// 			for(int i=x_start; i<x_start+row_per_block; i++) 
		// 				for(int j=y_start; j<y_start+row_per_block; j++)
		// 					L[i][j] = oned_mat[idx++];
		// 			cout << "gathered " << src << endl;
		// 		}
		// 		// get only U
		// 		else {
		// 			MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv L
		// 			idx = 0;
		// 			for(int i=x_start; i<x_start+row_per_block; i++) 
		// 				for(int j=y_start; j<y_start+row_per_block; j++)
		// 					U[i][j] = oned_mat2[idx++];
		// 			cout << "gathered " << src << endl;
		// 		}
		// 	}
		// }
	}

	// center processes 
	else if(x_start == y_start) {
		// initialize mul mat
		Data** mul = new Data*[row_per_block];
		for(int i=0; i<row_per_block; i++)
			mul[i] = new Data[row_per_block];

		// recv and subtract
		for(int i=0; i<rank/block_sq; i++) {
			int src1 = (rank / block_sq) * block_sq + i;	// L
			int src2 = (rank % block_sq) + (i * block_sq);	// U
			MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    // receive L
			MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, src2 ,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   // receive U
			idx = 0;
			for(int i=0; i<row_per_block; i++)
				for(int j=0; j<row_per_block; j++)
					small_L[i][j] = oned_mat[idx++];
			idx=0;
			for(int i=0; i<row_per_block; i++)
				for(int j=0; j<row_per_block; j++)
					small_U[i][j] = oned_mat2[idx++];

			// multiply two matrices from row and col, and subtract from original A
			Mul_mat(small_L, small_U, mul, row_per_block);
			Sub_mat(small_A, mul, row_per_block);	
		}

		// block LU
		for(int i=0; i<row_per_block; i++) {
			memset(small_L[i], 0, sizeof(Data) * row_per_block);
		}
		for(int i=0; i<row_per_block; i++)
			small_L[i][i] = 1;
		for(int i=0; i<row_per_block; i++) {
			memset(small_U[i], 0, sizeof(Data) * row_per_block);
		}

		block_LU(small_A, small_L, small_U, row_per_block);
			
		// send to master both L and U
		Copy_to_one(small_L, oned_mat, row_per_block);
		Copy_to_one(small_U, oned_mat2, row_per_block);
		MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);	// send final 'L' to master
		MPI_Send(&oned_mat2[0], oned_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);   // send final 'U' to master

		if(rank != process_size - 1) {
			// send to row the inverse of L and U
			Inverse_mat(small_L, row_per_block);
            Inverse_mat(small_U, row_per_block);
			Copy_to_one(small_L, oned_mat, row_per_block); 
			Copy_to_one(small_U, oned_mat2, row_per_block);
			int limit = (rank/block_sq+1) * block_sq;
			for(int dst = rank+1; dst<limit; dst++) {
				MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD); // send row the inverse of L
			}

			// send to col the inverse of U
			for(int dst = rank+block_sq; dst<process_size; dst+=block_sq) {
				MPI_Send(&oned_mat2[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD); // send col the inverse of U
			}
		}
		Free_mat(mul, row_per_block);
	}

	// rest blocks
	else {
		// initialize mul mat
		Data** mul = new Data*[row_per_block];
		for(int i=0; i<row_per_block; i++)
			mul[i] = new Data[row_per_block];

        int iter = (x_start > y_start) ? (rank%block_sq) : (rank/block_sq);
        // multiply two matrices from row and col, and subtract from original A
        for(int i=0; i<iter; i++) {
			int src1 = (rank / block_sq) * block_sq + i;	// L
			int src2 = (rank % block_sq) + (i * block_sq);	// U
            MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // recv L
            MPI_Recv(&oned_mat2[0], oned_size, MPI_DOUBLE, src2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv U
            idx = 0;
            for(int i=0; i<row_per_block; i++)
                for(int j=0; j<row_per_block; j++)
                    small_L[i][j] = oned_mat[idx++];
            idx=0;
            for(int i=0; i<row_per_block; i++)
                for(int j=0; j<row_per_block; j++)
                    small_U[i][j] = oned_mat2[idx++];

            Mul_mat(small_L, small_U, mul, row_per_block);
            Sub_mat(small_A, mul, row_per_block);   
        }
		if(x_start > y_start) {
            // Recv inverse of U from center process
			int src = (rank % block_sq) + block_sq * (rank % block_sq);
			MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv inverse U
			idx = 0;			
			for(int i=0; i<row_per_block; i++)
				for(int j=0; j<row_per_block; j++)
					small_U[i][j] = oned_mat[idx++];
            // get final L and now it's time to send to master
			Mul_mat(small_A, small_U, mul, row_per_block);
            // send final L to master
            Copy_to_one(mul, oned_mat, row_per_block);
            MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            // send final L to row
            for (int dst = rank+1; dst < (rank / block_sq + 1) * block_sq ; dst++) {
                MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD);
            }
		}
		else if (x_start < y_start) {
            // Recv inverse of L from center process
			int src = (rank / block_sq) + block_sq * ( rank / block_sq);
			MPI_Recv(&oned_mat[0], oned_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recv inverse L
			idx = 0;			
			for(int i=0; i<row_per_block; i++)
				for(int j=0; j<row_per_block; j++)
					small_L[i][j] = oned_mat[idx++];
            // get final U and now it's time to send to master
			Mul_mat(small_L, small_A, mul, row_per_block); 
            // send final U to master
            Copy_to_one(mul, oned_mat, row_per_block);
            MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            // send final U to col
            for(int dst = rank + block_sq; dst<process_size; dst+=block_sq) {
                MPI_Send(&oned_mat[0], oned_size, MPI_DOUBLE, dst, 1, MPI_COMM_WORLD);
            }
		}
		Free_mat(mul, row_per_block);
	}
		
	
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0; i<N; i++){
		MPI_Bcast(&L[i][0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&U[i][0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	if(rank == 0) {
		gettimeofday(&end_point, NULL);
		optime = (double)(end_point.tv_sec)+(double)(end_point.tv_usec)/1000000.0-(double)(start_point.tv_sec)-(double)(start_point.tv_usec)/1000000.0;
		cout << "Initialize and LU time: " << optime << " secs" << endl;
		//cout << "===== L =====" << endl;
		//Print_mat(L, N);
		//cout << "===== U =====" << endl;
		//Print_mat(U, N);

		// verify(A_original, L, U, N);
	}
	gettimeofday(&start_point, NULL);

	int row_per_process = N / process_size;		// 8000 / 64 = 125
	int my_start = rank * row_per_process;		// 0 125 250 ...
	int my_end = my_start + row_per_process;	// 125 250 375 ... 8000
	Data tmp = 0.0;
	Data norm = 0.0;
	Data recv_val = 0.0;
	Data errsum = 0.0;
	
	if(N % row_per_process != 0) {
		if(rank == process_size - 1) {
			my_end = N;
			// cout << "I am last rank: " << "start = " << my_start << " end = " << my_end << endl;
		}
	}

	if(rank == 0){
		for(int i=my_start; i<my_end; i++) {
			norm = 0;
			for(int j=0; j<N; j++) {
				tmp = 0;
				for(int k=0; k<N; k++) {
					tmp += L[i][k] * U[k][j];
				}
				norm += pow(A_original[i][j] - tmp, 2);
			}
			errsum += sqrt(norm);
		}
		for(int i=1; i<process_size; i++) {
			MPI_Recv(&recv_val, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			errsum += recv_val;
		}
	}
	else {
		for(int i=my_start; i<my_end; i++) {
			norm = 0;
			for(int j=0; j<N; j++) {
				tmp = 0;
				for(int k=0; k<N; k++) {
					tmp += L[i][k] * U[k][j];
				}
				norm += pow(A_original[i][j] - tmp, 2);
			}
			errsum += sqrt(norm);
		}
		MPI_Send(&errsum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	gettimeofday(&end_point, NULL);
	optime = (double)(end_point.tv_sec)+(double)(end_point.tv_usec)/1000000.0-(double)(start_point.tv_sec)-(double)(start_point.tv_usec)/1000000.0;

	if(rank == 0){
		printf("Errsum by parallelize: %.9lf\n", errsum);
		//cout << "Errsum: " << errsum << endl;
		cout << "Verify time: " << optime << " secs" << endl;
	}
	Free_mat(small_A, row_per_block);
	Free_mat(small_L, row_per_block);
	Free_mat(small_U, row_per_block);
	free(oned_mat);
	free(oned_mat2);
}

int main(int argc, char** argv){
    if(argc != 3 ) {
        cout << "Usage: ./proj3 N Seed" << endl;
        exit(-1);
    }
	gettimeofday(&start_point, NULL);
    N = atoi(argv[1]);
	original_N = atoi(argv[1]);
    int seed = atoi(argv[2]);
    srand(seed);
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_size);
	int sq_pr = sqrt(process_size);

	if(N % sq_pr != 0) {
		N += (sq_pr - (N % sq_pr));
		if(rank == 0) {
			cout << "Resizing original matrix from " << original_N << " to " << N << endl;
			cout << "original_N: " << original_N << endl;
			cout << "resized_N: " << N << endl;
		}
	}

    A_original = new Data*[N];    
    for(int i=0; i<N; i++) A_original[i] = new Data[N];


    Data** A = new Data*[N];
    for(int i=0; i<N; i++) {
        A[i] = new Data[N];
        memset(A[i], 0, sizeof(Data) * N);
    }

    Data** L = new Data*[N];
    for(int i=0; i<N; i++) {
        L[i] = new Data[N];
        memset(L[i], 0, sizeof(Data) * N);
    }

    Data** U = new Data*[N];
    for(int i=0; i<N; i++) {
        U[i] = new Data[N];
        memset(U[i], 0, sizeof(Data) * N);
    }

    for(int i=0; i<N; i++) L[i][i] = 1;

	// resized
	if(original_N != N) {
		for(int i=0; i<original_N; i++)
			for(int j=0; j<original_N; j++) 
				A[i][j] = (int)rand()%1000 + 1;
		for(int i=original_N; i<N; i++) {
			for(int j=original_N; j<N; j++) {
				if(i==j) A[i][j] = 1;
				else A[i][j] = 0;
			}
		}
	}

	// not resized
	else{
		for(int i=0; i<N; i++) {
			for(int j=0; j<N; j++) {
				A[i][j] = (int)rand() % 1000 + 1;
			}
		}
	}


	// if(rank == 0) {
	// 	Print_mat(A, N);
	// }
    Copy_mat(A, A_original, N);
	gettimeofday(&end_point, NULL);
	optime = (double)(end_point.tv_sec)+(double)(end_point.tv_usec)/1000000.0-(double)(start_point.tv_sec)-(double)(start_point.tv_usec)/1000000.0;
	if(rank == 0)
		cout << "Initialize matrix time: " << optime << " secs, now start LU" << endl;


    gettimeofday(&start_point, NULL);
    LU(A, L, U);

	Free_mat(A_original, N);
    Free_mat(A, N);
    Free_mat(L, N);
    Free_mat(U, N);
    MPI_Finalize();
    return 0;
}
