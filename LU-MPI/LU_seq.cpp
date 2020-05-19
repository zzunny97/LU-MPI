#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

int N;
typedef double Data;
using namespace std;

void Print_mat(Data** p) {
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            //cout << p[i][j] << "\t";
			printf("%.3lf\t", p[i][j]);
        }
        cout << endl;
    }
    cout << endl;
}

void Copy_mat(Data** src, Data** dst) {
    for(int i=0; i<N; i++) 
        memcpy(dst[i], src[i], sizeof(Data)* N);
}

void Free_mat(Data** p) {
    for(int i=0; i<N; i++)
        delete[] p[i];
    delete[] p;
}


void verify(Data** A, Data** L, Data** U) {
    double diffMax = 0.0;
    Data** C = new Data*[N];
    for(int i=0; i<N; i++) {
        C[i] = new Data[N];
        memset(C[i], 0, sizeof(Data) * N);
    }

    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            for(int k=0; k<N; k++) 
                C[i][j] += L[i][k] * U[k][j];

    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            double target = A[i][j] - C[i][j];
            if(target > diffMax)
                diffMax = target;
        }
    }
    cout << "diffMax: " << diffMax << endl;
    Free_mat(C);
}

void LU(Data** A, Data** L, Data** U) {

    for(int k=0; k<N; k++) {
        U[k][k] = A[k][k];
        for(int i=k+1; i<N; i++) {
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }
        for(int i=k+1; i<N; i++) 
            for(int j=k+1; j<N; j++) 
                A[i][j] -= L[i][k] * U[k][j];
    }

}

int main(int argc, char** argv){
    if(argc != 3 ) {
        cout << "Usage ./proj3 N Seed" << endl;
        exit(-1);
    }
	struct timeval start, end;
	double time;

	gettimeofday(&start, NULL);
    N = atoi(argv[1]);
    int seed = atoi(argv[2]);
    srand(seed);

    cout << "Initialize" << endl;
    Data** A_original = new Data*[N];    
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
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            A[i][j] = rand()%10 + 1;

    Copy_mat(A, A_original);
    
    cout << "LU decompose" << endl;
    LU(A, L, U);
	gettimeofday(&end, NULL);
	time = (double)(end.tv_sec)+(double)(end.tv_usec)/1000000.0-(double)(start.tv_sec)-(double)(start.tv_usec)/1000000.0;

	cout << "Initialize and LU time: " << time << endl;

	gettimeofday(&start, NULL);
    cout << "Verify" << endl;
    verify(A_original, L, U);
	gettimeofday(&end, NULL);
	time = (double)(end.tv_sec)+(double)(end.tv_usec)/1000000.0-(double)(start.tv_sec)-(double)(start.tv_usec)/1000000.0;
	cout << "Verify time: " << time << endl;

    Free_mat(A);
    Free_mat(L);
    Free_mat(U);
    return 0;
}
