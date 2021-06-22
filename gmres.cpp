#include <iostream>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

struct Pair {
	int rowIndex;
	int colIndex;
	double value;

	void setPair(int row, int col, double val){
		this->rowIndex = row;
		this->colIndex = col;
		this->value = val;
	}
};

double norm(vector<double> &a){
	double sumSqr = 0;
	for(int i = 0; i < a.size(); i++){
		sumSqr += a[i]*a[i];
	}
	return sqrt(sumSqr);
}

void printVector(string name, vector<double> &x){
	cout << name << ": ";
	for(int i = 0; i < x.size(); i++){	
		cout << " " << x[i];
	}
	cout << '\n';
}

void printMatrix(string name, vector<vector<double>> &A){
	cout << name << ": " << '\n';
	for(int i = 0; i < A.size(); i++){
		for(int j = 0; j < A[i].size(); j++){	
			cout << " " << A[i][j];
		}
		cout << '\n';
	}
	cout << '\n';
}

bool sortHelper(Pair a, Pair b){
	return a.rowIndex < b.rowIndex;
}

bool checkIsSymmetric(Pair *csrMatrix, Pair *originalMatrix, int L){
	for(int i=0; i < L; i++){
		if(csrMatrix[i].colIndex != originalMatrix[i].rowIndex || csrMatrix[i].value != originalMatrix[i].value){
			return false; 
		}
	}
	return true;
}

vector<double> matrixVectorProduct(vector<int> &I, vector<int> &J, vector<double> &V, vector<double> &x, bool isSymmetric){
	vector<double> y(x.size(), 0);
	for(int i = 0; i < I.size()-1; i++){
		int i1 = I[i]-1;
		int i2 = I[i+1]-1;
		for(int k = i1; k < i2; k++){
			y[i] += V[k]*x[J[k]-1];
			if(isSymmetric && J[k] != i+1){
				y[J[k]-1] += V[k]*x[i];
			}
		}
	}

	return y;
}

// void backwardsSubs(){
//     int m = 3;
//     double R[3][3] = {{1, 1, 1}, {0, 2, 2}, {0, 0, 3}};
//     double g[3] = {4, 5, 6};
//     vector<double> y1(m, 0);
//     for(int i = m-1; i >= 0; i--){
//         double sum = 0;
//         for(int j = m-1; j >= i+1; j--){
//             sum += R[i][j] * y1[j];
//         }
//         y1[i] = (g[i]-sum)/R[i][i];
//     }
//     printVector("y", y1);
// }

// void calculateXm(){
//     int m = 3;
//     double v[3][3] = {{1, 1, 1}, {0, 2, 2}, {0, 0, 3}};
//     double y[3] = {2, 0, 1};
//     vector<double> y1(m, 0);

//     vector<double> sum(m, 0);
//     for(int k = 0; k < m; k++){
//         for(int i = 0; i < m; i++){
//             sum[i] += y[k]*v[k][i];
//         }
//     }
//     printVector("sum", sum);
// }


void getKrylov(vector<int> &I, vector<int> &J, vector<double> &V, vector<vector<double>> &v, int j, vector<vector<double>> &h, int m, bool isSymmetric){
	vector<double> w = matrixVectorProduct(I, J, V, v[j], isSymmetric);
	vector<double> h_colj;
	for(int i = 0; i <= j; i++){
		h_colj.push_back(inner_product(w.begin(), w.end(), v[i].begin(), 0.0));
		for(int k = 0; k < w.size(); k++){
			w[k] = w[k] - h_colj[i]* v[i][k];
		}
	}
	h_colj.push_back(norm(w));
	vector<double> v_next;
	for(int l = 0; l < w.size(); l++){
		v_next.push_back(w[l]/h_colj[j+1]);
	}
	h_colj.resize(m+1, 0);
	h[j] = h_colj;
	v.push_back(v_next);
}

void GMRES(vector<int> &I, vector<int> &J, vector<double> &V, vector<double> &x0, vector<double> &b, int m, bool isSymmetric, double &rho){
    vector<double> r0(x0.size(), 0), g(m+1, 0);
	g[0] = 1;
	vector<vector<double>> v, h(m, vector<double> (m+1, 0));
    double R[m+2][m+1] = {0};
	double normr0 = 0;

	// r0 = b-A*x0
	vector<double> Ax = matrixVectorProduct(I, J, V, x0, isSymmetric);
	for(int i = 0; i < r0.size(); i++){
		r0[i] = b[i] - Ax[i];
		normr0 += r0[i]*r0[i];
	}
	// v1 = r0/||r0||
	normr0 = sqrt(normr0);
	for(int i = 0; i < r0.size(); i++){
		r0[i] = r0[i]/normr0;
	}
	v.push_back(r0);
	// g = |r0|e1
	for(int i = 0; i < g.size(); i++){
		g[i] = normr0*g[i];
	}

	vector<double> c(m+5, 0), s(m+5, 0);
    int iter = -1;
	for(int j = 0; j < m; j++){
		getKrylov(I, J, V, v, j, h, m, isSymmetric);
		for(int k = 1; k <= j; k++){
			double tempPrev = c[k-1]*h[j][k-1]+s[k-1]*h[j][k];
			h[j][k] = -s[k-1]*h[j][k-1]+c[k-1]*h[j][k];
			h[j][k-1] = tempPrev;
		}
		double denom = sqrt(h[j][j]*h[j][j]+h[j][j+1]*h[j][j+1]);
		c[j] = (double)h[j][j]/denom;
		s[j] = (double)h[j][j+1]/denom;
		h[j][j] = denom;
		g[j+1] = -s[j]*g[j];
        iter++;
        rho = abs(g[j+1]); // updated 
		g[j] = c[j]*g[j];
	}
	// printMatrix("H(m)", h);
    // Transform current H[m] to R[m] matrix
    for(int i = 0; i < h.size(); i++){
        for(int j = 0; j < h[i].size(); j++){
            // if(j == i+1) continue;
            R[j][i] = h[i][j];
        }
    }

    // for(int i = 0; i < h.size(); i++){
    //     for(int j = 0; j < h[i].size(); j++){
    //         cout << R[i][j] << ' ';
    //     }
    //     cout << '\n';   
    // }


    // // Backwards substition Ry = g
    vector<double> y(m, 0);
    for(int i = m-1; i >= 0; i--){
        double sum = 0;
        for(int j = m-1; j >= i+1; j--){
            sum += R[i][j] * y[j];
        }
        y[i] = (g[i]-sum)/R[i][i];
    }

    // // printVector("y", y);

    // // xm = x0 + yk*vk
    vector<double> xm(x0.size(), 0);
    vector<double> sum(x0.size(), 0);
    for(int k = 0; k < m; k++){
        for(int i = 0; i < v[k].size(); i++){
            sum[i] += y[k]*v[k][i];
        }
    }
    // printVector("sum", sum);
    for(int i = 0; i < x0.size(); i++){
        xm[i] = x0[i]+sum[i];
    }

    printVector("xm", xm);
    // x0 = xm;
}

void computeMatrixVectorProduct(){

	// Open the file:
	ifstream fin("orsirr_1.mtx");
	int M, N, L;
	while(fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> M >> N >> L;
	Pair csrMatrix[L+1], originalMatrix[L+1];

	for(int i = 0; i < L; i++){
		int m, n;
		double data;
		fin >> m >> n >> data;
		csrMatrix[i].setPair(m, n, data);
		originalMatrix[i].setPair(m, n, data);
	}

	fin.close();
	sort(csrMatrix, csrMatrix+L, sortHelper);
	bool isSymmetric = checkIsSymmetric(csrMatrix, originalMatrix, L);

	// Initialize your matrix;
	vector<int> I, J; // initialize row and col index array for CSR
	vector<double> V; // initialize values array that contain CSR matrix values

	// Assemble CSR row, col and V based on the symmetry
	int prevRow = -1;
	for(int i = 0; i < L; i++){
		if((isSymmetric == true && csrMatrix[i].colIndex <= csrMatrix[i].rowIndex) || isSymmetric == false){
			J.push_back(csrMatrix[i].colIndex);
			V.push_back(csrMatrix[i].value);
			if(prevRow != csrMatrix[i].rowIndex){
				I.push_back(J.size());
			}
			prevRow = csrMatrix[i].rowIndex;	
		}
	} 
	I.push_back(J.size()+1);

	// Solution X vector
	vector<double> xstar(I.size()-1, 1);

	// initialize y with zerosf
	vector<double> ystar;

	ystar = matrixVectorProduct(I, J, V, xstar, isSymmetric);

	// Initial Guess X
	vector<double> x(I.size()-1, 0);
    double rho = 1;
	int countIter = 0;

	// Restarted GMRES
    // while(rho > 0.0000001){
    //     GMRES(I, J, V, x, ystar, 40, isSymmetric, rho);
	// 	countIter++;
    // }

	// Full GMRES
	GMRES(I, J, V, x, ystar, 500, isSymmetric, rho);

    cout << "countIter: " << countIter*40;
    // printVector("xstar", xstar);
	// printVector("ystar", ystar);
	cout << '\n';
}

int main(){
	computeMatrixVectorProduct();
	return 0;
}