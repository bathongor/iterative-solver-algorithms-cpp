#include <iostream>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
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

template <class T>	
double norm(T &a){
	double sumSqr = 0;
	for(int i = 0; i < a.size(); i++){
		sumSqr += a[i]*a[i];
	}
	return sqrt(sumSqr);
}

template <class T>	
void printVector(string name, T &x){
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
			if(isSymmetric == true && J[k]-1 != i){
				y[J[k]-1] += V[k]*x[i];
			}
		}
	}

	return y;
}

void CG(vector<int> &I, vector<int> &J, vector<double> &V, vector<double> &x0, vector<double> &b, int m, bool isSymmetric){
    vector<vector<double>> x(m+1, vector<double> (b.size(), 0)), r(m+1, vector<double> (b.size(), 0)), p(m+1, vector<double> (b.size(), 0));
    // r0 = b-A*x0
    // p0 = r0
    vector<double> Ax0 = matrixVectorProduct(I, J, V, x0, isSymmetric);
    for(int i = 0; i < b.size(); i++){
		r[0][i] = b[i] - Ax0[i];
        p[0][i] = r[0][i];
	}
	double normr0 = norm<vector<double>>(r[0]);
	int i = 0;
    while(i < m){
        vector<double> Apm = matrixVectorProduct(I, J, V, p[i], isSymmetric);
        double alpha = inner_product(r[i].begin(), r[i].end(), r[i].begin(), 0.0)/inner_product(Apm.begin(), Apm.end(), p[i].begin(), 0.0);
        // x[m+1] = xm + alpha*pm
        for(int j = 0; j < x[i].size(); j++){
            x[i+1][j] = x[i][j] + alpha*p[i][j];
            r[i+1][j] = r[i][j] - alpha*Apm[j];
        }
        double beta = inner_product(r[i+1].begin(), r[i+1].begin(), r[i+1].begin(), 0.0) / inner_product(r[i].begin(), r[i].end(), r[i].begin(), 0.0);
        // p[m+1] = r[m+1] + beta*pm
        for(int j = 0; j < x[i].size(); j++){
            p[i+1][j] = r[i+1][j] + beta*p[i][j];
        }
		double normri = norm<vector<double>>(r[i]);
		double convergence = normri/normr0;
		cout << "convergence: " << convergence << '\n';
		if(convergence < 0.00000001) {
			cout << "Converged!: took " << i << " iterations" << '\n';
			printVector<vector<double>>("xm",x[i]);
			return ;
		}
		i++;
    }
	cout << "Couldn't converge with given m!" << '\n';
	printVector<vector<double>>("xm", x[i]);
}

void computeMatrixVectorProduct(){

	// Open the file:
	ifstream fin("s3rmt3m3.mtx");
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
	bool isSymmetric = true;

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

	// printVector<vector<int>>("I", I);
	// printVector<vector<int>>("J", J);
	// printVector<vector<double>>("V", V);

	// Solution X vector
	vector<double> xstar(I.size()-1, 1);

	// initialize y with zerosf
	vector<double> ystar = matrixVectorProduct(I, J, V, xstar, isSymmetric);

	// Initial Guess X
	vector<double> x0(I.size()-1, 0);

    CG(I, J, V, x0, ystar, 5000, true);
	// cout << '\n';
}

int main(){
	computeMatrixVectorProduct();
	return 0;
}