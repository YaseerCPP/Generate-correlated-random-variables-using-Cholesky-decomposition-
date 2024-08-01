#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

typedef vector<vector<double>> Matrix;

// Function to perform Cholesky Decomposition
void CholeskyDecomposition(const Matrix& Sigma, Matrix& L) {
    int n = Sigma.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k)
                sum += L[i][k] * L[j][k];

            if (i == j)
                L[i][j] = sqrt(Sigma[i][i] - sum);
            else
                L[i][j] = (Sigma[i][j] - sum) / L[j][j];
        }
    }
}

// Function to generate correlated normal random variables
vector<double> generateCorrelatedNormal(const Matrix& L, const vector<double>& Z) {
    int n = L.size();
    vector<double> Y(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j <= i; ++j) {
            sum += L[i][j] * Z[j];
        }
        Y[i] = sum;
    }

    return Y;
}

// Function to print a matrix
void printMatrix(const Matrix& M) {
    int n = M.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(10) << M[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    // Example: Correlation matrix Sigma
    Matrix Sigma = {{1.0, 0.5},
                    {0.5, 1.0}};

    int n = Sigma.size();
    Matrix L(n, vector<double>(n, 0.0));

    // Perform Cholesky decomposition on Sigma
    CholeskyDecomposition(Sigma, L);

    cout << "Cholesky Lower Triangular Matrix L:" << endl;
    printMatrix(L);

    // Generate independent standard normal variables Z
    vector<double> Z(n, 0.0);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < n; ++i) {
        Z[i] = dist(gen);
    }

    // Generate correlated normal variables Y = L * Z
    vector<double> Y = generateCorrelatedNormal(L, Z);

    cout << "Generated Correlated Normal Variables Y:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "Y[" << i << "] = " << Y[i] << endl;
    }

    return 0;
}
