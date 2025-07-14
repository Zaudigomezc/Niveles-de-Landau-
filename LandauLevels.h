#ifndef LANDAU_LEVELS_H
#define LANDAU_LEVELS_H

#include <vector>
#include <complex>
#include <string>
#include <utility> // Para std::pair

// Incluir solo lo necesario de Eigen en el encabezado
// Dense para VectorXd, MatrixXd, Dynamic, Dynamic
// Sparse para SparseMatrix
// Eigenvalues para resolver autovalores (aunque la clase solver se usa en .cpp)
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Usar declaraciones 'using' de forma limitada o específica
// o calificarlas con Eigen:: y std::
using namespace Eigen;
using namespace std; // Se puede ser más específico (ej. std::vector, std::complex)

// Definiciones de tipos para mayor claridad
typedef SparseMatrix<complex<double>> SpMatrixComplex;
typedef Matrix<complex<double>, Dynamic, Dynamic> MatrixComplex;
typedef VectorXcd VectorComplex;

class LandauLevels {
private:
    int Nx, Ny;
    double Lx, Ly, dx, dy;
    VectorXd x, y;
    MatrixXd X_grid, Y_grid;

    // Constantes físicas (pueden ser constexpr o static const si es necesario)
    static constexpr double hbar = 1.0545718e-34; // Constante de Planck reducida
    static constexpr double m = 9.10938356e-31;    // Masa del electrón
    static constexpr double e = 1.60217662e-19;    // Carga elemental
    double B; // Campo magnético

    SpMatrixComplex H; // Hamiltoniano disperso
    MatrixXd A_y;     // Componente Y del potencial vectorial

    // Métodos privados auxiliares
    SparseMatrix<double> createLaplacian1D(int n, double h);
    SparseMatrix<double> createIdentity(int n);
    SpMatrixComplex kronecker(const SparseMatrix<double>& A, const SparseMatrix<double>& B);

public:
    // Constructor con valores por defecto
    LandauLevels(int nx = 100, int ny = 100, double lx = 200e-10, double ly = 200e-10, double b = 5.0);

    // Métodos públicos principales
    void buildHamiltonian();
    bool isHermitian(double tolerance = 1e-9);
    pair<VectorXd, MatrixComplex> solveEigenvalues(int k_states = 20);
    void saveWavefunction(const VectorComplex& psi, int state_idx, const string& filename);
    void run();
};

#endif // LANDAU_LEVELS_H
