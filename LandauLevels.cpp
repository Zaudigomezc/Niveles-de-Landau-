#include "LandauLevels.h" // Incluye la declaración de tu clase
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>

// Eigenvalues.h se incluye aquí porque ComplexEigenSolver se usa en el cuerpo de la función.
#include <Eigen/Eigenvalues>

using namespace std::chrono; // Para los cálculos de tiempo

// Implementación del constructor
LandauLevels::LandauLevels(int nx, int ny, double lx, double ly, double b)
    : Nx(nx), Ny(ny), Lx(lx), Ly(ly), B(b) {

    dx = Lx / Nx;
    dy = Ly / Ny;

    x = VectorXd::LinSpaced(Nx, -Lx/2, Lx/2);
    y = VectorXd::LinSpaced(Ny, -Ly/2, Ly/2);

    X_grid.resize(Nx, Ny);
    Y_grid.resize(Nx, Ny);
    A_y.resize(Nx, Ny);

    auto start_grid = high_resolution_clock::now();
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            X_grid(i, j) = x(i);
            Y_grid(i, j) = y(j);
            A_y(i, j) = B * x(i); // Potencial vectorial en Gauge de Landau (Ay = Bx)
        }
    }
    auto end_grid = high_resolution_clock::now();
    cout << "Tiempo creacion de mallas: "
         << duration_cast<milliseconds>(end_grid - start_grid).count()
         << " ms" << endl;
}

// Implementación de createLaplacian1D
SparseMatrix<double> LandauLevels::createLaplacian1D(int n, double h) {
    auto start = high_resolution_clock::now();
    SparseMatrix<double> T(n, n);
    vector<Triplet<double>> triplets;

    for(int i = 0; i < n; i++) {
        triplets.push_back(Triplet<double>(i, i, -2.0 / (h*h)));
        if(i > 0) triplets.push_back(Triplet<double>(i, i-1, 1.0 / (h*h)));
        if(i < n-1) triplets.push_back(Triplet<double>(i, i+1, 1.0 / (h*h)));
    }

    T.setFromTriplets(triplets.begin(), triplets.end());
    auto end = high_resolution_clock::now();
    cout << "  Tiempo Laplaciano 1D: "
         << duration_cast<microseconds>(end - start).count()
         << " μs" << endl;
    return T;
}

// Implementación de createIdentity
SparseMatrix<double> LandauLevels::createIdentity(int n) {
    auto start = high_resolution_clock::now();
    SparseMatrix<double> I(n, n);
    I.setIdentity();
    auto end = high_resolution_clock::now();
    cout << "  Tiempo matriz identidad: "
         << duration_cast<microseconds>(end - start).count()
         << " μs" << endl;
    return I;
}

// Implementación de kronecker
SpMatrixComplex LandauLevels::kronecker(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    auto start = high_resolution_clock::now();
    int m1 = A.rows(), n1 = A.cols();
    int m2 = B.rows(), n2 = B.cols();

    SpMatrixComplex result(m1 * m2, n1 * n2);
    vector<Triplet<complex<double>>> triplets;

    for(int k_A = 0; k_A < A.outerSize(); ++k_A) {
        for(SparseMatrix<double>::InnerIterator it_A(A, k_A); it_A; ++it_A) {
            int i1 = it_A.row(), j1 = it_A.col();
            double val1 = it_A.value();

            for(int k_B = 0; k_B < B.outerSize(); ++k_B) {
                for(SparseMatrix<double>::InnerIterator it_B(B, k_B); it_B; ++it_B) {
                    int i2 = it_B.row(), j2 = it_B.col();
                    double val2 = it_B.value();

                    triplets.push_back(Triplet<complex<double>>(i1 * m2 + i2, j1 * n2 + j2, val1 * val2));
                }
            }
        }
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    auto end = high_resolution_clock::now();
    cout << "  Tiempo producto Kronecker: "
         << duration_cast<milliseconds>(end - start).count()
         << " ms" << endl;
    return result;
}

// Implementación de buildHamiltonian
void LandauLevels::buildHamiltonian() {
    auto start_total = high_resolution_clock::now();
    cout << "Construyendo Hamiltoniano..." << endl;

    auto start_laplacian = high_resolution_clock::now();
    SparseMatrix<double> T_x = createLaplacian1D(Nx, dx);
    SparseMatrix<double> T_y = createLaplacian1D(Ny, dy);
    auto end_laplacian = high_resolution_clock::now();

    auto start_identity = high_resolution_clock::now();
    SparseMatrix<double> Ix = createIdentity(Nx);
    SparseMatrix<double> Iy = createIdentity(Ny);
    auto end_identity = high_resolution_clock::now();

    auto start_kronecker = high_resolution_clock::now();
    SpMatrixComplex Tx_full = kronecker(T_x, Iy);
    SpMatrixComplex Ty_full = kronecker(Ix, T_y);
    auto end_kronecker = high_resolution_clock::now();

    auto start_H0 = high_resolution_clock::now();
    SpMatrixComplex H0 = -(hbar*hbar/(2.0*m)) * (Tx_full + Ty_full);
    auto end_H0 = high_resolution_clock::now();

    // Términos con campo magnético
    // auto start_B_terms = high_resolution_clock::now(); // Esta línea fue eliminada/comentada para resolver la advertencia

    VectorXd A_y_flat(Nx * Ny);
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            A_y_flat(i * Ny + j) = A_y(i, j);
        }
    }

    // Término cuadrático: (e*e * Ay^2) / (2m)
    auto start_quad = high_resolution_clock::now();
    SpMatrixComplex H_B_quadratic(Nx * Ny, Nx * Ny);
    vector<Triplet<complex<double>>> quad_triplets;
    for(int i = 0; i < Nx * Ny; i++) {
        quad_triplets.push_back(Triplet<complex<double>>(i, i, (e*e)*A_y_flat(i)*A_y_flat(i)/(2.0*m)));
    }
    H_B_quadratic.setFromTriplets(quad_triplets.begin(), quad_triplets.end());
    auto end_quad = high_resolution_clock::now();

    // Término lineal: -(e/m) * Ay * Py
    auto start_linear = high_resolution_clock::now();
    SparseMatrix<double> D_y_1D(Ny, Ny); // Operador de derivada parcial en y (aproximación de diferencia central)
    vector<Triplet<double>> dy_triplets;
    for(int i = 0; i < Ny; i++) {
        if(i > 0) dy_triplets.push_back(Triplet<double>(i, i-1, -1.0/(2.0*dy)));
        if(i < Ny-1) dy_triplets.push_back(Triplet<double>(i, i+1, 1.0/(2.0*dy)));
    }
    D_y_1D.setFromTriplets(dy_triplets.begin(), dy_triplets.end());

    // P_y = -i * hbar * d/dy. Aquí representamos d/dy y luego multiplicamos por -i*hbar
    SpMatrixComplex P_y_operator = kronecker(Ix, D_y_1D) * complex<double>(0, -hbar);

    // Matriz diagonal de Ay_flat
    SpMatrixComplex Ay_matrix(Nx * Ny, Nx * Ny);
    vector<Triplet<complex<double>>> ay_triplets;
    for(int i = 0; i < Nx * Ny; i++) {
        ay_triplets.push_back(Triplet<complex<double>>(i, i, A_y_flat(i)));
    }
    Ay_matrix.setFromTriplets(ay_triplets.begin(), ay_triplets.end());

    // El término lineal completo: -(e/m) * Ay_matrix * P_y_operator
    SpMatrixComplex H_B_linear = -(e/m) * (Ay_matrix * P_y_operator);
    auto end_linear = high_resolution_clock::now();

    // Suma de todos los términos del Hamiltoniano
    H = H0 + H_B_quadratic + H_B_linear;
    auto end_total = high_resolution_clock::now();

    cout << "Tiempos parciales:\n"
         << "  Laplacianos: " << duration_cast<milliseconds>(end_laplacian - start_laplacian).count() << " ms\n"
         << "  Identidades: " << duration_cast<milliseconds>(end_identity - start_identity).count() << " ms\n"
         << "  Kronecker: " << duration_cast<milliseconds>(end_kronecker - start_kronecker).count() << " ms\n"
         << "  H0: " << duration_cast<milliseconds>(end_H0 - start_H0).count() << " ms\n"
         << "  Termino cuadratico: " << duration_cast<milliseconds>(end_quad - start_quad).count() << " ms\n"
         << "  Termino lineal: " << duration_cast<milliseconds>(end_linear - start_linear).count() << " ms\n"
         << "Tiempo total construccion Hamiltoniano: "
         << duration_cast<milliseconds>(end_total - start_total).count()
         << " ms" << endl;
}

// Implementación de isHermitian
bool LandauLevels::isHermitian(double tolerance) {
    auto start = high_resolution_clock::now();
    SpMatrixComplex H_conj_T = H.adjoint(); // H.adjoint() calcula la traspuesta conjugada

    // Comprobación de la diferencia entre H y H_conj_T
    // Convertir a densa para una resta más sencilla y comparación elemento a elemento
    MatrixComplex diff = MatrixComplex(H) - MatrixComplex(H_conj_T);

    double max_diff = 0.0;
    for (int i = 0; i < diff.rows(); ++i) {
        for (int j = 0; j < diff.cols(); ++j) {
            max_diff = max(max_diff, abs(diff(i, j)));
        }
    }

    auto end = high_resolution_clock::now();
    cout << "Tiempo verificacion hermiticidad: "
         << duration_cast<milliseconds>(end - start).count()
         << " ms" << endl;

    return max_diff < tolerance;
}

// Implementación de solveEigenvalues
pair<VectorXd, MatrixComplex> LandauLevels::solveEigenvalues(int k_states) {
    auto start_total = high_resolution_clock::now();

    cout << "\n--- Verificacion de Hermiticidad ---" << endl;
    bool hermitian = isHermitian();
    cout << "Es el Hamiltoniano hermitiano? " << (hermitian ? "Si" : "No") << endl;

    int matrix_size = H.rows();
    cout << "\n--- Resolviendo eigenvalores ---" << endl;
    cout << "Tamano de matriz: " << matrix_size << "x" << matrix_size << endl;

    auto start_convert = high_resolution_clock::now();
    MatrixComplex H_dense = MatrixComplex(H); // Convierte la matriz dispersa a densa
    auto end_convert = high_resolution_clock::now();

    auto start_solve = high_resolution_clock::now();
    // ComplexEigenSolver calcula los autovalores y autovectores de una matriz compleja densa
    ComplexEigenSolver<MatrixComplex> solver(H_dense);
    auto end_solve = high_resolution_clock::now();

    if(solver.info() != Success) {
        cerr << "Error en el calculo de eigenvalores" << endl;
        return make_pair(VectorXd(), MatrixComplex());
    }

    VectorXd eigenvals_real = solver.eigenvalues().real(); // Los autovalores deben ser reales para un Hamiltoniano Hermitiano
    MatrixComplex eigenvecs = solver.eigenvectors();

    // Ordenar los autovalores y seleccionar los k_states más bajos
    vector<pair<double, int>> indexed_eigenvals;
    for(int i = 0; i < eigenvals_real.size(); i++) {
        indexed_eigenvals.push_back(make_pair(eigenvals_real(i), i));
    }
    sort(indexed_eigenvals.begin(), indexed_eigenvals.end());

    int num_states = min(k_states, (int)eigenvals_real.size());
    VectorXd sorted_eigenvals(num_states);
    MatrixComplex sorted_eigenvecs(eigenvecs.rows(), num_states);

    for(int i = 0; i < num_states; i++) {
        sorted_eigenvals(i) = indexed_eigenvals[i].first;
        sorted_eigenvecs.col(i) = eigenvecs.col(indexed_eigenvals[i].second);
    }

    auto end_total = high_resolution_clock::now();

    cout << "Tiempos de solucion:\n"
         << "  Conversion a densa: " << duration_cast<milliseconds>(end_convert - start_convert).count() << " ms\n"
         << "  Calculo eigenvalores: " << duration_cast<milliseconds>(end_solve - start_solve).count() << " ms\n"
         << "Tiempo total solveEigenvalues: "
         << duration_cast<milliseconds>(end_total - start_total).count()
         << " ms" << endl;

    return make_pair(sorted_eigenvals, sorted_eigenvecs);
}

// Implementación de saveWavefunction
void LandauLevels::saveWavefunction(const VectorComplex& psi, int state_idx, const string& filename) {
    auto start = high_resolution_clock::now();
    ofstream file(filename);
    if(!file.is_open()) {
        cerr << "Error abriendo archivo " << filename << endl;
        return;
    }

    file << "# x(A) y(A) Re(psi) Im(psi) |psi|^2 phase(rad)\n";
    file << scientific << setprecision(10); // Formato científico y precisión

    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            int idx = i * Ny + j; // Mapeo 2D a 1D
            complex<double> val = psi(idx);
            file << x(i)*1e10 << " " << y(j)*1e10 << " " // Convertir a Angstroms
                 << val.real() << " " << val.imag() << " "
                 << norm(val) << " " << arg(val) << "\n"; // norm(val) es |psi|^2, arg(val) es la fase
        }
        file << "\n"; // Línea vacía para separar "bloques" para Gnuplot o Matplotlib 3D
    }

    file.close();
    auto end = high_resolution_clock::now();
    cout << "  Tiempo guardado estado " << state_idx << ": "
         << duration_cast<milliseconds>(end - start).count()
         << " ms" << endl;
}

// Implementación de run
void LandauLevels::run() {
    auto start_total = high_resolution_clock::now();
    cout << "=== Simulacion de Niveles de Landau ===" << endl;
    cout << "Parametros:\n"
         << "  Nx = " << Nx << ", Ny = " << Ny << "\n"
         << "  Lx = " << Lx*1e10 << " A, Ly = " << Ly*1e10 << " A\n" // Muestra en Angstroms
         << "  B = " << B << " T" << endl;

    buildHamiltonian(); // Construye el operador Hamiltoniano
    auto [energies, states] = solveEigenvalues(20); // Resuelve los autovalores y autovectores

    if(energies.size() == 0) {
        cerr << "Error: No se calcularon eigenvalores" << endl;
        return;
    }

    VectorXd energies_eV = energies / e; // Convierte energías a eV
    cout << "\nEnergias (eV):\n";
    for(int i = 0; i < energies_eV.size(); i++) {
        cout << "  Estado " << i << ": " << fixed << setprecision(6) << energies_eV(i) << " eV\n";
    }

    ofstream energy_file("energies.txt");
    if(energy_file.is_open()) {
        energy_file << scientific << setprecision(10);
        for(int i = 0; i < energies_eV.size(); i++) {
            energy_file << energies_eV(i) << "\n";
        }
        energy_file.close();
        cout << "Energias guardadas en energies.txt" << endl;
    }

    cout << "\nGuardando funciones de onda..." << endl;
    for(int i = 0; i < energies.size(); i++) {
        string filename = "wavefunction_state_" + to_string(i) + ".dat";
        saveWavefunction(states.col(i), i, filename);
    }

    auto end_total = high_resolution_clock::now();
    cout << "\n=== Simulacion completada ===" << endl;
    cout << "Tiempo total de ejecucion: "
         << duration_cast<milliseconds>(end_total - start_total).count()
         << " ms" << endl;

    // Script Python para visualización (se mantiene aquí por conveniencia)
    cout << "\nScript Python para visualizacion:\n"
         << "```python\n"
         << "import numpy as np\n"
         << "import matplotlib.pyplot as plt\n"
         << "from matplotlib.animation import FuncAnimation\n\n"
         << "# Parametros (deben coincidir con la simulacion C++)\n"
         << "Nx, Ny = " << Nx << ", " << Ny << "\n"
         << "num_states = " << energies.size() << "\n\n"
         << "def load_wavefunction(state_idx):\n"
         << "    filename = f'wavefunction_state_{state_idx}.dat'\n"
         << "    data = np.loadtxt(filename, comments='#')\n"
         << "    x = data[:, 0].reshape(Nx, Ny)\n"
         << "    y = data[:, 1].reshape(Nx, Ny)\n"
         << "    psi_sq = data[:, 4].reshape(Nx, Ny)  # |ψ|²\n"
         << "    return x, y, psi_sq\n\n"
         << "def plot_state(state_idx):\n"
         << "    x, y, psi_sq = load_wavefunction(state_idx)\n"
         << "    plt.clf()\n"
         << "    plt.imshow(psi_sq.T, extent=[x.min(), x.max(), y.min(), y.max()],\n"
         << "               origin='lower', aspect='auto', cmap='inferno')\n"
         << "    plt.colorbar(label='|ψ|²')\n"
         << "    plt.title(f'Nivel de Landau {state_idx}')\n"
         << "    plt.xlabel('x (Å)')\n"
         << "    plt.ylabel('y (Å)')\n\n"
         << "# Animacion\n"
         << "fig = plt.figure(figsize=(8, 6))\n"
         << "energies = np.loadtxt('energies.txt')\n\n"
         << "def update(frame):\n"
         << "    plot_state(frame)\n"
         << "    plt.title(f'Nivel {frame} - Energia: {energies[frame]:.4f} eV')\n"
         << "    return fig,\n\n"
         << "ani = FuncAnimation(fig, update, frames=num_states, interval=500, blit=True)\n"
         << "plt.show()\n"
         << "# ani.save('landau_levels.gif', writer='pillow', fps=2)  # Descomentar para guardar\n"
         << "```" << endl;
}
