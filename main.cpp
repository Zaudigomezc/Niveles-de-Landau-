#include "LandauLevels.h" // Incluye la declaración de tu clase
#include <iostream>       // Para cerr

int main() {
    try {
        // Puedes cambiar los parámetros aquí si lo deseas
        LandauLevels sim(40, 40, 200e-10, 200e-10, 5.0);
        sim.run();
    } catch(const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
