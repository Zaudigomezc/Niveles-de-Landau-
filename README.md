# Niveles-de-Landau
Este proyecto simula los niveles de energía y las funciones de onda de un electrón en un campo magnético uniforme (efecto Landau) utilizando el método numérico de diferencias finitas y la biblioteca Eigen para la diagonalización de la matriz.

Descripción General
El proyecto consta de un programa principal en C++ que construye y diagonaliza el Hamiltoniano del sistema en una grilla 2D, y scripts en Python para la visualización de los resultados. La simulación calcula las energías propias y las funciones de onda correspondientes a los primeros estados de Landau.

Componentes principales:
LandauLevels.h / LandauLevels.cpp: Implementación de la clase LandauLevels que contiene la lógica para construir el Hamiltoniano cuántico del electrón en un campo magnético utilizando un gauge específico (Gauge de Landau: 
vecA=(0,Bx,0)). Resuelve el problema de autovalores para encontrar las energías y funciones de onda.

main.cpp: Punto de entrada del programa C++. Configura los parámetros de la simulación (tamaño de la grilla, dimensiones espaciales, campo magnético) y ejecuta el cálculo.

build.sh: Script de shell para compilar el código C++ y ejecutar la simulación. Está configurado para usar g++ y la biblioteca Eigen.

plot_landau.py: Script de Python que lee los datos de las funciones de onda y las energías generados por el programa C++. Crea una animación GIF que visualiza la función de onda compleja (codificando fase como color y amplitud como brillo) para diferentes niveles de Landau.

graficas.py: Un script de Python más completo para análisis y visualización de los resultados. Genera varias gráficas:

Espectro de energías comparado con la teoría analítica.

Representaciones 2D de la densidad de probabilidad, parte real, parte imaginaria y fase de las funciones de onda.

Visualización 3D de la densidad de probabilidad y la parte real.

Animación de la evolución de las funciones de onda a través de los estados.

Análisis teórico de cómo las energías y la longitud magnética cambian con el campo magnético.

Instalación y Requisitos
Requisitos de C++:
Compilador C++17: Se recomienda g++ (versión 7 o superior).

Eigen3: Biblioteca de álgebra lineal para C++. Asegúrate de tenerla instalada. En sistemas basados en Debian/Ubuntu, puedes instalarla con:

Bash

sudo apt-get install libeigen3-dev
Si Eigen está en una ubicación no estándar, deberás ajustar la ruta -I/usr/include/eigen3 en build.sh a la ruta correcta (por ejemplo, -I/opt/local/include/eigen3 o -I./eigen si la has descargado en el directorio del proyecto).

Requisitos de Python:
Python 3.x

Bibliotecas de Python:

numpy

matplotlib

Pillow (para guardar animaciones GIF)

Puedes instalar estas dependencias con pip:

Bash

pip install numpy matplotlib pillow
