#!/bin/bash

# Compilador C++
CXX=g++

# Banderas del compilador
# -std=c++17: Utiliza el estándar C++17
# -Wall: Habilita todas las advertencias comunes
# -O2: Nivel de optimización 2 (buen equilibrio entre velocidad y tamaño)
# -I/usr/include/eigen3: Ruta a los encabezados de Eigen.
#     AJUSTA ESTO SI TUS ENCABEZADOS DE EIGEN ESTÁN EN UNA RUTA DIFERENTE.
#     Ejemplos:
#       - Si están en /opt/local/include/eigen3: -I/opt/local/include/eigen3
#       - Si descargaste Eigen y está en una carpeta 'eigen' en el mismo directorio que tu proyecto: -I./eigen
# -march=native: Optimiza el código para la arquitectura de tu CPU (puede mejorar el rendimiento)
CXXFLAGS="-std=c++17 -Wall -O2 -I/usr/include/eigen3 -march=native"

# Archivos fuente (todos los .cpp que componen tu programa)
SOURCES="LandauLevels.cpp main.cpp"

# Nombre del ejecutable de salida
EXECUTABLE="landau_sim"

echo "Iniciando compilación..."
# Compilar los archivos fuente y enlazarlos para crear el ejecutable
$CXX $CXXFLAGS $SOURCES -o $EXECUTABLE

# Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Compilación exitosa. Ejecutando $EXECUTABLE..."
    ./$EXECUTABLE
else
    echo "¡Error de compilación! Por favor, revisa los mensajes anteriores."
fi
