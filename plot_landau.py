import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import hsv_to_rgb # Importar la función de conversión HSV a RGB

# Parámetros de la simulación (ajustar si cambian en C++)
# Asegúrate de que Nx, Ny, num_states coincidan con la salida de tu C++
# Por ejemplo, si usaste LandauLevels simulation(20, 20, ...);
Nx = 70
Ny = 70
num_states = 20 # Si pediste 20 autovalores

def plot_complex_wavefunction_frame(state_idx):
    filename = f'wavefunction_state_{state_idx}.dat'
    try:
        data = np.loadtxt(filename, comments='#')
    except IOError:
        print(f"Error: No se encontró el archivo {filename}. Asegúrate de que el programa C++ se ejecutó correctamente.")
        return None, None, None, None # Retorna valores nulos para indicar un error

    # Las columnas son: x, y, Re(psi), Im(psi), |psi|^2, phase
    # Reconstruimos psi_complex a partir de las partes real e imaginaria guardadas
    real_psi = data[:, 2].reshape(Nx, Ny)
    imag_psi = data[:, 3].reshape(Nx, Ny)
    psi_complex_grid = real_psi + 1j * imag_psi

    # Obtener las coordenadas x e y para el extent del imshow
    x_coords_flat = data[:, 0]
    y_coords_flat = data[:, 1]
    
    unique_x = np.unique(x_coords_flat)
    unique_y = np.unique(y_coords_flat)
    min_x, max_x = unique_x.min(), unique_x.max()
    min_y, max_y = unique_y.min(), unique_y.max()

    # --- Lógica de codificación de fase-amplitud ---
    amplitude = np.abs(psi_complex_grid)
    
    # Normalizar la amplitud para el brillo (Value en HSV)
    max_amplitude = np.max(amplitude)
    if max_amplitude > 1e-10: # Evitar división por cero o valores muy pequeños
        normalized_amplitude = amplitude / max_amplitude
    else:
        normalized_amplitude = np.zeros_like(amplitude)

    # Calcular la fase y normalizarla para el tono (Hue en HSV)
    # np.angle devuelve el ángulo en radianes entre -pi y pi
    phase = np.angle(psi_complex_grid)
    hue = (phase + np.pi) / (2 * np.pi) # Normalizar de 0 a 1

    # Crear la imagen HSV (Hue, Saturation, Value)
    # La saturación (Saturation) la mantenemos en un valor alto (ej. 0.9 o 1.0) para colores vibrantes
    # El valor (Value) es la amplitud normalizada
    hsv_image = np.stack([hue, np.ones_like(hue) * 0.9, normalized_amplitude], axis=-1)
    
    # Convertir de HSV a RGB
    rgb_image = hsv_to_rgb(hsv_image)
    # -------------------------------------------------

    plt.clf() # Limpia la figura antes de dibujar el nuevo frame

    # Usamos imshow para mostrar la imagen RGB
    # Asegúrate de que el orden de las dimensiones de rgb_image sea (Ny, Nx, 3) si origin='lower' y extent es [xmin, xmax, ymin, ymax]
    # La reshape(Nx, Ny) y stack hacen que sea (Nx, Ny, 3). imshow espera (M, N, 3) donde M es filas (y) y N es columnas (x).
    # Por lo tanto, si tus datos x van con columnas e y con filas, puede ser que necesites transponer rgb_image antes de imshow
    # o ajustar la forma en la que se reconstruye psi_complex_grid.
    # Con el data[:,0].reshape(Nx,Ny) y data[:,1].reshape(Nx,Ny) y los bucles en C++, x varía en la primera dimensión (filas Nx)
    # y y en la segunda (columnas Ny). Por lo tanto, psi_complex_grid[i,j] corresponde a x[i], y[j].
    # Para imshow, esto es directamente mapeable a (row, col) si row=x, col=y.
    # PERO, imshow por defecto asume que la primera dimensión es el eje Y y la segunda el eje X.
    # Por eso, si nuestros datos están como (Nx, Ny) y queremos que Nx sea X e Ny sea Y, a menudo necesitamos transponer.
    # Para la coherencia con tu imagen, asumiendo que es un mapa 2D normal, `imshow` con `origin='lower'`
    # y datos `rgb_image` en forma `(Ny, Nx, 3)` sería ideal, o `(Nx, Ny, 3)` con una transposición.
    # La forma en que reconstruyes `psi_complex_grid` es `(Nx, Ny)`.
    # Si quieres que el eje X sea la primera dimensión de psi_complex_grid y el eje Y la segunda,
    # y imshow mapea la primera dimensión a Y y la segunda a X, entonces necesitas transponer `rgb_image`.

    # Intentemos con la transposición para asegurar que los ejes coinciden con x, y
    img = plt.imshow(rgb_image.transpose(1, 0, 2), origin='lower',
                     extent=[min_x, max_x, min_y, max_y], aspect='equal')
    
    plt.xlabel('x (Å)', fontsize=10)
    plt.ylabel('y (Å)', fontsize=10)
    plt.title(f'Nivel de Landau {state_idx} (Energía: {energies_eV_str_fmt(state_idx)})', fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=8)

    # Añadir texto explicativo para la codificación de color
    plt.text(max_x * 0.95, max_y * 0.95, 'Color: Fase\nBrillo: Amplitud',
             horizontalalignment='right', verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'),
             fontsize=8)


def load_energies():
    try:
        energies_eV = np.loadtxt('energies.txt')
        return energies_eV
    except IOError:
        print("Error: No se encontró el archivo energies.txt.")
        return np.zeros(num_states) # Retorna un array de ceros si el archivo no existe

energies_eV_loaded = load_energies()

def energies_eV_str_fmt(index):
    if index < len(energies_eV_loaded):
        return f"{energies_eV_loaded[index]:.6f} eV"
    return "N/A"


fig = plt.figure(figsize=(7, 6)) # Ajustar un poco el tamaño de la figura para que quepa todo

def update(frame):
    print(f'Renderizando frame {frame+1}/{num_states}')
    # plot_complex_wavefunction_frame ahora solo dibuja y no retorna artistas específicos para blit
    plot_complex_wavefunction_frame(frame)
    return [] # Con blit=False, solo necesitamos devolver una lista vacía

# CAMBIO CLAVE: blit=False (como antes, para evitar el error de 'unhashable type: list')
ani = FuncAnimation(fig, update, frames=num_states, blit=False, interval=500)

print(f"Guardando animación en landau_levels_animation_phase_amplitude.gif con {num_states} frames...")
ani.save('landau_levels_animation_phase_amplitude.gif', writer='pillow', fps=2)
print("Animación guardada.")

plt.show()
