import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
from matplotlib.patches import Circle
import os
import glob

# Configuración de matplotlib para mejor calidad
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10

class LandauVisualization:
    def __init__(self, Nx=40, Ny=40, B=5.0):
        """
        Inicializar el visualizador de niveles de Landau
        
        Parámetros:
        Nx, Ny: Dimensiones de la grilla (deben coincidir con la simulación C++)
        B: Campo magnético en Tesla
        """
        self.Nx = Nx
        self.Ny = Ny
        self.B = B
        
        # Constantes físicas
        self.hbar = 1.0545718e-34  # J·s
        self.m = 9.10938356e-31    # kg (masa del electrón)
        self.e = 1.60217662e-19    # C (carga elemental)
        
        # Longitud magnética
        self.l_B = np.sqrt(self.hbar / (self.e * self.B))
        
        # Cargar datos
        self.load_data()
        
    def load_data(self):
        """Cargar energías y verificar archivos de funciones de onda"""
        try:
            self.energies_eV = np.loadtxt('energies.txt')
            self.num_states = len(self.energies_eV)
            print(f"Cargadas {self.num_states} energías")
        except FileNotFoundError:
            print("Error: No se encontró 'energies.txt'. Ejecuta primero la simulación C++.")
            self.energies_eV = np.array([])
            self.num_states = 0
            
        # Verificar archivos de funciones de onda
        self.wavefunction_files = glob.glob('wavefunction_state_*.dat')
        self.wavefunction_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
        print(f"Encontrados {len(self.wavefunction_files)} archivos de funciones de onda")
        
    def load_wavefunction(self, state_idx):
        """Cargar una función de onda específica"""
        filename = f'wavefunction_state_{state_idx}.dat'
        if not os.path.exists(filename):
            print(f"Advertencia: No se encontró {filename}")
            return None, None, None, None, None, None
            
        data = np.loadtxt(filename, comments='#')
        
        # Las columnas son: x, y, Re(psi), Im(psi), |psi|^2, phase
        x_coords = data[:, 0].reshape(self.Nx, self.Ny)
        y_coords = data[:, 1].reshape(self.Nx, self.Ny)
        real_part = data[:, 2].reshape(self.Nx, self.Ny)
        imag_part = data[:, 3].reshape(self.Nx, self.Ny)
        amplitude_sq = data[:, 4].reshape(self.Nx, self.Ny)
        phase = data[:, 5].reshape(self.Nx, self.Ny)
        
        return x_coords, y_coords, real_part, imag_part, amplitude_sq, phase
    
    def plot_energy_spectrum(self):
        """Graficar el espectro de energías vs número cuántico"""
        if len(self.energies_eV) == 0:
            print("No hay datos de energías para graficar")
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Espectro de energías
        ax1.plot(range(len(self.energies_eV)), self.energies_eV, 'bo-', markersize=6)
        ax1.set_xlabel('Número de estado')
        ax1.set_ylabel('Energía (eV)')
        ax1.set_title(f'Espectro de Niveles de Landau\nB = {self.B} T')
        ax1.grid(True, alpha=0.3)
        
        # Comparación con teoría analítica
        # Para niveles de Landau: E_n = ℏωc(n + 1/2)
        omega_c = self.e * self.B / self.m  # Frecuencia ciclotrón
        n_values = np.arange(len(self.energies_eV))
        E_theory = self.hbar * omega_c * (n_values + 0.5) / self.e  # en eV
        
        ax2.plot(n_values, self.energies_eV, 'bo-', label='Simulación', markersize=6)
        ax2.plot(n_values, E_theory, 'r--', label='Teoría: ℏωc(n+1/2)', linewidth=2)
        ax2.set_xlabel('Número cuántico n')
        ax2.set_ylabel('Energía (eV)')
        ax2.set_title('Comparación con Teoría Analítica')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('energy_spectrum.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Imprimir información
        print(f"\n=== Análisis del Espectro ===")
        print(f"Longitud magnética: {self.l_B*1e10:.2f} Å")
        print(f"Frecuencia ciclotrón: {omega_c/2/np.pi/1e12:.2f} THz")
        print(f"Energía cuántica ℏωc: {self.hbar*omega_c/self.e*1000:.2f} meV")
        
    def plot_wavefunction_2d(self, state_idx):
        """Graficar función de onda 2D con múltiples representaciones"""
        x, y, real_part, imag_part, amplitude_sq, phase = self.load_wavefunction(state_idx)
        
        if x is None:
            return
            
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Configurar extensión para imshow
        extent = [x.min(), x.max(), y.min(), y.max()]
        
        # 1. Densidad de probabilidad |ψ|²
        im1 = axes[0,0].imshow(amplitude_sq.T, extent=extent, origin='lower', 
                              aspect='auto', cmap='hot')
        axes[0,0].set_title(f'|ψ|² - Estado {state_idx}')
        axes[0,0].set_xlabel('x (Å)')
        axes[0,0].set_ylabel('y (Å)')
        plt.colorbar(im1, ax=axes[0,0])
        
        # 2. Parte real
        im2 = axes[0,1].imshow(real_part.T, extent=extent, origin='lower', 
                              aspect='auto', cmap='RdBu_r')
        axes[0,1].set_title('Re(ψ)')
        axes[0,1].set_xlabel('x (Å)')
        axes[0,1].set_ylabel('y (Å)')
        plt.colorbar(im2, ax=axes[0,1])
        
        # 3. Parte imaginaria
        im3 = axes[0,2].imshow(imag_part.T, extent=extent, origin='lower', 
                              aspect='auto', cmap='RdBu_r')
        axes[0,2].set_title('Im(ψ)')
        axes[0,2].set_xlabel('x (Å)')
        axes[0,2].set_ylabel('y (Å)')
        plt.colorbar(im3, ax=axes[0,2])
        
        # 4. Fase
        im4 = axes[1,0].imshow(phase.T, extent=extent, origin='lower', 
                              aspect='auto', cmap='hsv')
        axes[1,0].set_title('Fase(ψ)')
        axes[1,0].set_xlabel('x (Å)')
        axes[1,0].set_ylabel('y (Å)')
        plt.colorbar(im4, ax=axes[1,0])
        
        # 5. Contornos de |ψ|²
        axes[1,1].contour(x, y, amplitude_sq, levels=10, colors='black', alpha=0.7)
        axes[1,1].contourf(x, y, amplitude_sq, levels=20, cmap='viridis', alpha=0.8)
        axes[1,1].set_title('Contornos de |ψ|²')
        axes[1,1].set_xlabel('x (Å)')
        axes[1,1].set_ylabel('y (Å)')
        axes[1,1].set_aspect('equal')
        
        # 6. Corte transversal en el centro
        center_y_idx = self.Ny // 2
        center_x_idx = self.Nx // 2
        
        axes[1,2].plot(x[:, center_y_idx], amplitude_sq[:, center_y_idx], 'b-', 
                      label=f'Corte en y={y[center_x_idx, center_y_idx]:.1f} Å')
        axes[1,2].plot(y[center_x_idx, :], amplitude_sq[center_x_idx, :], 'r-', 
                      label=f'Corte en x={x[center_x_idx, center_y_idx]:.1f} Å')
        axes[1,2].set_xlabel('Posición (Å)')
        axes[1,2].set_ylabel('|ψ|²')
        axes[1,2].set_title('Cortes Transversales')
        axes[1,2].legend()
        axes[1,2].grid(True, alpha=0.3)
        
        plt.suptitle(f'Estado {state_idx} - Energía: {self.energies_eV[state_idx]:.6f} eV', 
                    fontsize=14)
        plt.tight_layout()
        plt.savefig(f'wavefunction_2d_state_{state_idx}.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_wavefunction_3d(self, state_idx):
        """Graficar función de onda en 3D"""
        x, y, real_part, imag_part, amplitude_sq, phase = self.load_wavefunction(state_idx)
        
        if x is None:
            return
            
        fig = plt.figure(figsize=(12, 5))
        
        # Superficie 3D de |ψ|²
        ax1 = fig.add_subplot(121, projection='3d')
        surf1 = ax1.plot_surface(x, y, amplitude_sq, cmap='hot', alpha=0.8)
        ax1.set_xlabel('x (Å)')
        ax1.set_ylabel('y (Å)')
        ax1.set_zlabel('|ψ|²')
        ax1.set_title(f'|ψ|² - Estado {state_idx}')
        
        # Superficie 3D de la parte real
        ax2 = fig.add_subplot(122, projection='3d')
        surf2 = ax2.plot_surface(x, y, real_part, cmap='RdBu_r', alpha=0.8)
        ax2.set_xlabel('x (Å)')
        ax2.set_ylabel('y (Å)')
        ax2.set_zlabel('Re(ψ)')
        ax2.set_title(f'Re(ψ) - Estado {state_idx}')
        
        plt.suptitle(f'Visualización 3D - Energía: {self.energies_eV[state_idx]:.6f} eV')
        plt.tight_layout()
        plt.savefig(f'wavefunction_3d_state_{state_idx}.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def create_animation(self):
        """Crear animación de todos los estados"""
        if len(self.wavefunction_files) == 0:
            print("No hay archivos de funciones de onda para animar")
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        def update(frame):
            ax1.clear()
            ax2.clear()
            
            x, y, real_part, imag_part, amplitude_sq, phase = self.load_wavefunction(frame)
            
            if x is None:
                return
                
            extent = [x.min(), x.max(), y.min(), y.max()]
            
            # Densidad de probabilidad
            im1 = ax1.imshow(amplitude_sq.T, extent=extent, origin='lower', 
                           aspect='auto', cmap='hot', vmin=0, vmax=np.max(amplitude_sq))
            ax1.set_title(f'|ψ|² - Estado {frame}')
            ax1.set_xlabel('x (Å)')
            ax1.set_ylabel('y (Å)')
            
            # Parte real
            max_real = np.max(np.abs(real_part))
            im2 = ax2.imshow(real_part.T, extent=extent, origin='lower', 
                           aspect='auto', cmap='RdBu_r', vmin=-max_real, vmax=max_real)
            ax2.set_title(f'Re(ψ) - Estado {frame}')
            ax2.set_xlabel('x (Å)')
            ax2.set_ylabel('y (Å)')
            
            # Añadir información de energía
            energy_text = f'Energía: {self.energies_eV[frame]:.6f} eV'
            fig.suptitle(energy_text, fontsize=12)
            
            print(f'Renderizando frame {frame+1}/{len(self.wavefunction_files)}')
            
        ani = FuncAnimation(fig, update, frames=len(self.wavefunction_files), 
                          interval=800, repeat=True)
        
        # Guardar animación
        print("Guardando animación...")
        ani.save('landau_levels_animation.gif', writer='pillow', fps=1.25)
        print("Animación guardada como 'landau_levels_animation.gif'")
        plt.show()
        
    def magnetic_field_comparison(self, B_values=[1, 3, 5, 7, 10]):
        """
        Comparar cómo cambian las energías con diferentes campos magnéticos
        (Requiere ejecutar la simulación C++ para cada valor de B)
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Gráfica teórica de cómo varían las energías con B
        n_max = 10
        for n in range(n_max):
            energies_theory = []
            for B in B_values:
                omega_c = self.e * B / self.m
                E_n = self.hbar * omega_c * (n + 0.5) / self.e * 1000  # en meV
                energies_theory.append(E_n)
            
            ax1.plot(B_values, energies_theory, 'o-', label=f'n = {n}')
        
        ax1.set_xlabel('Campo magnético B (T)')
        ax1.set_ylabel('Energía (meV)')
        ax1.set_title('Niveles de Landau vs Campo Magnético\n(Predicción Teórica)')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # Longitud magnética vs campo
        l_B_values = [np.sqrt(self.hbar / (self.e * B)) * 1e10 for B in B_values]
        ax2.plot(B_values, l_B_values, 'ro-', linewidth=2)
        ax2.set_xlabel('Campo magnético B (T)')
        ax2.set_ylabel('Longitud magnética (Å)')
        ax2.set_title('Longitud Magnética vs Campo')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('magnetic_field_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def comprehensive_analysis(self):
        """Realizar un análisis completo de los datos"""
        print("\n" + "="*50)
        print("ANÁLISIS COMPLETO DE NIVELES DE LANDAU")
        print("="*50)
        
        if len(self.energies_eV) == 0:
            print("Error: No se encontraron datos. Ejecuta primero la simulación C++.")
            return
            
        # 1. Espectro de energías
        print("\n1. Generando espectro de energías...")
        self.plot_energy_spectrum()
        
        # 2. Análisis de degeneración
        print("\n2. Analizando degeneración...")
        energy_gaps = np.diff(self.energies_eV) * 1000  # en meV
        print(f"Separaciones entre niveles (meV): {energy_gaps[:5]}")
        print(f"Separación promedio: {np.mean(energy_gaps):.2f} ± {np.std(energy_gaps):.2f} meV")
        
        # 3. Funciones de onda de algunos estados
        states_to_plot = min(4, len(self.wavefunction_files))
        print(f"\n3. Graficando funciones de onda de los primeros {states_to_plot} estados...")
        
        for i in range(states_to_plot):
            print(f"   Graficando estado {i}...")
            self.plot_wavefunction_2d(i)
            
        # 4. Visualización 3D del estado fundamental
        if len(self.wavefunction_files) > 0:
            print("\n4. Visualización 3D del estado fundamental...")
            self.plot_wavefunction_3d(0)
            
        # 5. Comparación con diferentes campos magnéticos
        print("\n5. Análisis teórico para diferentes campos...")
        self.magnetic_field_comparison()
        
        # 6. Crear animación
        print("\n6. Creando animación...")
        self.create_animation()
        
        print("\n" + "="*50)
        print("ANÁLISIS COMPLETADO")
        print("="*50)
        print("Archivos generados:")
        print("- energy_spectrum.png")
        print("- wavefunction_2d_state_*.png")
        print("- wavefunction_3d_state_*.png") 
        print("- magnetic_field_analysis.png")
        print("- landau_levels_animation.gif")

# Función principal para usar el visualizador
def main():
    """Función principal para ejecutar todas las visualizaciones"""
    # Parámetros que deben coincidir con la simulación C++
    Nx, Ny = 40, 40  # Ajustar según tu simulación
    B = 5.0          # Tesla
    
    # Crear el visualizador
    visualizer = LandauVisualization(Nx=Nx, Ny=Ny, B=B)
    
    # Ejecutar análisis completo
    visualizer.comprehensive_analysis()

if __name__ == "__main__":
    main()
