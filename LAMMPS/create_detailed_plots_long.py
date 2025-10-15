#!/usr/bin/env python3
"""
Create detailed plots from LAMMPS hard disk simulation
Shows wall position, velocity, particle counts, and FFT analysis
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy import signal

# Load data
print("Loading data...")
data = np.loadtxt('wall_long.dat', skiprows=1)
step = data[:, 0]
wall_x = data[:, 1]
wall_vx = data[:, 2]
n_left = data[:, 3]
n_right = data[:, 4]

# Time array
dt_sampling = 20 * 0.0001  # 20 steps * 0.0001 timestep
time = (step - step[0]) * 0.0001

print(f"Loaded {len(wall_x)} data points")
print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} LJ units")
print(f"Total simulation steps: {step[-1] - step[0]:.0f}")

# ============================================================================
# FIGURE 1: Comprehensive Dynamics
# ============================================================================

fig = plt.figure(figsize=(16, 12))

# 1. Wall position
ax1 = plt.subplot(4, 1, 1)
ax1.plot(time, wall_x, 'b-', linewidth=1, alpha=0.8)
ax1.axhline(0, color='k', linestyle='--', alpha=0.3)
ax1.set_ylabel('Wall Position (σ)', fontsize=12, fontweight='bold')
ax1.set_title('2D Hard Disk Gas: Wall Oscillation Dynamics', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(time[0], time[-1])

# Add statistics
mean_x = np.mean(wall_x)
std_x = np.std(wall_x)
ax1.text(0.02, 0.98, f'Mean: {mean_x:.3f} σ\nStd: {std_x:.3f} σ',
         transform=ax1.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
         fontsize=10)

# 2. Wall velocity
ax2 = plt.subplot(4, 1, 2)
ax2.plot(time, wall_vx, 'r-', linewidth=1, alpha=0.8)
ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
ax2.set_ylabel('Wall Velocity (σ/τ)', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(time[0], time[-1])

mean_vx = np.mean(wall_vx)
std_vx = np.std(wall_vx)
ax2.text(0.02, 0.98, f'Mean: {mean_vx:.4f} σ/τ\nStd: {std_vx:.4f} σ/τ',
         transform=ax2.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5),
         fontsize=10)

# 3. Particle counts
ax3 = plt.subplot(4, 1, 3)
ax3.plot(time, n_left, 'b-', linewidth=1.5, label='Left chamber', alpha=0.8)
ax3.plot(time, n_right, 'r-', linewidth=1.5, label='Right chamber', alpha=0.8)
ax3.plot(time, n_left + n_right, 'k--', linewidth=1, label='Total', alpha=0.5)
ax3.set_ylabel('Particle Count', fontsize=12, fontweight='bold')
ax3.legend(loc='upper right', fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(time[0], time[-1])

# 4. Phase space (position vs velocity)
ax4 = plt.subplot(4, 1, 4)
# Color by time for trajectory visualization
colors = plt.cm.viridis(np.linspace(0, 1, len(time)))
ax4.scatter(wall_x, wall_vx, c=time, s=1, cmap='viridis', alpha=0.5)
ax4.set_xlabel('Wall Position (σ)', fontsize=12, fontweight='bold')
ax4.set_ylabel('Wall Velocity (σ/τ)', fontsize=12, fontweight='bold')
ax4.set_title('Phase Space Trajectory', fontsize=11, fontweight='bold')
ax4.grid(True, alpha=0.3)
cbar = plt.colorbar(ax4.collections[0], ax=ax4)
cbar.set_label('Time (LJ units)', fontsize=10)

plt.tight_layout()
plt.savefig('wall_dynamics_long.png', dpi=200, bbox_inches='tight')
plt.savefig('wall_dynamics_long.pdf', bbox_inches='tight')
print("\nSaved: wall_dynamics_long.png and wall_dynamics_long.pdf")
plt.close()

# ============================================================================
# FIGURE 2: FFT Analysis
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Remove DC component
signal_x = wall_x - np.mean(wall_x)
signal_vx = wall_vx - np.mean(wall_vx)

# Window and FFT
window = np.hanning(len(signal_x))
signal_x_win = signal_x * window
signal_vx_win = signal_vx * window

N = len(signal_x)
freq = fftfreq(N, dt_sampling)[:N//2]

fft_x = fft(signal_x_win)
power_x = 2.0/N * np.abs(fft_x[:N//2])

fft_vx = fft(signal_vx_win)
power_vx = 2.0/N * np.abs(fft_vx[:N//2])

# Find peaks
idx_peak_x = np.argmax(power_x[1:]) + 1
freq_peak_x = freq[idx_peak_x]

idx_peak_vx = np.argmax(power_vx[1:]) + 1
freq_peak_vx = freq[idx_peak_vx]

# Plot 1: Position signal
axes[0, 0].plot(time, wall_x, 'b-', linewidth=0.8, alpha=0.7)
axes[0, 0].set_xlabel('Time (LJ units)', fontsize=11)
axes[0, 0].set_ylabel('Wall Position (σ)', fontsize=11)
axes[0, 0].set_title('Wall Position Time Series', fontsize=12, fontweight='bold')
axes[0, 0].grid(True, alpha=0.3)

# Plot 2: Position FFT
axes[0, 1].semilogy(freq, power_x, 'b-', linewidth=1)
axes[0, 1].axvline(freq_peak_x, color='r', linestyle='--', linewidth=2,
                   label=f'Peak: {freq_peak_x:.5f} Hz')
axes[0, 1].set_xlabel('Frequency (ω/2π, LJ units)', fontsize=11)
axes[0, 1].set_ylabel('Power (log)', fontsize=11)
axes[0, 1].set_title('Position Power Spectrum', fontsize=12, fontweight='bold')
axes[0, 1].legend(fontsize=10)
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_xlim(0, freq[len(freq)//10])

# Plot 3: Velocity signal
axes[1, 0].plot(time, wall_vx, 'r-', linewidth=0.8, alpha=0.7)
axes[1, 0].set_xlabel('Time (LJ units)', fontsize=11)
axes[1, 0].set_ylabel('Wall Velocity (σ/τ)', fontsize=11)
axes[1, 0].set_title('Wall Velocity Time Series', fontsize=12, fontweight='bold')
axes[1, 0].grid(True, alpha=0.3)

# Plot 4: Velocity FFT
axes[1, 1].semilogy(freq, power_vx, 'r-', linewidth=1)
axes[1, 1].axvline(freq_peak_vx, color='b', linestyle='--', linewidth=2,
                   label=f'Peak: {freq_peak_vx:.5f} Hz')
axes[1, 1].set_xlabel('Frequency (ω/2π, LJ units)', fontsize=11)
axes[1, 1].set_ylabel('Power (log)', fontsize=11)
axes[1, 1].set_title('Velocity Power Spectrum', fontsize=12, fontweight='bold')
axes[1, 1].legend(fontsize=10)
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_xlim(0, freq[len(freq)//10])

plt.tight_layout()
plt.savefig('wall_fft_long.png', dpi=200, bbox_inches='tight')
plt.savefig('wall_fft_long.pdf', bbox_inches='tight')
print("Saved: wall_fft_long.png and wall_fft_long.pdf")
plt.close()

# ============================================================================
# FIGURE 3: Zoomed view of oscillations
# ============================================================================

fig, axes = plt.subplots(2, 1, figsize=(14, 8))

# Show first 5 periods
period = 1.0 / freq_peak_x
t_zoom = 5 * period
mask = time <= t_zoom

axes[0].plot(time[mask], wall_x[mask], 'b-', linewidth=2)
axes[0].set_ylabel('Wall Position (σ)', fontsize=12, fontweight='bold')
axes[0].set_title(f'First 5 Oscillation Periods (T ≈ {period:.2f} LJ units)',
                  fontsize=13, fontweight='bold')
axes[0].grid(True, alpha=0.3)

axes[1].plot(time[mask], wall_vx[mask], 'r-', linewidth=2)
axes[1].set_xlabel('Time (LJ units)', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Wall Velocity (σ/τ)', fontsize=12, fontweight='bold')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('wall_oscillations_long_zoom.png', dpi=200, bbox_inches='tight')
plt.savefig('wall_oscillations_long_zoom.pdf', bbox_inches='tight')
print("Saved: wall_oscillations_long_zoom.png and wall_oscillations_long_zoom.pdf")
plt.close()

# ============================================================================
# Print Summary
# ============================================================================

print("\n" + "="*60)
print("SIMULATION ANALYSIS SUMMARY")
print("="*60)

print(f"\nSimulation Parameters:")
print(f"  Total steps: {step[-1] - step[0]:.0f}")
print(f"  Time span: {time[-1]:.2f} LJ time units")
print(f"  Data points: {len(time)}")
print(f"  Sampling interval: {dt_sampling:.6f} LJ units")

print(f"\nWall Dynamics:")
print(f"  Mean position: {mean_x:.4f} ± {std_x:.4f} σ")
print(f"  Position range: [{np.min(wall_x):.4f}, {np.max(wall_x):.4f}] σ")
print(f"  Mean velocity: {mean_vx:.6f} ± {std_vx:.4f} σ/τ")
print(f"  Velocity range: [{np.min(wall_vx):.4f}, {np.max(wall_vx):.4f}] σ/τ")

print(f"\nParticle Distribution:")
print(f"  Mean left: {np.mean(n_left):.1f} particles")
print(f"  Mean right: {np.mean(n_right):.1f} particles")
print(f"  Total particles: {int(n_left[0] + n_right[0])}")

print(f"\nOscillation Analysis:")
print(f"  Fundamental frequency: {freq_peak_x:.6f} (ω/2π)")
print(f"  Angular frequency: {2*np.pi*freq_peak_x:.6f} (ω)")
print(f"  Period: {period:.4f} LJ time units")
print(f"  Number of cycles in data: {time[-1]/period:.1f}")

# Speed of sound estimates
L0 = 10.0
wavelength = 2 * L0
c_s_measured = 2 * np.pi * freq_peak_x * wavelength
c_s_theory_T1 = np.sqrt(2.0)
c_s_theory_T05 = np.sqrt(2 * 0.5)

print(f"\nSpeed of Sound Estimates:")
print(f"  Theoretical (2D ideal gas, T=1.0): {c_s_theory_T1:.4f} σ/τ")
print(f"  Theoretical (2D ideal gas, T=0.5): {c_s_theory_T05:.4f} σ/τ")
print(f"  Measured from oscillation: {c_s_measured:.4f} σ/τ")
print(f"  Ratio (measured/expected): {c_s_measured/c_s_theory_T05:.2f}")

print(f"\n" + "="*60)
print("PLOTS GENERATED:")
print("="*60)
print("  1. wall_analysis_long.png/pdf - Quick FFT analysis")
print("  2. wall_dynamics_long.png/pdf - Comprehensive dynamics")
print("  3. wall_fft_long.png/pdf - Detailed FFT analysis")
print("  4. wall_oscillations_long_zoom.png/pdf - First 5 periods")
print("="*60)

# Save summary to text file
with open('simulation_summary_long.txt', 'w') as f:
    f.write("="*60 + "\n")
    f.write("LAMMPS 2D HARD DISK SIMULATION - ANALYSIS SUMMARY\n")
    f.write("="*60 + "\n\n")

    f.write(f"Simulation Parameters:\n")
    f.write(f"  Total steps: {step[-1] - step[0]:.0f}\n")
    f.write(f"  Time span: {time[-1]:.2f} LJ time units\n")
    f.write(f"  Data points: {len(time)}\n\n")

    f.write(f"Wall Dynamics:\n")
    f.write(f"  Mean position: {mean_x:.4f} ± {std_x:.4f} σ\n")
    f.write(f"  Mean velocity: {mean_vx:.6f} ± {std_vx:.4f} σ/τ\n\n")

    f.write(f"Oscillation Analysis:\n")
    f.write(f"  Frequency: {freq_peak_x:.6f} (ω/2π)\n")
    f.write(f"  Period: {period:.4f} LJ time units\n")
    f.write(f"  Measured c_s: {c_s_measured:.4f} σ/τ\n\n")

    f.write(f"Particle Distribution:\n")
    f.write(f"  Mean left: {np.mean(n_left):.1f}\n")
    f.write(f"  Mean right: {np.mean(n_right):.1f}\n")

print("\nSaved: simulation_summary_long.txt")
print("\nAll plots generated successfully!")
