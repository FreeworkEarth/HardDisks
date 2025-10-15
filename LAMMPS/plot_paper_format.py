#!/usr/bin/env python3
"""
Plot wall oscillations in paper format
Frequency range: 0 to 0.5 Hz as requested
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Load data
print("Loading wall_long.dat...")
data = np.loadtxt('wall_long.dat', skiprows=1)
step = data[:, 0]
wall_x = data[:, 1]
wall_vx = data[:, 2]
n_left = data[:, 3]
n_right = data[:, 4]

# Time array
dt_sampling = 20 * 0.0001  # 20 steps × 0.0001 timestep
time = step * 0.0001

print(f"Loaded {len(wall_x)} data points")
print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} LJ units")

# ============================================================================
# Figure 1: Wall Position and FFT (0-0.5 Hz)
# ============================================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Wall position time series
ax1.plot(time, wall_x, 'b-', linewidth=0.8)
ax1.set_xlabel('Time (LJ units)', fontsize=14)
ax1.set_ylabel('Wall X Position (σ)', fontsize=14)
ax1.set_title('Wall Oscillation - 2D Hard Disk Gas (Paper Configuration)',
              fontsize=15, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.tick_params(labelsize=12)

# Add info box
info_text = f'Simulation time: {time[-1]:.0f} LJ units\n'
info_text += f'Particles: ~144 (72 left, 72 right)\n'
info_text += f'Wall mass: M = 200m'
ax1.text(0.02, 0.98, info_text,
         transform=ax1.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
         fontsize=11)

# Plot 2: FFT Power Spectrum (0 to 0.5 Hz)
signal = wall_x - np.mean(wall_x)
window = np.hanning(len(signal))
signal_windowed = signal * window

N = len(signal)
yf = fft(signal_windowed)
xf = fftfreq(N, dt_sampling)[:N//2]
power = 2.0/N * np.abs(yf[:N//2])

# Find peak frequency
idx_peak = np.argmax(power[1:]) + 1
freq_peak = xf[idx_peak]

# Plot only 0 to 0.5 Hz as requested
freq_max = 0.5
mask = xf <= freq_max

ax2.semilogy(xf[mask], power[mask], 'r-', linewidth=1.5)
ax2.axvline(freq_peak, color='k', linestyle='--', linewidth=2.5,
            label=f'Peak: f = {freq_peak:.5f} Hz')
ax2.set_xlabel('Frequency f (ω/2π, LJ units)', fontsize=14)
ax2.set_ylabel('Power (log scale)', fontsize=14)
ax2.set_title('FFT Power Spectrum (0 to 0.5 Hz)', fontsize=15, fontweight='bold')
ax2.set_xlim(0, freq_max)
ax2.legend(fontsize=12, loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('paper_wall_analysis.png', dpi=300, bbox_inches='tight')
plt.savefig('paper_wall_analysis.pdf', bbox_inches='tight')
print("\nSaved: paper_wall_analysis.png and paper_wall_analysis.pdf")
plt.close()

# ============================================================================
# Figure 2: Detailed 4-Panel View
# ============================================================================

fig = plt.figure(figsize=(16, 12))

# Panel 1: Position
ax1 = plt.subplot(4, 1, 1)
ax1.plot(time, wall_x, 'b-', linewidth=1)
ax1.set_ylabel('Wall X (σ)', fontsize=13, fontweight='bold')
ax1.set_title('Wall Dynamics - Paper Configuration', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.tick_params(labelsize=11)

# Panel 2: Velocity
ax2 = plt.subplot(4, 1, 2)
ax2.plot(time, wall_vx, 'r-', linewidth=1)
ax2.set_ylabel('Wall Velocity (σ/τ)', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.tick_params(labelsize=11)

# Panel 3: Particle distribution
ax3 = plt.subplot(4, 1, 3)
ax3.plot(time, n_left, 'b-', linewidth=1.5, label='Left', alpha=0.8)
ax3.plot(time, n_right, 'r-', linewidth=1.5, label='Right', alpha=0.8)
ax3.set_ylabel('Particle Count', fontsize=13, fontweight='bold')
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.tick_params(labelsize=11)

# Panel 4: FFT (0-0.5 Hz)
ax4 = plt.subplot(4, 1, 4)
ax4.semilogy(xf[mask], power[mask], 'purple', linewidth=1.5)
ax4.axvline(freq_peak, color='k', linestyle='--', linewidth=2,
            label=f'f = {freq_peak:.5f} Hz')
ax4.set_xlabel('Frequency (Hz)', fontsize=13, fontweight='bold')
ax4.set_ylabel('Power', fontsize=13, fontweight='bold')
ax4.set_xlim(0, freq_max)
ax4.legend(fontsize=11)
ax4.grid(True, alpha=0.3)
ax4.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('paper_full_analysis.png', dpi=300, bbox_inches='tight')
plt.savefig('paper_full_analysis.pdf', bbox_inches='tight')
print("Saved: paper_full_analysis.png and paper_full_analysis.pdf")
plt.close()

# ============================================================================
# Print Summary
# ============================================================================

print("\n" + "="*60)
print("PAPER FORMAT ANALYSIS")
print("="*60)

print(f"\nSimulation:")
print(f"  Total time: {time[-1]:.1f} LJ time units")
print(f"  Data points: {len(time):,}")

print(f"\nParticles:")
print(f"  Mean left: {np.mean(n_left):.1f}")
print(f"  Mean right: {np.mean(n_right):.1f}")
print(f"  Total: {int(n_left[0] + n_right[0])}")

print(f"\nWall Oscillation:")
print(f"  Peak frequency: f = {freq_peak:.6f} (ω/2π)")
print(f"  Angular frequency: ω = {2*np.pi*freq_peak:.6f} rad/τ")
print(f"  Period: T = {1/freq_peak:.2f} LJ time units")

# Speed of sound estimate
L0 = 10.0  # Half box length from the working simulation
wavelength = 2 * L0
c_s = 2 * np.pi * freq_peak * wavelength

print(f"\nSpeed of Sound Estimate:")
print(f"  Measured: c_s = {c_s:.4f} σ/τ")
print(f"  (For 2D ideal gas at T=0.5: c_s ≈ 1.0 σ/τ)")

print(f"\nPlots saved with frequency range 0 to 0.5 Hz as requested")
print("="*60)

# Save summary
with open('paper_summary.txt', 'w') as f:
    f.write("PAPER FORMAT ANALYSIS - LAMMPS 2D HARD DISK GAS\n")
    f.write("="*60 + "\n\n")
    f.write(f"Simulation time: {time[-1]:.1f} LJ units\n")
    f.write(f"Particles: ~{int(n_left[0] + n_right[0])} ({np.mean(n_left):.0f} left, {np.mean(n_right):.0f} right)\n")
    f.write(f"\nPeak frequency: {freq_peak:.6f} Hz (ω/2π)\n")
    f.write(f"Period: {1/freq_peak:.2f} LJ time units\n")
    f.write(f"Speed of sound: {c_s:.4f} σ/τ\n")

print("\nSaved: paper_summary.txt")
