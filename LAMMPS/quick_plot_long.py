#!/usr/bin/env python3
"""Quick plot of wall oscillations without full dependencies"""
import sys

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from scipy.fft import fft, fftfreq
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("\nInstall with:")
    print("  pip3 install --break-system-packages numpy scipy matplotlib")
    print("  OR: conda install numpy scipy matplotlib")
    sys.exit(1)

# Load data
data = np.loadtxt('wall_long.dat', skiprows=1)
step = data[:, 0]
wall_x = data[:, 1]

# Time array
dt = 20 * 0.0001  # sampling_interval * timestep
time = (step - step[0]) * 0.0001

print(f"Loaded {len(wall_x)} data points")
print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} LJ units")

# Create figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

# Plot 1: Wall position vs time
ax1.plot(time, wall_x, 'b-', linewidth=0.5)
ax1.set_xlabel('Time (LJ units)', fontsize=12)
ax1.set_ylabel('Wall X Position (σ)', fontsize=12)
ax1.set_title('Wall Oscillations in 2D Hard Disk Gas', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)

# Plot 2: FFT Power Spectrum
signal = wall_x - np.mean(wall_x)
window = np.hanning(len(signal))
signal_windowed = signal * window

N = len(signal)
yf = fft(signal_windowed)
xf = fftfreq(N, dt)[:N//2]
power = 2.0/N * np.abs(yf[:N//2])

# Find peak frequency
idx_peak = np.argmax(power[1:]) + 1
freq_peak = xf[idx_peak]

ax2.semilogy(xf, power, 'r-', linewidth=1)
ax2.axvline(freq_peak, color='k', linestyle='--', linewidth=2,
            label=f'Peak: {freq_peak:.4f} (ω/2π)')
ax2.set_xlabel('Frequency (ω/2π, LJ units)', fontsize=12)
ax2.set_ylabel('Power (log scale)', fontsize=12)
ax2.set_title('FFT Power Spectrum', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, xf[len(xf)//10])

plt.tight_layout()
plt.savefig('wall_analysis_long.png', dpi=150, bbox_inches='tight')
plt.savefig('wall_analysis_long.pdf', bbox_inches='tight')
print(f"\nSaved: wall_analysis_long.png and wall_analysis_long.pdf")

# Print results
print(f"\n=== Analysis Results ===")
print(f"Dominant frequency: {freq_peak:.6f} (ω/2π)")
print(f"Angular frequency: {2*np.pi*freq_peak:.6f} (ω)")
print(f"Period: {1.0/freq_peak:.6f} LJ time units")

# Speed of sound estimate
L0 = 10.0
wavelength = 2 * L0
c_s_measured = 2 * np.pi * freq_peak * wavelength
c_s_theory = np.sqrt(2.0)  # 2D ideal gas at T=1

print(f"\nSpeed of sound (theoretical, T=1): {c_s_theory:.6f}")
print(f"Speed of sound (measured): {c_s_measured:.6f}")
print(f"Note: Simulation ran at T=0.5, so expect c_s ~ {np.sqrt(2*0.5):.6f}")

# plt.show()
