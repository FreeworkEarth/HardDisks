#!/usr/bin/env python3
"""
Analyze wall oscillations from LAMMPS speed of sound experiment
Performs FFT to extract oscillation frequency
Replicates analysis from hspist3/experiments_speed_of_sound
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq
import argparse

def load_wall_data(filename):
    """Load wall position time series from LAMMPS output"""
    data = np.loadtxt(filename, skiprows=1)

    step = data[:, 0]
    wall_x = data[:, 1]
    wall_vx = data[:, 2] if data.shape[1] > 2 else None

    return step, wall_x, wall_vx

def analyze_frequency(time, signal_data, dt):
    """Perform FFT analysis to extract dominant frequency"""

    # Remove DC component (mean)
    signal_data = signal_data - np.mean(signal_data)

    # Apply windowing to reduce spectral leakage
    window = np.hanning(len(signal_data))
    signal_windowed = signal_data * window

    # Compute FFT
    N = len(signal_data)
    yf = fft(signal_windowed)
    xf = fftfreq(N, dt)[:N//2]

    # Power spectrum
    power = 2.0/N * np.abs(yf[:N//2])

    # Find dominant frequency (exclude DC component)
    idx_peak = np.argmax(power[1:]) + 1
    freq_dominant = xf[idx_peak]

    return xf, power, freq_dominant

def plot_analysis(time, wall_x, freq, power, freq_dominant, output_prefix):
    """Create analysis plots"""

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Time series
    axes[0].plot(time, wall_x, 'b-', linewidth=0.5)
    axes[0].set_xlabel('Time (LJ units)')
    axes[0].set_ylabel('Wall Position (σ)')
    axes[0].set_title('Wall Position vs Time')
    axes[0].grid(True, alpha=0.3)

    # Power spectrum
    axes[1].semilogy(freq, power, 'r-', linewidth=1)
    axes[1].axvline(freq_dominant, color='k', linestyle='--',
                    label=f'Peak: {freq_dominant:.4f} (ω/2π)')
    axes[1].set_xlabel('Frequency (ω/2π, LJ units)')
    axes[1].set_ylabel('Power')
    axes[1].set_title('Power Spectrum (FFT)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xlim(0, freq[len(freq)//10])  # Show first 10% of spectrum

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_analysis.png', dpi=150)
    plt.savefig(f'{output_prefix}_analysis.pdf')
    print(f"Saved plots: {output_prefix}_analysis.png/pdf")

    plt.show()

def calculate_speed_of_sound(freq_dominant, L0, wall_mass_factor):
    """
    Calculate speed of sound from oscillation frequency

    For a piston oscillating in a 2D ideal gas:
    ω = sqrt(γ * P * A / (M * L))

    where γ = 2 (2D), P = pressure, A = cross-section, M = wall mass, L = box half-width

    Speed of sound: c_s = sqrt(γ * P / ρ) = sqrt(γ * k_B * T / m)
    """

    # In LJ units with T=1, m=1, k_B=1, gamma=2:
    # c_s = sqrt(2 * 1 / 1) = sqrt(2) ≈ 1.414

    c_s_theoretical = np.sqrt(2.0)

    # Oscillation frequency relates to sound speed:
    # For fundamental mode: f ∝ c_s / L
    wavelength = 2 * L0  # Fundamental mode
    c_s_measured = 2 * np.pi * freq_dominant * wavelength

    return c_s_theoretical, c_s_measured

def main():
    parser = argparse.ArgumentParser(description='Analyze LAMMPS wall oscillation data')
    parser.add_argument('filename', help='Wall oscillation data file (wall_oscillation.dat)')
    parser.add_argument('--dt', type=float, default=10, help='Sampling interval (timesteps)')
    parser.add_argument('--dt-timestep', type=float, default=0.00005,
                        help='LAMMPS timestep (LJ units)')
    parser.add_argument('--L0', type=float, default=7.5, help='Box half-width (sigma)')
    parser.add_argument('--wall-mass-factor', type=float, default=20.0,
                        help='Wall mass / particle mass')
    parser.add_argument('--output', default='wall_oscillation', help='Output file prefix')
    parser.add_argument('--skip', type=int, default=0,
                        help='Skip initial timesteps (equilibration)')

    args = parser.parse_args()

    # Load data
    print(f"Loading data from {args.filename}...")
    step, wall_x, wall_vx = load_wall_data(args.filename)

    # Skip equilibration phase
    if args.skip > 0:
        mask = step >= args.skip
        step = step[mask]
        wall_x = wall_x[mask]
        if wall_vx is not None:
            wall_vx = wall_vx[mask]

    # Compute time array
    dt_real = args.dt * args.dt_timestep  # Real time between samples
    time = (step - step[0]) * args.dt_timestep

    print(f"Data points: {len(wall_x)}")
    print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} (LJ units)")
    print(f"Sampling dt: {dt_real:.6f} (LJ units)")

    # Frequency analysis
    freq, power, freq_dominant = analyze_frequency(time, wall_x, dt_real)

    print(f"\n=== Frequency Analysis ===")
    print(f"Dominant frequency: {freq_dominant:.6f} (ω/2π, LJ units)")
    print(f"Angular frequency: {2*np.pi*freq_dominant:.6f} (ω, LJ units)")
    print(f"Period: {1.0/freq_dominant:.6f} (LJ time units)")

    # Speed of sound estimate
    c_s_theory, c_s_measured = calculate_speed_of_sound(freq_dominant, args.L0,
                                                         args.wall_mass_factor)

    print(f"\n=== Speed of Sound ===")
    print(f"Theoretical (2D ideal gas, T=1): {c_s_theory:.6f} (LJ units)")
    print(f"Measured (from oscillation): {c_s_measured:.6f} (LJ units)")
    print(f"Ratio (measured/theory): {c_s_measured/c_s_theory:.4f}")

    # Statistics
    print(f"\n=== Wall Motion Statistics ===")
    print(f"Mean position: {np.mean(wall_x):.6f} σ")
    print(f"Std deviation: {np.std(wall_x):.6f} σ")
    print(f"Amplitude (peak-to-peak): {np.ptp(wall_x):.6f} σ")

    # Create plots
    plot_analysis(time, wall_x, freq, power, freq_dominant, args.output)

    # Save results
    results = {
        'frequency_hz': freq_dominant,
        'angular_frequency': 2*np.pi*freq_dominant,
        'period': 1.0/freq_dominant,
        'c_s_theoretical': c_s_theory,
        'c_s_measured': c_s_measured,
        'mean_position': np.mean(wall_x),
        'amplitude_std': np.std(wall_x),
        'amplitude_p2p': np.ptp(wall_x)
    }

    with open(f'{args.output}_results.txt', 'w') as f:
        f.write("=== LAMMPS Speed of Sound Analysis ===\n\n")
        for key, val in results.items():
            f.write(f"{key}: {val:.8f}\n")

    print(f"\nResults saved to {args.output}_results.txt")

if __name__ == '__main__':
    main()
