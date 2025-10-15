#!/usr/bin/env python3
"""
Analyze energy transfer experiment from LAMMPS
Track energy flow: KE_right → SpringPE → KE_left
Compute efficiency and transfer dynamics
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

def load_energy_data(filename):
    """Load energy transfer data from LAMMPS output"""
    data = np.loadtxt(filename, skiprows=1)

    ke_particles = data[:, 0]
    ke_wall = data[:, 1]
    ke_total = data[:, 2]
    temp = data[:, 3]
    spring_pe = data[:, 4]
    total_e = data[:, 5]
    ke_left = data[:, 6]
    ke_right = data[:, 7]

    return {
        'ke_particles': ke_particles,
        'ke_wall': ke_wall,
        'ke_total': ke_total,
        'temp': temp,
        'spring_pe': spring_pe,
        'total_e': total_e,
        'ke_left': ke_left,
        'ke_right': ke_right
    }

def load_wall_data(filename):
    """Load wall position data"""
    data = np.loadtxt(filename, skiprows=1)

    step = data[:, 0]
    wall_x = data[:, 1]
    wall_disp = data[:, 2]
    n_left = data[:, 3]
    n_right = data[:, 4]

    return step, wall_x, wall_disp, n_left, n_right

def analyze_energy_transfer(energy_data, time):
    """Analyze energy transfer dynamics"""

    # Energy transferred from right to left
    ke_initial_right = energy_data['ke_right'][0]
    ke_final_left = energy_data['ke_left'][-1]

    # Maximum spring potential energy (peak energy storage)
    spring_pe_max = np.max(energy_data['spring_pe'])
    spring_pe_max_idx = np.argmax(energy_data['spring_pe'])
    time_max_spring = time[spring_pe_max_idx]

    # Energy efficiency: how much energy made it to left chamber
    efficiency = ke_final_left / ke_initial_right if ke_initial_right > 0 else 0

    # Energy conservation check
    e_drift = np.std(energy_data['total_e']) / np.mean(energy_data['total_e'])

    results = {
        'ke_initial_right': ke_initial_right,
        'ke_final_left': ke_final_left,
        'spring_pe_max': spring_pe_max,
        'time_max_spring': time_max_spring,
        'efficiency': efficiency,
        'energy_drift_relative': e_drift,
        'total_energy_mean': np.mean(energy_data['total_e']),
        'total_energy_std': np.std(energy_data['total_e'])
    }

    return results

def plot_energy_transfer(time, energy_data, wall_x, n_left, n_right, output_prefix):
    """Create comprehensive energy transfer plots"""

    fig = plt.figure(figsize=(16, 12))

    # 1. Energy components vs time
    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(time, energy_data['ke_left'], 'b-', label='KE Left', linewidth=1.5)
    ax1.plot(time, energy_data['ke_right'], 'r-', label='KE Right', linewidth=1.5)
    ax1.plot(time, energy_data['spring_pe'], 'g-', label='Spring PE', linewidth=1.5)
    ax1.plot(time, energy_data['ke_wall'], 'm-', label='KE Wall', linewidth=1)
    ax1.set_xlabel('Time (LJ units)')
    ax1.set_ylabel('Energy (ε)')
    ax1.set_title('Energy Transfer Dynamics')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Total energy conservation
    ax2 = plt.subplot(3, 2, 2)
    ax2.plot(time, energy_data['total_e'], 'k-', linewidth=1)
    ax2.set_xlabel('Time (LJ units)')
    ax2.set_ylabel('Total Energy (ε)')
    ax2.set_title(f'Energy Conservation (σ/μ = {np.std(energy_data["total_e"])/np.mean(energy_data["total_e"]):.2e})')
    ax2.grid(True, alpha=0.3)

    # 3. Wall position
    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(time, wall_x, 'b-', linewidth=1)
    ax3.set_xlabel('Time (LJ units)')
    ax3.set_ylabel('Wall Position (σ)')
    ax3.set_title('Divider Wall Position')
    ax3.grid(True, alpha=0.3)

    # 4. Particle counts
    ax4 = plt.subplot(3, 2, 4)
    ax4.plot(time, n_left, 'b-', label='Left', linewidth=1.5)
    ax4.plot(time, n_right, 'r-', label='Right', linewidth=1.5)
    ax4.set_xlabel('Time (LJ units)')
    ax4.set_ylabel('Particle Count')
    ax4.set_title('Particles per Chamber')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Temperature
    ax5 = plt.subplot(3, 2, 5)
    ax5.plot(time, energy_data['temp'], 'orange', linewidth=1)
    ax5.set_xlabel('Time (LJ units)')
    ax5.set_ylabel('Temperature (kB units)')
    ax5.set_title('System Temperature')
    ax5.grid(True, alpha=0.3)

    # 6. Energy flow diagram (snapshot at different times)
    ax6 = plt.subplot(3, 2, 6)

    # Select 4 time snapshots
    n_snapshots = 4
    indices = np.linspace(0, len(time)-1, n_snapshots, dtype=int)

    bar_width = 0.2
    x_pos = np.arange(n_snapshots)

    ke_left_snap = energy_data['ke_left'][indices]
    ke_right_snap = energy_data['ke_right'][indices]
    spring_pe_snap = energy_data['spring_pe'][indices]
    ke_wall_snap = energy_data['ke_wall'][indices]

    ax6.bar(x_pos - 1.5*bar_width, ke_left_snap, bar_width, label='KE Left', color='b')
    ax6.bar(x_pos - 0.5*bar_width, ke_right_snap, bar_width, label='KE Right', color='r')
    ax6.bar(x_pos + 0.5*bar_width, spring_pe_snap, bar_width, label='Spring PE', color='g')
    ax6.bar(x_pos + 1.5*bar_width, ke_wall_snap, bar_width, label='KE Wall', color='m')

    ax6.set_xlabel('Time Snapshot')
    ax6.set_ylabel('Energy (ε)')
    ax6.set_title('Energy Distribution Evolution')
    ax6.set_xticks(x_pos)
    ax6.set_xticklabels([f't={time[i]:.1f}' for i in indices], rotation=45)
    ax6.legend()
    ax6.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_analysis.png', dpi=150)
    plt.savefig(f'{output_prefix}_analysis.pdf')
    print(f"Saved plots: {output_prefix}_analysis.png/pdf")

    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Analyze LAMMPS energy transfer experiment')
    parser.add_argument('--energy-file', default='energy_transfer.dat',
                        help='Energy data file')
    parser.add_argument('--wall-file', default='wall_energy_transfer.dat',
                        help='Wall position data file')
    parser.add_argument('--dt', type=float, default=1000, help='Logging interval (timesteps)')
    parser.add_argument('--dt-timestep', type=float, default=0.000125,
                        help='LAMMPS timestep (LJ units)')
    parser.add_argument('--output', default='energy_transfer',
                        help='Output file prefix')

    args = parser.parse_args()

    # Load data
    print(f"Loading energy data from {args.energy_file}...")
    energy_data = load_energy_data(args.energy_file)

    print(f"Loading wall data from {args.wall_file}...")
    step, wall_x, wall_disp, n_left, n_right = load_wall_data(args.wall_file)

    # Compute time array
    dt_real = args.dt * args.dt_timestep
    time = step * args.dt_timestep

    print(f"Data points: {len(time)}")
    print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} (LJ units)")

    # Analyze transfer
    results = analyze_energy_transfer(energy_data, time)

    print(f"\n=== Energy Transfer Analysis ===")
    print(f"Initial KE (right chamber): {results['ke_initial_right']:.6f} ε")
    print(f"Final KE (left chamber): {results['ke_final_left']:.6f} ε")
    print(f"Maximum spring PE: {results['spring_pe_max']:.6f} ε")
    print(f"Time of max spring PE: {results['time_max_spring']:.6f} LJ units")
    print(f"Transfer efficiency: {results['efficiency']*100:.2f}%")

    print(f"\n=== Energy Conservation ===")
    print(f"Total energy (mean): {results['total_energy_mean']:.6f} ε")
    print(f"Total energy (std): {results['total_energy_std']:.6e} ε")
    print(f"Relative drift: {results['energy_drift_relative']:.6e}")

    # Particle transfer
    print(f"\n=== Particle Transfer ===")
    print(f"Initial left: {n_left[0]:.0f}, right: {n_right[0]:.0f}")
    print(f"Final left: {n_left[-1]:.0f}, right: {n_right[-1]:.0f}")
    print(f"Net transfer (right→left): {n_left[-1] - n_left[0]:.0f} particles")

    # Create plots
    plot_energy_transfer(time, energy_data, wall_x, n_left, n_right, args.output)

    # Save results
    with open(f'{args.output}_results.txt', 'w') as f:
        f.write("=== LAMMPS Energy Transfer Analysis ===\n\n")
        for key, val in results.items():
            f.write(f"{key}: {val:.8e}\n")

        f.write(f"\n=== Particle Transfer ===\n")
        f.write(f"initial_left: {n_left[0]:.0f}\n")
        f.write(f"initial_right: {n_right[0]:.0f}\n")
        f.write(f"final_left: {n_left[-1]:.0f}\n")
        f.write(f"final_right: {n_right[-1]:.0f}\n")
        f.write(f"net_transfer: {n_left[-1] - n_left[0]:.0f}\n")

    print(f"\nResults saved to {args.output}_results.txt")

if __name__ == '__main__':
    main()
