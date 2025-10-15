#!/bin/bash
# Master script to run all LAMMPS hard disk experiments
# Usage: ./run_experiments.sh

set -e  # Exit on error

echo "========================================"
echo "LAMMPS Hard Disk Gas Experiments"
echo "========================================"
echo ""

# Check if LAMMPS is installed
LAMMPS_CMD=""
if command -v lammps &> /dev/null; then
    LAMMPS_CMD="lammps"
elif command -v lmp_serial &> /dev/null; then
    LAMMPS_CMD="lmp_serial"
elif command -v lmp_mpi &> /dev/null; then
    LAMMPS_CMD="lmp_mpi"
else
    echo "ERROR: LAMMPS not found in PATH"
    echo "Please install LAMMPS:"
    echo "  conda install -c conda-forge lammps"
    echo "  OR brew install lammps"
    exit 1
fi

echo "Using LAMMPS command: $LAMMPS_CMD"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python3 not found"
    exit 1
fi

# Create output directory
mkdir -p results
cd results

echo "1. Running basic hard disk simulation..."
echo "----------------------------------------"
$LAMMPS_CMD -in ../in.harddisk -log log.harddisk
echo "✓ Basic simulation complete"
echo "  Output: traj.lammpstrj, wall_position.dat, energy.dat"
echo ""

echo "2. Running speed of sound experiment..."
echo "----------------------------------------"
$LAMMPS_CMD -in ../in.speed_of_sound -log log.speed_of_sound
echo "✓ Speed of sound simulation complete"
echo ""

echo "3. Analyzing wall oscillations..."
echo "----------------------------------------"
python3 ../analyze_oscillations.py wall_oscillation.dat \
    --dt 10 \
    --dt-timestep 0.00005 \
    --L0 7.5 \
    --wall-mass-factor 20.0 \
    --output speed_of_sound

echo "✓ Oscillation analysis complete"
echo "  Results: speed_of_sound_results.txt, speed_of_sound_analysis.png"
echo ""

echo "4. Running energy transfer experiment..."
echo "----------------------------------------"
$LAMMPS_CMD -in ../in.energy_transfer -log log.energy_transfer
echo "✓ Energy transfer simulation complete"
echo ""

echo "5. Analyzing energy transfer..."
echo "----------------------------------------"
python3 ../analyze_energy_transfer.py \
    --energy-file energy_transfer.dat \
    --wall-file wall_energy_transfer.dat \
    --dt 1000 \
    --dt-timestep 0.000125 \
    --output energy_transfer

echo "✓ Energy transfer analysis complete"
echo "  Results: energy_transfer_results.txt, energy_transfer_analysis.png"
echo ""

echo "========================================"
echo "All experiments complete!"
echo "========================================"
echo ""
echo "Visualization:"
echo "  ovito traj.lammpstrj               # Basic simulation"
echo "  ovito traj_speed_of_sound.lammpstrj   # Speed of sound"
echo "  ovito traj_energy_transfer.lammpstrj  # Energy transfer"
echo ""
echo "Results summary:"
echo "  cat speed_of_sound_results.txt"
echo "  cat energy_transfer_results.txt"
echo ""
echo "Plots generated:"
echo "  speed_of_sound_analysis.png"
echo "  energy_transfer_analysis.png"
