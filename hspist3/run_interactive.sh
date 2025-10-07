#!/bin/bash
# Interactive Simulation Launcher

echo "Hard Disk Simulation - Interactive Mode"
echo "======================================"
echo ""
echo "Available modes:"
echo "1. Basic interactive simulation"
echo "2. Custom parameters"
echo "3. GUI configuration"
echo "4. Exit"
echo ""

read -p "Choose mode (1-4): " choice

case $choice in
    1)
        echo "Starting basic interactive simulation..."
        ./00ALLINONE --no-experiments
        ;;
    2)
        echo "Custom parameters mode"
        read -p "Number of particles (default 1000): " particles
        read -p "L0 units (default 20.0): " l0
        read -p "Height units (default 10.0): " height
        read -p "Number of walls (default 1): " walls
        
        # Set defaults if empty
        particles=${particles:-1000}
        l0=${l0:-20.0}
        height=${height:-10.0}
        walls=${walls:-1}
        
        echo "Starting simulation with: particles=$particles, l0=$l0, height=$height, walls=$walls"
        ./00ALLINONE --no-experiments --particles=$particles --l0=$l0 --height=$height --num-walls=$walls
        ;;
    3)
        echo "Starting GUI configuration..."
        if command -v python3 &> /dev/null; then
            python3 experiment_gui.py
        else
            echo "Python3 not found. Building C GUI..."
            make -f Makefile.gui run
        fi
        ;;
    4)
        echo "Goodbye!"
        exit 0
        ;;
    *)
        echo "Invalid choice. Starting basic simulation..."
        ./00ALLINONE --no-experiments
        ;;
esac
