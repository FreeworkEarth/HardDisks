#!/usr/bin/env python3
"""
Experiment Configuration GUI for Hard Disk Simulation
Provides a graphical interface for setting up complex experiments
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
import subprocess
import os

class ExperimentConfigGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Hard Disk Simulation - Experiment Configuration")
        self.root.geometry("800x600")
        
        # Configuration variables
        self.config = {
            'experiment_type': 'speed_of_sound',
            'num_walls': 1,
            'wall_positions': [0.0],
            'wall_masses': [200.0],
            'protocol': 'sigmoidal',
            'particles': 1000,
            'l0_units': 20.0,
            'height_units': 10.0,
            'wall_thickness': 1.0,
            'repeats': 1,
            'steps': 3000
        }
        
        self.create_widgets()
        
    def create_widgets(self):
        # Main notebook for tabs
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Experiment Type Tab
        self.create_experiment_tab(notebook)
        
        # Wall Configuration Tab
        self.create_wall_tab(notebook)
        
        # Protocol Tab
        self.create_protocol_tab(notebook)
        
        # Simulation Parameters Tab
        self.create_simulation_tab(notebook)
        
        # Control buttons
        self.create_control_buttons()
        
    def create_experiment_tab(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Experiment Type")
        
        # Experiment type selection
        ttk.Label(frame, text="Experiment Type:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.experiment_var = tk.StringVar(value='speed_of_sound')
        experiment_combo = ttk.Combobox(frame, textvariable=self.experiment_var, width=30)
        experiment_combo['values'] = [
            'speed_of_sound',
            'energy_transfer', 
            'multi_wall',
            'szilard_engine',
            'atp_synthase',
            'protocol_optimization',
            'dynamic_ib',
            'wall_mid'
        ]
        experiment_combo.grid(row=0, column=1, sticky='w', padx=5, pady=5)
        
        # Experiment descriptions
        descriptions = {
            'speed_of_sound': 'Standard speed of sound experiment with single wall',
            'energy_transfer': 'Energy transfer efficiency from right piston to left wall',
            'multi_wall': 'Multi-wall system analysis (1-5 walls)',
            'szilard_engine': 'Szilard engine entropy experiments',
            'atp_synthase': 'ATP synthase analog experiments',
            'protocol_optimization': 'Find optimal piston protocols',
            'dynamic_ib': 'Dynamic information bottleneck analysis',
            'wall_mid': 'Single wall in middle position'
        }
        
        ttk.Label(frame, text="Description:").grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self.desc_text = tk.Text(frame, height=4, width=50)
        self.desc_text.grid(row=1, column=1, sticky='w', padx=5, pady=5)
        
        # Update description when experiment type changes
        def update_description(*args):
            exp_type = self.experiment_var.get()
            self.desc_text.delete(1.0, tk.END)
            self.desc_text.insert(1.0, descriptions.get(exp_type, ''))
        
        self.experiment_var.trace_add('write', lambda *args: update_description())
        update_description()  # Initial description
        
    def create_wall_tab(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Wall Configuration")
        
        # Number of walls
        ttk.Label(frame, text="Number of Walls (1-5):").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.num_walls_var = tk.IntVar(value=1)
        num_walls_spin = ttk.Spinbox(frame, from_=1, to=5, textvariable=self.num_walls_var, width=10)
        num_walls_spin.grid(row=0, column=1, sticky='w', padx=5, pady=5)
        
        # Wall positions and masses
        ttk.Label(frame, text="Wall Positions (left to right):").grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self.positions_text = tk.Text(frame, height=3, width=40)
        self.positions_text.grid(row=1, column=1, sticky='w', padx=5, pady=5)
        self.positions_text.insert(1.0, "0.0")
        
        ttk.Label(frame, text="Wall Mass Factors (left to right):").grid(row=2, column=0, sticky='nw', padx=5, pady=5)
        self.masses_text = tk.Text(frame, height=3, width=40)
        self.masses_text.grid(row=2, column=1, sticky='w', padx=5, pady=5)
        self.masses_text.insert(1.0, "200.0")
        
        # Auto-generate positions button
        def auto_generate():
            num_walls = self.num_walls_var.get()
            positions = []
            masses = []
            for i in range(num_walls):
                positions.append(f"{i * 5.0:.1f}")  # 5 sigma units apart
                masses.append("200.0")
            
            self.positions_text.delete(1.0, tk.END)
            self.positions_text.insert(1.0, ", ".join(positions))
            self.masses_text.delete(1.0, tk.END)
            self.masses_text.insert(1.0, ", ".join(masses))
        
        ttk.Button(frame, text="Auto-generate", command=auto_generate).grid(row=3, column=1, sticky='w', padx=5, pady=5)
        
    def create_protocol_tab(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Piston Protocol")
        
        # Protocol type
        ttk.Label(frame, text="Protocol Type:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        self.protocol_var = tk.StringVar(value='sigmoidal')
        protocol_combo = ttk.Combobox(frame, textvariable=self.protocol_var, width=20)
        protocol_combo['values'] = ['step', 'sigmoidal', 'linear', 'sinusoidal', 'optimal']
        protocol_combo.grid(row=0, column=1, sticky='w', padx=5, pady=5)
        
        # Protocol parameters
        ttk.Label(frame, text="Max Velocity:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        self.max_velocity_var = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self.max_velocity_var, width=10).grid(row=1, column=1, sticky='w', padx=5, pady=5)
        
        ttk.Label(frame, text="Protocol Duration:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        self.duration_var = tk.DoubleVar(value=10.0)
        ttk.Entry(frame, textvariable=self.duration_var, width=10).grid(row=2, column=1, sticky='w', padx=5, pady=5)
        
        # Protocol descriptions
        protocol_desc = {
            'step': 'Instant acceleration (unphysical)',
            'sigmoidal': 'Smooth sigmoidal acceleration (recommended)',
            'linear': 'Linear acceleration ramp',
            'sinusoidal': 'Sinusoidal velocity profile',
            'optimal': 'AI-optimized protocol (experimental)'
        }
        
        ttk.Label(frame, text="Description:").grid(row=3, column=0, sticky='nw', padx=5, pady=5)
        self.protocol_desc_text = tk.Text(frame, height=3, width=50)
        self.protocol_desc_text.grid(row=3, column=1, sticky='w', padx=5, pady=5)
        
        def update_protocol_desc(*args):
            protocol = self.protocol_var.get()
            self.protocol_desc_text.delete(1.0, tk.END)
            self.protocol_desc_text.insert(1.0, protocol_desc.get(protocol, ''))
        
        self.protocol_var.trace_add('write', lambda *args: update_protocol_desc())
        update_protocol_desc()
        
    def create_simulation_tab(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Simulation Parameters")
        
        # Basic parameters
        params = [
            ('Particles:', 'particles', 1000, 'Total number of particles'),
            ('L0 Units:', 'l0_units', 20.0, 'Box half-length in sigma units'),
            ('Height Units:', 'height_units', 10.0, 'Box height in sigma units'),
            ('Wall Thickness:', 'wall_thickness', 1.0, 'Wall thickness in sigma units'),
            ('Repeats:', 'repeats', 1, 'Number of experiment repeats'),
            ('Steps:', 'steps', 3000, 'Simulation steps')
        ]
        
        self.param_vars = {}
        for i, (label, key, default, tooltip) in enumerate(params):
            ttk.Label(frame, text=label).grid(row=i, column=0, sticky='w', padx=5, pady=5)
            var = tk.DoubleVar(value=default) if isinstance(default, float) else tk.IntVar(value=default)
            self.param_vars[key] = var
            entry = ttk.Entry(frame, textvariable=var, width=15)
            entry.grid(row=i, column=1, sticky='w', padx=5, pady=5)
            
            # Add tooltip
            ttk.Label(frame, text=f"({tooltip})", font=('Arial', 8)).grid(row=i, column=2, sticky='w', padx=5, pady=5)
        
    def create_control_buttons(self):
        button_frame = ttk.Frame(self.root)
        button_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Button(button_frame, text="Generate Command", command=self.generate_command).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Save Configuration", command=self.save_config).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Load Configuration", command=self.load_config).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Run Simulation", command=self.run_simulation).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Exit", command=self.root.quit).pack(side='right', padx=5)
        
        # Command output
        ttk.Label(button_frame, text="Generated Command:").pack(side='left', padx=5)
        self.command_text = tk.Text(button_frame, height=2, width=60)
        self.command_text.pack(side='left', padx=5)
        
    def generate_command(self):
        """Generate the CLI command based on current configuration"""
        cmd_parts = ["./00ALLINONE"]
        
        # Experiment type
        cmd_parts.append(f"--experiment={self.experiment_var.get()}")
        
        # Wall configuration
        num_walls = self.num_walls_var.get()
        cmd_parts.append(f"--num-walls={num_walls}")
        
        # Wall positions
        positions = self.positions_text.get(1.0, tk.END).strip()
        if positions:
            cmd_parts.append(f"--wall-positions={positions}")
        
        # Wall masses
        masses = self.masses_text.get(1.0, tk.END).strip()
        if masses:
            cmd_parts.append(f"--wall-mass-factors={masses}")
        
        # Protocol
        cmd_parts.append(f"--protocol={self.protocol_var.get()}")
        
        # Simulation parameters
        for key, var in self.param_vars.items():
            if key == 'particles':
                cmd_parts.append(f"--particles={int(var.get())}")
            elif key == 'l0_units':
                cmd_parts.append(f"--l0={var.get()}")
            elif key == 'height_units':
                cmd_parts.append(f"--height={var.get()}")
            elif key == 'wall_thickness':
                cmd_parts.append(f"--wall-thickness={var.get()}")
            elif key == 'repeats':
                cmd_parts.append(f"--repeats={int(var.get())}")
            elif key == 'steps':
                cmd_parts.append(f"--steps={int(var.get())}")
        
        command = " ".join(cmd_parts)
        self.command_text.delete(1.0, tk.END)
        self.command_text.insert(1.0, command)
        
        return command
        
    def save_config(self):
        """Save current configuration to file"""
        config = {
            'experiment_type': self.experiment_var.get(),
            'num_walls': self.num_walls_var.get(),
            'wall_positions': self.positions_text.get(1.0, tk.END).strip(),
            'wall_masses': self.masses_text.get(1.0, tk.END).strip(),
            'protocol': self.protocol_var.get(),
            'max_velocity': self.max_velocity_var.get(),
            'duration': self.duration_var.get(),
            'parameters': {key: var.get() for key, var in self.param_vars.items()}
        }
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            with open(filename, 'w') as f:
                json.dump(config, f, indent=2)
            messagebox.showinfo("Success", f"Configuration saved to {filename}")
    
    def load_config(self):
        """Load configuration from file"""
        filename = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'r') as f:
                    config = json.load(f)
                
                # Update GUI with loaded config
                self.experiment_var.set(config.get('experiment_type', 'speed_of_sound'))
                self.num_walls_var.set(config.get('num_walls', 1))
                self.positions_text.delete(1.0, tk.END)
                self.positions_text.insert(1.0, config.get('wall_positions', '0.0'))
                self.masses_text.delete(1.0, tk.END)
                self.masses_text.insert(1.0, config.get('wall_masses', '200.0'))
                self.protocol_var.set(config.get('protocol', 'sigmoidal'))
                
                # Update parameters
                params = config.get('parameters', {})
                for key, var in self.param_vars.items():
                    if key in params:
                        var.set(params[key])
                
                messagebox.showinfo("Success", f"Configuration loaded from {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load configuration: {e}")
    
    def run_simulation(self):
        """Run the simulation with current configuration"""
        command = self.generate_command()
        
        try:
            # Change to the correct directory
            os.chdir('/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3')
            
            # Run the simulation
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                messagebox.showinfo("Success", "Simulation completed successfully!")
            else:
                messagebox.showerror("Error", f"Simulation failed:\n{result.stderr}")
                
        except Exception as e:
            messagebox.showerror("Error", f"Failed to run simulation: {e}")

def main():
    root = tk.Tk()
    app = ExperimentConfigGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
