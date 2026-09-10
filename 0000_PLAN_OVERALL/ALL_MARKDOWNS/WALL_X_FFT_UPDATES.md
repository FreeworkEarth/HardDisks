# wall_x_FFT.py Updates

## Changes Made

Updated the analysis script to work with new timestamped output folders.

### Before (Broken):
```python
folder_path = "experiments_speed_of_sound/mode1_normalized_units/mode1_normalized_units_1particle_1_pixel_1_1_factor/"
```

This hardcoded path no longer exists!

### After (Fixed):
```python
# Auto-detect latest simulation folder
base_path = "experiments_speed_of_sound/mode1_normalized_units/"
folder_path = find_latest_simulation_folder(base_path, use_specific_timestamp)
```

## New Features

### 1. Auto-Detect Latest Run (Default)

Set `use_specific_timestamp = None` at the top of the script:

```python
use_specific_timestamp = None  # Auto-detect latest
```

The script will automatically find and analyze the most recent `simulation_*` folder.

**Example output:**
```
📁 Auto-detected latest run: simulation_05_11_24_14_30
```

### 2. Analyze Specific Run (Optional)

To analyze a specific timestamped run, set:

```python
use_specific_timestamp = "simulation_05_11_24_14_30"
```

**Example output:**
```
📁 Using specified run: simulation_05_11_24_14_30
```

### 3. Backwards Compatibility

If no timestamped folders exist, the script falls back to the base directory:

**Example output:**
```
ℹ️  No timestamped folders found, using base path: experiments_speed_of_sound/mode1_normalized_units/
```

## Usage

### Run New Experiments

```bash
cd hspist3
./00ALLINONE --mode=experiments --particles=100 --temperature=1000 --kbt1 \
  --experiment-lengths 7.5,10,15,20,25,30,35 \
  --experiment-wall-masses 200 \
  --experiment-repeats 3
```

This creates:
```
experiments_speed_of_sound/mode1_normalized_units/simulation_05_11_24_14_30/
  ├── wall_x_positions_L0_75_wallmassfactor_200_run0.csv
  ├── wall_x_positions_L0_75_wallmassfactor_200_run1.csv
  ├── ...
```

### Analyze Latest Run

```bash
cd hspist3
source venv/bin/activate
python3 wall_x_FFT.py
```

Output:
```
📁 Auto-detected latest run: simulation_05_11_24_14_30
Processing files from: .../simulation_05_11_24_14_30/
...
```

### Analyze Specific Run

Edit `wall_x_FFT.py`:
```python
use_specific_timestamp = "simulation_05_11_24_12_30"  # Morning run
```

Then run:
```bash
python3 wall_x_FFT.py
```

Output:
```
📁 Using specified run: simulation_05_11_24_12_30
Processing files from: .../simulation_05_11_24_12_30/
...
```

## List Available Runs

```bash
ls -lt experiments_speed_of_sound/mode1_normalized_units/
```

Output example:
```
drwxr-xr-x  simulation_05_11_24_16_45/  # Latest (4:45 PM)
drwxr-xr-x  simulation_05_11_24_14_30/  # Earlier (2:30 PM)
drwxr-xr-x  simulation_05_11_24_12_15/  # Morning (12:15 PM)
```

## Important Notes

### Your Data Is Now Safe!

The old `rm -f *.csv` commands have been removed. Each experiment run gets its own timestamped folder, so you can:
- Compare results across different runs
- Go back to previous data anytime
- Keep a full history of experiments

### Disk Space

Old runs are NOT automatically deleted. If you run many experiments, clean up manually:

```bash
# List all runs with sizes
du -sh experiments_speed_of_sound/mode1_normalized_units/simulation_*/

# Delete old runs you don't need
rm -rf experiments_speed_of_sound/mode1_normalized_units/simulation_04_11_24_*
```

### Timestamp Format

- Format: `simulation_DD_MM_YY_HH_MM`
- Example: `simulation_05_11_24_14_30`
  - 05: Day (5th)
  - 11: Month (November)
  - 24: Year (2024)
  - 14: Hour (2 PM)
  - 30: Minute (30)

## Testing

Test the auto-detection:

```bash
cd hspist3
python3 -c "
import os
import glob

base_path = 'experiments_speed_of_sound/mode1_normalized_units/'
pattern = os.path.join(base_path, 'simulation_*')
folders = glob.glob(pattern)

if folders:
    latest = max(folders, key=os.path.getmtime)
    print(f'Latest run: {os.path.basename(latest)}')
    print(f'Full path: {latest}')
    print(f'CSV files: {len(glob.glob(os.path.join(latest, \"*.csv\")))}')
else:
    print('No timestamped folders found!')
"
```

Expected output:
```
Latest run: simulation_05_11_24_14_30
Full path: experiments_speed_of_sound/mode1_normalized_units/simulation_05_11_24_14_30
CSV files: 21
```

