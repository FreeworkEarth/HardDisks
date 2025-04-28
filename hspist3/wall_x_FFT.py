import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.markers as mmarkers
from scipy.signal import spectrogram
from scipy.ndimage import gaussian_filter1d
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.signal import welch
from scipy.optimize import root_scalar
from matplotlib.cm import viridis
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
# Simulated structure — replace with your real aggregated data
from collections import defaultdict
stats = defaultdict(list)  # key: (L0, M), value: list of peak frequencies (or amps)
# Then, once you gather all values into stats[(L0, M)], compute:
from statistics import mean, stdev
import re
import os
import glob

## source venv/bin/activate
### --- Configurable flags ---
enable_batch_processing =False
     # Set to False to analyze one file only
show_plots = False                  # Set to False to disable plot display
FFT_mode = "windowed_hanning_SCIPY"       # Options: "raw", "centered", "windowed_hanning","windowed_hanning_SCIPY", "welchs_smoothing_on_windowed_signal_hanning"
# --- Choose which peak you want to return ---

peak_selection_mode = "raw_near_smoothed_lowfreq"  # Options: 'smoothed_max', 'raw_max', 'raw_near_smoothed_lowfreq'← your default

# Constants from simulation in C
radius = 1
N = 50  # per side
A = 10
K = np.pi/2
pi = np.pi
k_B = 1.380649e-23  # Boltzmann constant
T = 300  # Temperature in Kelvin
mass_particle = 1.67e-27  # Example mass (adjust as appropriate)
n_value = 1  # Assuming fundamental mode for closed tube



print(os.getcwd())  # this is where the script expects to find the folder
# Check what files exist in the expected folder path

folder_path = os.path.abspath("experiments_speed_of_sound/mode0_real_units/")
single_plot_filename = "wall_x_positions_L0_75_wallmassfactor_100_run0.csv"  # Leave as None to plot all

assert os.path.isdir(folder_path), f"❌ Folder path does not exist: {folder_path}"
plot_title = "Peak Frequency Shift vs Wall Mass and Box Length" 

# Show what folder we are looking in
print("🔍 Current working directory:", os.getcwd())
print("🔍 Full search path:", os.path.join(folder_path, "*.csv"))

# Debug what files it actually sees
filepaths = sorted(glob.glob(os.path.join(folder_path, "*.csv")))
print(f"📂 Files found: {filepaths}")

### --- Storage for results ---
results = []






################################################## SIGMA ESTIMATOR (dominant frequency from raw wall signal)
def estimate_sigma_from_dominant_freq(signal, dt, min_sigma=1.0, max_sigma=10.0, freq_floor=0.01, freq_ceiling=2.0):
    """
    Estimate smoothing sigma based on dominant frequency of signal.
    Faster oscillation → higher frequency → higher sigma.
    Slower oscillation → lower frequency → lower sigma.
   
    Estimate smoothing sigma based on dominant frequency of signal.
    Returns both sigma and dominant frequency.
    """
    # Center signal
    signal_centered = signal - np.mean(signal)

    ############### FFT

    # ON CENTERED SIGNAL
    #fft_result = np.fft.fft(signal_centered)

    # ON WINDOWED SIGNAL
    window = np.hanning(len(signal))  # Or any other window
    windowed_signal = signal * window
    fft_result = np.fft.fft(windowed_signal)


    #freqs = np.fft.fftfreq(len(signal_centered), d=dt)
    freqs = np.fft.fftfreq(len(signal_centered), d=dt)

    power = np.abs(fft_result)**2

    # Use only positive frequencies
    pos_mask = freqs > 0
    freqs = freqs[pos_mask]
    power = power[pos_mask]

    # Dominant frequency
    dominant_freq = freqs[np.argmax(power)]

    # Normalize and map to sigma
    norm_freq = np.clip((dominant_freq - freq_floor) / (freq_ceiling - freq_floor), 0, 1)
    sigma = min_sigma + (max_sigma - min_sigma) * norm_freq

    print(f"🎯 Dominant frequency of wall signal: {dominant_freq:.4f} Hz")
    print(f"📏 Estimated smoothing sigma: {sigma:.2f}")

    return sigma, dominant_freq
    

def estimate_sigma_from_spectral_centroid(signal, dt, min_sigma=1.0, max_sigma=10.0, freq_floor=0.01, freq_ceiling=2.0):
    """
    Estimate smoothing sigma using spectral centroid (frequency center of mass).
    This is more robust against noise and wild jitter.

    What is the spectral centroid?

        The spectral centroid is a measure of the “center of mass” of a signal’s power spectrum. In simple terms, it tells you where the bulk of the frequency energy is located. 
        It’s computed as:

        \text{Centroid} = \frac{\sum f_i \cdot A(f_i)}{\sum A(f_i)}

        Where:
            •	f_i = frequency bin
            •	A(f_i) = amplitude at that frequency


    That’s the spectral centroid as power of A.
        Mathematically:
        f_{\text{centroid}} = \frac{\sum f_i \cdot P(f_i)}{\sum P(f_i)}

        Where:
            •	f_i = frequency bin
            •	P(f_i) = power (squared amplitude) at that frequency

    In audio and music:
    •	High centroid = bright or noisy sound (e.g., cymbals, white noise)
    •	Low centroid = dull, bassy sound (e.g., cello, drum kick)

    In speed of sound case:
    •	Fast, jittery wall = energy spread over many higher freqs → high centroid → high sigma
    •	Heavy, sluggish wall = energy bunched in low freqs → low centroid → low sigma

    Why it works better than just taking max freq?

        Because the maximum frequency power is often unreliable:
            •	Jitter/noise can spike individual bins
            •	Peaks can bounce around from run to run
            •	Doesn’t tell you if most energy is low or high

        Centroid fixes that:
            •	It integrates over the whole FFT curve
            •	Smooths out anomalies
            •	Feels the true average behavior of the signal

    Why it’s perfect for your smoothing sigma?

        Smoothing = blur the signal so we only care about long-term trends.
        But how much to blur depends on how “twitchy” the signal is.

        So we:
            1.	Compute the centroid
            2.	Normalize it from “slow” to “fast”
            3.	Map that to a sigma range: slow = low sigma, fast = high sigma

    """
    signal_centered = signal - np.mean(signal)
    fft_result = np.fft.fft(signal_centered)
    freqs = np.fft.fftfreq(len(signal_centered), d=dt)
    power = np.abs(fft_result)**2

    # Positive frequencies only
    pos_mask = freqs > 0
    freqs = freqs[pos_mask]
    power = power[pos_mask]

    # Spectral centroid: weighted average of frequencies
    spectral_centroid = np.sum(freqs * power) / np.sum(power)

    # Normalize and map to sigma
    norm_freq = np.clip((spectral_centroid - freq_floor) / (freq_ceiling - freq_floor), 0, 1)
    sigma = min_sigma + (max_sigma - min_sigma) * norm_freq

    print(f"🎯 Spectral centroid of wall signal: {spectral_centroid:.4f} Hz")
    print(f"📏 Estimated smoothing sigma: {sigma:.2f}")

    return sigma, spectral_centroid



def annotate_peaks_smart1(freqs, amplitudes, ax, label_prefix="Peak",
                         yscale='linear', min_freq=0.0, smooth_sigma=2,
                         prominence=None, color='red', max_peaks=3,
                         show_smoothed_curve=True, smoothed_color='red'):
    """
    Detects top N peaks using smoothed signal, then refines location from raw signal.
    Shows:
    - Smoothed peak
    - Raw value at smoothed peak
    - Refined raw peak in a window
    - Global raw peak
    """

    freqs_full = freqs
    amplitudes_full = amplitudes
    smoothed_full = gaussian_filter1d(amplitudes_full, sigma=smooth_sigma)

    # ⚠️ Include 0 Hz exactly
    valid_idx = freqs_full >= min_freq
    freqs = freqs_full[valid_idx]
    amplitudes = amplitudes_full[valid_idx]
    smoothed_amp = smoothed_full[valid_idx]

    print(f"🔬 Searching in frequency range {freqs[0]:.4f} to {freqs[-1]:.4f}")
    print(f"🔬 Max smoothed amplitude: {np.max(smoothed_amp):.2e}")
    print(f"🔬 Median: {np.median(smoothed_amp):.2e}, Std: {np.std(smoothed_amp):.2e}")

    # === Plot smoothed curve ===
    if show_smoothed_curve:
        ax.plot(freqs, smoothed_amp, label=f"Smoothed FFT (σ={smooth_sigma})",
                color=smoothed_color, linestyle='--', alpha=0.6)

    # === Adaptive prominence or override ===
    if prominence is None:
        max_amp = np.max(smoothed_amp)
        median_amp = np.median(smoothed_amp)
        std_amp = np.std(smoothed_amp)
        base_prominence = max(1e-12, 0.01 * max_amp, median_amp + 0.5 * std_amp)
    else:
        base_prominence = prominence

    print(f"🔍 Prominence base: {base_prominence:.2e}")

    # === Peak search ===
    for scale in [1.0, 1.5, 2.0, 3.0, 5.0]:
        prom = base_prominence * scale
        peaks, properties = find_peaks(smoothed_amp, prominence=prom)
        if 0 < len(peaks) <= max_peaks:
            break
    else:
        print("⚠️ No peaks found.")
        print("⚠️ Smart prominence-based peak finder failed. Using basic fallback peak detection.")
        
        # Fallback: naive peak finder on smoothed curve with no prominence, just local maxima
        basic_peaks, _ = find_peaks(smoothed_amp)
        
        if len(basic_peaks) == 0:
            print("❌ Even fallback peak finder found no peaks.")
            return None
        
        # Pick the peak with max amplitude in the fallback result
        fallback_idx = basic_peaks[np.argmax(smoothed_amp[basic_peaks])]
        fallback_freq = freqs[fallback_idx]
        fallback_amp = amplitudes_full[np.argmin(np.abs(freqs_full - fallback_freq))]

        ax.plot(fallback_freq, fallback_amp, 'o', color='cyan', label=f"Fallback Peak @ {fallback_freq:.2f} Hz")
        ax.annotate(f"{fallback_freq:.2f} Hz\nFallback", 
                    xy=(fallback_freq, fallback_amp), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8, color='cyan')

        print(f"📍 Fallback peak: {fallback_freq:.3f} Hz, Raw Amp: {fallback_amp:.2e}")
        
        return fallback_freq, fallback_amp
        

    print(f"✅ Found {len(peaks)} peaks (prominence: {prom:.2e})")

    prominences = properties["prominences"]
    top_peaks = np.argsort(prominences)[::-1][:max_peaks]
    selected_peaks = peaks[top_peaks]

    highest_peak = None
    highest_amp = -np.inf

    for i, idx in enumerate(selected_peaks):
        freq = freqs[idx]
        amp_smoothed = smoothed_amp[idx]

        # Refined peak search in raw signal near smoothed peak
        bandwidth = 0.25  # Hz window
        window_mask = (freqs_full >= freq - bandwidth) & (freqs_full <= freq + bandwidth)
        if np.any(window_mask):
            window_freqs = freqs_full[window_mask]
            window_amps = amplitudes_full[window_mask]
            refined_idx = np.argmax(window_amps)
            raw_peak_freq = window_freqs[refined_idx]
            raw_peak_amp = window_amps[refined_idx]
        else:
            raw_peak_freq = freq
            raw_peak_amp = amplitudes_full[np.argmin(np.abs(freqs_full - freq))]

        # Plot peaks
        ax.plot(freq, amp_smoothed, 'o', color='red', label=f"Smoothed Peak @ {freq:.2f} Hz")
        ax.plot(raw_peak_freq, raw_peak_amp, 'o', color='blue', label=f"Refined Raw Peak @ {raw_peak_freq:.2f} Hz")
        ax.axvline(freq, linestyle='--', color='gray', alpha=0.5)

        ax.annotate(f"{freq:.2f} Hz\nSmoothed: {amp_smoothed:.2e}",
                    xy=(freq, amp_smoothed), xytext=(5, 10),
                    textcoords='offset points', fontsize=8, color='red')

        ax.annotate(f"{raw_peak_freq:.2f} Hz\nRaw: {raw_peak_amp:.2e}",
                    xy=(raw_peak_freq, raw_peak_amp), xytext=(5, -15),
                    textcoords='offset points', fontsize=8, color='blue')

        print(f"📍 Peak {i+1} — Smoothed @ {freq:.3f} Hz, Raw refined @ {raw_peak_freq:.3f} Hz")

        if raw_peak_amp > highest_amp:
            highest_amp = raw_peak_amp
            highest_peak = (raw_peak_freq, raw_peak_amp)

    # Global raw max
    raw_peak_idx = np.argmax(amplitudes_full)
    raw_peak_freq = freqs_full[raw_peak_idx]
    raw_peak_amp = amplitudes_full[raw_peak_idx]
    ax.plot(raw_peak_freq, raw_peak_amp, 'X', color='blue', markersize=8, label=f"Raw Highest Peak @ {raw_peak_freq:.2f} Hz")
    ax.annotate(f"{raw_peak_freq:.2f} Hz\n{raw_peak_amp:.2e}",
                xy=(raw_peak_freq, raw_peak_amp),
                xytext=(5, -10), textcoords='offset points',
                fontsize=8, color='blue')

    if yscale == 'log':
        ax.set_yscale('log')

    # Deduplicate legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=9)

    if highest_peak:
        print(f"📍 ANNOTATE_PEAKS_SMART - Highest Smoothed Peak @ {highest_peak[0]:.3f} Hz — Raw: {highest_peak[1]:.2e}")
        return highest_peak
    else:
        print("⚠️ INSIDE ANNOTATE_PEAKS_SMART - No highest peak found.")
        return None



def annotate_peaks_smart(freqs, amplitudes, ax, label_prefix="Peak",
                         yscale='linear', min_freq=0.0, max_freq=None,
                         smooth_sigma=2, prominence=None, color='red',
                         max_peaks=3, show_smoothed_curve=True,
                         smoothed_color='red', return_mode="refined",
                         fallback=True, bandwidth=0.25):
    """
    Detects and plots peaks using smoothed spectrum.
    Also refines peak position from raw signal nearby.
    return_mode: 'smoothed', 'refined', 'raw'
    fallback: enable if no peak found
    """

    freqs_full = np.asarray(freqs)
    amplitudes_full = np.asarray(amplitudes)

    if max_freq is None:
        max_freq = np.max(freqs_full)
    valid_idx = (freqs_full >= min_freq) & (freqs_full <= max_freq)
    freqs = freqs_full[valid_idx]
    amplitudes = amplitudes_full[valid_idx]

    smoothed_amp = gaussian_filter1d(amplitudes, sigma=smooth_sigma)

    print(f"🔬 Searching from {freqs[0]:.4f} Hz to {freqs[-1]:.4f} Hz")
    print(f"🔬 Max smoothed: {np.max(smoothed_amp):.2e}, Median: {np.median(smoothed_amp):.2e}, Std: {np.std(smoothed_amp):.2e}")

    if show_smoothed_curve:
        ax.plot(freqs, smoothed_amp, label=f"Smoothed FFT (σ={smooth_sigma})",
                color=smoothed_color, linestyle='--', alpha=0.6)

    if prominence is None:
        base_prominence = max(
            1e-12,
            0.01 * np.max(smoothed_amp),
            np.median(smoothed_amp) + 0.5 * np.std(smoothed_amp)
        )
    else:
        base_prominence = prominence

    print(f"🔍 Prominence base: {base_prominence:.2e}")

    peaks, properties = find_peaks(smoothed_amp, prominence=base_prominence)
    if len(peaks) == 0:
        print("⚠️ No peaks found on smoothed signal.")
        if fallback:
            fallback_peaks, _ = find_peaks(smoothed_amp)
            if len(fallback_peaks) > 0:
                idx = fallback_peaks[np.argmax(smoothed_amp[fallback_peaks])]
                fallback_freq = freqs[idx]
                raw_idx = np.argmin(np.abs(freqs_full - fallback_freq))
                fallback_amp = amplitudes_full[raw_idx]
                ax.plot(fallback_freq, fallback_amp, 'o', color='cyan', label=f"Fallback Peak @ {fallback_freq:.2f} Hz")
                ax.annotate(f"{fallback_freq:.2f} Hz\nFallback", (fallback_freq, fallback_amp),
                            textcoords='offset points', xytext=(5, 5), fontsize=8, color='cyan')
                return fallback_freq, fallback_amp
            else:
                print("❌ No fallback peaks found either.")
                return None
        else:
            return None

    print(f"✅ Found {len(peaks)} peaks")
    top_peaks = np.argsort(properties["prominences"])[::-1][:max_peaks]
    selected_peaks = peaks[top_peaks]

    # Refine from raw
    refined_peaks = []
    for i, idx in enumerate(selected_peaks):
        f_smooth = freqs[idx]
        a_smooth = smoothed_amp[idx]

        # refine raw peak
        mask = (freqs_full >= f_smooth - bandwidth) & (freqs_full <= f_smooth + bandwidth)
        if np.any(mask):
            f_local = freqs_full[mask]
            a_local = amplitudes_full[mask]
            best_idx = np.argmax(a_local)
            f_refined = f_local[best_idx]
            a_refined = a_local[best_idx]
        else:
            f_refined = f_smooth
            a_refined = amplitudes_full[np.argmin(np.abs(freqs_full - f_smooth))]

        # Plot both
        ax.plot(f_smooth, a_smooth, 'o', color='red', label=f"Smoothed Peak @ {f_smooth:.2f} Hz")
        ax.plot(f_refined, a_refined, 'o', color='blue', label=f"Refined Raw Peak @ {f_refined:.2f} Hz")
        ax.axvline(f_smooth, linestyle='--', color='gray', alpha=0.4)

        ax.annotate(f"{f_smooth:.2f} Hz\nS: {a_smooth:.2e}", (f_smooth, a_smooth),
                    textcoords='offset points', xytext=(5, 10), fontsize=8, color='red')
        ax.annotate(f"{f_refined:.2f} Hz\nR: {a_refined:.2e}", (f_refined, a_refined),
                    textcoords='offset points', xytext=(5, -15), fontsize=8, color='blue')

        refined_peaks.append((f_refined, a_refined))

    # Global raw peak
    global_idx = np.argmax(amplitudes_full)
    global_freq = freqs_full[global_idx]
    global_amp = amplitudes_full[global_idx]
    ax.plot(global_freq, global_amp, 'X', color='navy', markersize=8, label=f"Raw Global Peak @ {global_freq:.2f} Hz")
    ax.annotate(f"{global_freq:.2f} Hz\n{global_amp:.2e}",
                (global_freq, global_amp),
                textcoords='offset points', xytext=(5, -10), fontsize=8, color='navy')

    if yscale == 'log':
        ax.set_yscale('log')

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=9)



    # --- Decide what to return ---
    if return_mode == "smoothed":
        return freqs[selected_peaks[0]], smoothed_amp[selected_peaks[0]]
    elif return_mode == "refined":
        return refined_peaks[0]
    elif return_mode == "raw":
        return global_freq, global_amp
    else:
        print(f"⚠️ Unknown return_mode '{return_mode}', defaulting to refined.")
        return refined_peaks[0]


######### FFT Theory
'''
Why and what is the FFT = Fast Fourier Transform?
We want to understand what frequencies are present in a time-domain signal (like oscillating wall position).
It’s an algorithm that efficiently computes the Discrete Fourier Transform (DFT) of a signal.

Step-by-step: What the FFT is doing

Let’s say you have a discrete time signal It could be Wall_X, sampled every dt seconds.
: x = [x_0, x_1, x_2, ..., x_{N-1}]

1. Assume periodicity
FFT assumes that your signal repeats itself every N samples — which is why sharp edges cause problems (discontinuity leads to noise — → spectral leakage).

This is called the “periodic extension” of the signal.
2. Decompose into sinusoids
The DFT (and FFT) says:
Any signal can be represented as a sum of sines and cosines of various frequencies, amplitudes, and phases.

Mathematically:
X_k = \sum_{n=0}^{N-1} x_n \cdot e^{-j 2\pi k n / N}

Where:
	•	X_k = frequency component at frequency bin k
	•	x_n = time-domain value at time sample n
	•	e^{-j 2πkn/N} = a complex sinusoid (Euler’s formula)

You get:
	•	A list of complex numbers X = [X_0, X_1, ..., X_{N-1}]
	•	Each complex value tells you: how much of that frequency is present (and with what phase).



What do the results mean?

1. Frequency bins
freqs = np.fft.fftfreq(N, d=dt)
These give you the frequency in Hz of each bin k.

2. Amplitude
amplitude = np.abs(X) / N
Gives you the magnitude of the frequency component — how strong that frequency is in your signal.

3. Power Spectrum square of the amplitude:
power = amplitude**2

Why use power?
	•	Physically, power ∝ energy in that frequency.
	•	Power doesn’t care about phase.
	•	Peaks in power spectrum often stand out more clearly (it exaggerates strong components and suppresses small noise).
	•	In physics and engineering, power spectrum = how energy is distributed over frequency.

FFT
    Fast algorithm to compute frequency components of a signal
    X_k
    Complex number: how much of that frequency is present
    np.abs(X_k)
    Amplitude at frequency f_k
    np.abs(X_k)**2
    Power at frequency f_k
    np.angle(X_k)
    Phase of that frequency component
    np.fft.fftfreq()
    Gives you the real Hz associated with each index k

 Important Notes
	•	FFT gives both positive and negative frequencies (use freqs >= 0 mask).
	•	The DC component (0 Hz) is the mean — you often subtract it before FFT.
	•	Windowing reduces “leakage” from discontinuities at the signal edges.

⸻

'''




def extract_params_from_filename(filename):
    match = re.search(r'L0_(\d+)_wallmassfactor_(\d+)_run(\d+)', filename)
    if match:
        L0 = int(match.group(1)) / 10.0  # Turn 100 into 10.0
        wall_mass_factor = int(match.group(2))
        run_number = int(match.group(3))
        return L0, wall_mass_factor, run_number
    else:
        return None, None, None



#### SOLVE TRANSCENDENTAL EQUATION
def transcendental_eq(x, wall_mass_factor, N):
    return 1 / np.tan(x) - (wall_mass_factor / (2 * N)) * x


def find_transcendental_root(wall_mass_factor, N, x_guess_min=0.01, x_guess_max=5):
    def eq(x):
        return 1/np.tan(x) - (wall_mass_factor / (2 * N)) * x
    sol = root_scalar(eq, bracket=[x_guess_min, x_guess_max], method='brentq')
    if sol.converged:
        return sol.root
    else:
        return None

# ✅ Plotting function
def plot_transcendental(wall_mass_factor, N, x_min=0.01, x_max=5):
    x = np.linspace(x_min, x_max, 1000)
    cot_x = 1 / np.tan(x)
    line = (wall_mass_factor / (2 * N)) * x

    plt.figure(figsize=(10, 6))
    plt.plot(x, cot_x, label=r"$\cot(x)$", color='royalblue')
    plt.plot(x, line, label=fr"$\frac{{\text{{wall_mass_factor}}}}{{2N}} x$" + f" (wall_mass_factor={wall_mass_factor})", color='darkorange')
    plt.ylim(-5, 5)
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel("x")
    plt.ylabel("Value")
    plt.title("Transcendental Equation: Intersection gives the root")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Example usage:
wall_mass_factor = 20
N = 50

# Solve
K_root = find_transcendental_root(wall_mass_factor, N)
print(f"✅ Found root K: {K_root:.6f}")

# Plot
#plot_transcendental(wall_mass_factor, N)


# --- Process a single file ---
def process_file(filepath):
    df = pd.read_csv(filepath)
    df.columns = df.columns.str.strip()

    if 'Time' not in df.columns or 'Wall_X' not in df.columns:
        print(f"❌ Skipping {filepath} — missing columns")
        return None

    time = df['Time'].values
    wall_x = df['Wall_X'].values
    dt = np.mean(np.diff(time))

    # center wall position to remove DC offset
    wall_x_centered = wall_x - np.mean(wall_x)

    # --- Plot: Wall X Position ---
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(time, wall_x, label="Wall X Position", color='royalblue')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [units]")
    ax.set_title("Wall X Position Over Time")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()

    if enable_batch_processing and not show_plots:
        plt.close(fig)  # Don’t show intermediate plots
    else:
        plt.show()


    ###### FFT ANALYSIS ######
    if FFT_mode == "raw":
        # ON RAW SIGNAL
        fft_result = np.fft.fft(wall_x)
        freqs = np.fft.fftfreq(len(wall_x), d=dt)
        amplitude = np.abs(fft_result) / len(wall_x)
        
    elif FFT_mode == "centered":
        # --- FFT: Wall X Position CENTERED ---
        fft_result = np.fft.fft(wall_x_centered)
        freqs = np.fft.fftfreq(len(wall_x_centered), d=dt)
        amplitude = np.abs(fft_result) / len(wall_x_centered)

    elif FFT_mode == "windowed_hanning":
        # ON WINDOWED SIGNAL
        window = np.hanning(len(wall_x))  # Or any other window
        windowed_signal = wall_x * window
        fft_result = np.fft.fft(windowed_signal)
        freqs = np.fft.fftfreq(len(wall_x), d=dt)
        amplitude = np.abs(fft_result) / len(wall_x)

    elif FFT_mode == "windowed_hanning_SCIPY":
        window = np.hanning(len(wall_x_centered))
        windowed_signal = wall_x_centered * window
        fft_result = fft(windowed_signal)
        freqs = fftfreq(len(windowed_signal), d=dt)
        amplitude = np.abs(fft_result) / len(windowed_signal)
        power = amplitude**2
    
    elif FFT_mode == "welchs_smoothing_on_windowed_signal_hanning":
        ####### WELCHS smoothing 
        freqs, power = welch(wall_x, fs=1/dt, window='hann', nperseg=1024)
        amplitude = np.sqrt(power)  # Convert power to amplitude


    #### POWER SPECTRUM
    power = amplitude**2
    amplitude[0] = 0
    power[0] = 0


    # Positive frequencies only
    pos_mask = freqs >= 0
    positive_freqs = freqs[pos_mask]
    positive_amplitude = amplitude[pos_mask]
    positive_power = power[pos_mask]




    # --- FFT: Wall X Velocity with Hanning window (test) ---
    # --- Compute Velocity ---
    velocity = np.diff(wall_x) / dt
    velocity_centered = velocity - np.mean(velocity)
    # --- Apply Hann window to velocity ---
    window = np.hanning(len(velocity_centered))
    velocity_windowed = velocity_centered * window
    # --- FFT: Velocity (raw and windowed) ---
    fft_vel = np.fft.fft(velocity_centered)
    fft_vel_windowed = np.fft.fft(velocity_windowed)
    freqs_vel = np.fft.fftfreq(len(velocity_centered), d=dt)
    amplitude_vel = np.abs(fft_vel) / len(velocity_centered)
    amplitude_vel_windowed = np.abs(fft_vel_windowed) / len(velocity_windowed)
    power_vel = amplitude_vel**2
    power_vel_windowed = amplitude_vel_windowed**2
    amplitude_vel[0] = 0
    amplitude_vel_windowed[0] = 0
    power_vel[0] = 0
    power_vel_windowed[0] = 0
    # Positive Vel frequencies only
    pos_mask_vel = freqs_vel >= 0
    positive_freqs_vel = freqs_vel[pos_mask_vel]
    positive_amplitude_vel = amplitude_vel[pos_mask_vel]
    positive_amplitude_vel_windowed = amplitude_vel_windowed[pos_mask_vel]
    positive_power_vel = power_vel[pos_mask_vel]
    positive_power_vel_windowed = power_vel_windowed[pos_mask_vel]


    # --- savgol smoothing Smooth Amplitude for Peak Finding Only ---
    window_length_savgol = 20
    polyorder_savgol = 10
    smoothed_amplitude = savgol_filter(positive_amplitude, window_length=window_length_savgol, polyorder=polyorder_savgol)
    smootthing_sigma_ampl, dominant_freq = estimate_sigma_from_dominant_freq(wall_x, dt)
    smootthing_sigma_ampl, centroid = estimate_sigma_from_spectral_centroid(wall_x, dt, min_sigma=1.0, max_sigma=10.0, freq_floor=0.01, freq_ceiling=2.0)
    print(f"Estimated smoothing sigma for amplitude: {smootthing_sigma_ampl:.2f}")



    ########### PLOTS #############

    L0, wall_mass_factor, run = extract_params_from_filename(filepath)
    if L0 is None:
        print(f"⚠️ Could not extract parameters from {filepath}")
    # --- Plot: Wall X Velocity ---
    # --- Plot: Original Amplitude Spectrum with Smart Annotated Peaks ---
    fig, ax = plt.subplots(figsize=(12, 5))
    # Plot original amplitude spectrum and Gaussian smoothed curve
    ax.plot(positive_freqs, positive_amplitude, label=f"Original FFT Amplitude with FFT_mode: {FFT_mode} ", color='dodgerblue')
    # Plot smoothed Savgol amplitude for visual reference
    ax.plot(positive_freqs, smoothed_amplitude, label=f"Smoothed FFT SavgolFilter (window={window_length_savgol},polyorder={polyorder_savgol})", color='black', alpha=0.5)
    # Add spectral centroid line
    ax.axvline(centroid, color='orange', linestyle='--', linewidth=2, label=f"Spectral Centroid: {centroid:.2f} Hz")
    #Shaded area around the centroid ± 0.1 Hz
    band_width = 0.1
    ax.axvspan(centroid - band_width, centroid + band_width, color='grey', alpha=0.1, label="Centroid Region ±0.1 Hz")

    annotate_peaks_smart(
        positive_freqs,
        positive_amplitude,
        ax,
        return_mode="refined",  # or "raw" or "smoothed" , "refined is smooth and then insode bandwidth on raw signal "
        label_prefix="Amp Peak",
        yscale='log',
        min_freq=0.001,         # Ignores below 0.5 Hz
        smooth_sigma=smootthing_sigma_ampl,       # Smoothing strength
        prominence = None,  # <<< manually override,       # Tweak depending on your signal
        color='green',
        show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
        smoothed_color='red',         # ← customize color if you want
    )
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude")
    ax.set_title(f"Wall X Amplitude Spectrum with Smart Peak Annotation with Boxlength/2:{L0} and Wall mass factor to particle mass:{wall_mass_factor} (run {run})")
    ax.set_xlim(0, 10)
    ax.set_ylim(auto=True)  # Let matplotlib handle it
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()

    if enable_batch_processing and not show_plots:
        plt.close(fig)  # Don’t show intermediate plots
    else:
        plt.show()
       


    ## PEAK FIND AMPLITUDE
    highest_peak_position = annotate_peaks_smart(
        positive_freqs,
        positive_amplitude,
        ax,
        return_mode="raw",  # or "raw" or "smoothed"
        label_prefix="Amp Peak",
        yscale='log',
        min_freq=0.0,         # Ignores below 0.5 Hz
        smooth_sigma=smootthing_sigma_ampl,       # Smoothing strength
        prominence = None,  # <<< manually override,       # Tweak depending on your signal
        color='green',
        show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
        smoothed_color='red',         # ← customize color if you want)
    )
    if highest_peak_position:
        print(f"Highest Peak Frequency: {highest_peak_position[0]} Hz, Amplitude: {highest_peak_position[1]:.2e}")

    

    #################### --- Plot: Power Spectrum with Peak Detection ---
    smoothed_power = savgol_filter(positive_power, window_length=window_length_savgol, polyorder=polyorder_savgol)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(positive_freqs, positive_power, label="Power Spectrum", color='purple')
    ax.plot(positive_freqs, smoothed_power, label="Smoothed Power Spectrum", color='black', alpha=0.5)
    ax.axvline(centroid, color='orange', linestyle='--', linewidth=2, label=f"Spectral Centroid: {centroid:.2f} Hz")
    ax.axvspan(centroid - band_width, centroid + band_width, color='grey', alpha=0.1, label="Centroid Region ±0.1 Hz")
    

    ## PEAK FIND POWER

    #  Focus only on 2–3 Hz region
    band_min, band_max = 2.0, 3.0
    band_mask = (positive_freqs >= band_min) & (positive_freqs <= band_max)
    diagnostic_freqs = positive_freqs[band_mask]
    diagnostic_power = positive_power[band_mask]

    # Optional: log-transform for better contrast
    # diagnostic_power = np.log10(diagnostic_power + 1e-20)

    highest_peak_power = annotate_peaks_smart(
        positive_freqs,
        positive_power,
        ax,
        return_mode="smoothed",  # or "raw" or "smoothed"
        yscale='linear',  # or 'log' if your data is too compressed
        min_freq=0.0,
        smooth_sigma=1,
        prominence=None,
        color='red',
        show_smoothed_curve=True,
        smoothed_color='darkred'
    )

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Power")
    ax.set_title("Wall X Power Spectrum with Smart Peak Detection")
    ax.set_xlim(0, 10)
    ax.set_ylim(auto=True)
    ax.set_xscale('linear')
    ax.set_yscale('linear')  # ← let matplotlib do the log transform
    ax.grid(True)
    ax.legend()
    plt.tight_layout()

    if enable_batch_processing and not show_plots:
        plt.close(fig)
    else:
        plt.show()

    if highest_peak_power:
        print(f"Highest Power Spectrum Peak: {highest_peak_power[0]:.4f} Hz, Power: {highest_peak_power[1]:.2e}")

    # --- Plot: Velocity FFT ---
    # --- savgol smoothing Smooth Amplitude for Peak Finding Only ---
    window_length_savgol = 20
    polyorder_savgol = 10
    smoothed_velocity = savgol_filter(positive_amplitude_vel, window_length=window_length_savgol, polyorder=polyorder_savgol)

    # --- Plot: Velocity FFT with Smart Peaks ---
    fig, ax = plt.subplots(figsize=(12, 5))
    # Plot original amplitude spectrum and Gaussian smoothed curve
    ax.plot(positive_freqs_vel, positive_amplitude_vel, label="Raw Velocity FFT", color='dodgerblue')
    # Plot smoothed Savgol amplitude for visual reference
    ax.plot(positive_freqs_vel, smoothed_velocity, label=f"Smoothed FFT SavgolFilter (window={window_length_savgol},polyorder={polyorder_savgol})", color='black', alpha=0.5)
    # Add spectral centroid line
    ax.axvline(centroid, color='orange', linestyle='--', linewidth=2, label=f"Spectral Centroid: {centroid:.2f} Hz")
    #Shaded area around the centroid ± 0.1 Hz
    band_width = 0.1
    ax.axvspan(centroid - band_width, centroid + band_width, color='grey', alpha=0.1, label="Centroid Region ±0.1 Hz")

    annotate_peaks_smart(
        positive_freqs_vel,
        positive_amplitude_vel,
        ax,
        label_prefix="Velocity Peak",
        yscale='log',
        min_freq=0.001,
        smooth_sigma=smootthing_sigma_ampl,
        prominence=None,
        color='green',
        show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
        smoothed_color='red'          # ← customize color if you want
    )
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude")
    ax.set_title("Velocity FFT with Smart Peak Detection")
    ax.set_xlim(0, 5)
    ax.set_ylim(auto=True)  # Let matplotlib handle it
    ax.set_xscale('linear')
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()

    if enable_batch_processing and not show_plots:
        plt.close(fig)  # Don’t show intermediate plots
    else:
        plt.show()

    highest_peak_vel = annotate_peaks_smart(
        positive_freqs_vel,
        positive_amplitude_vel,
        ax,
        label_prefix="Velocity Peak",
        yscale='log',
        min_freq=0.001,
        smooth_sigma=smootthing_sigma_ampl,
        prominence=None,
        color='green',
        show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
        smoothed_color='red'          # ← customize color if you want
    )
    if highest_peak_vel:
        print(f"Highest Peak Vel Frequency: {highest_peak_vel[0]} Hz, Amplitude: {highest_peak_vel[1]:.2e}")



    # --- Return the highest peak frequency and amplitude ---
    
    # 🛡️ Safe check
    if highest_peak_position is not None:
        return highest_peak_position[0], highest_peak_position[1]  # freq, amp
    else:
        return None








#### MAIN 
def main():
    if enable_batch_processing:
        filepaths = sorted(glob.glob(os.path.join(folder_path, "*.csv")))
        print("📂 Files found:", filepaths)


        for path in filepaths:
            L0, wall_mass_factor, run = extract_params_from_filename(path)
            if L0 is None:
                print(f"⚠️ Could not extract parameters from {path}")
                continue
            peak = process_file(path)
            print(f"PEAK: {peak},  Processing {path} with L0={L0}, M={wall_mass_factor}, run={run}")
            if peak:
                freq, amp = peak
                results.append((L0, wall_mass_factor, freq, amp))
                print(f"✅ {path} → L0={L0}, M={wall_mass_factor}, freq={freq:.4f} Hz, amp={amp:.2e}")
            else:
                print(f"❌ No peak found for {path}")

        # --- Plot all peaks together ---
        if results:
            results.sort()
            print(results)

            L0_vals, M_vals, freqs, amps = zip(*results)
            # Compute eta
            eta_vals = [np.pi * radius**2 * N / (4 * A * L0) for L0 in L0_vals]
            



            # #############=== Speed of Sound Table ===
            print("\n📊 Calculating Speed of Sound Comparison Table:")
            results_speed = []

            for i in range(len(freqs)):
                L0 = L0_vals[i]
                eta = eta_vals[i]
                freq = freqs[i]
                paper_example_freq = 0.06 # Hz from paper example power peak
                n_value_eq_20 = 2 

                # Equation 20
                c_s_eq_20 = 2 * L0 * freq / n_value_eq_20

                # Equation 21
                c_s_eq21 = 4 * L0 * freq / n_value

                c_21_paper_example = paper_example_freq * 4 * L0 / n_value # Hz from paper example power peak

                # Equation 19
                sqrt_term = np.sqrt((2 * N * mass_particle) / M_vals[i] * mass_particle)
                c_eq19 = freq * (2 * np.pi * L0) / sqrt_term
                
                sqrt_term_paper_example_freq = np.sqrt((2 * N * 1) / M_vals[i] * 1)
                c_eq19_paper_ex_freq =  paper_example_freq * (2 * np.pi * L0) / sqrt_term



                # Theoretical from Equation 26
                def theoretical_cs(a_val):
                    numerator = (1 + eta + 3 * a_val * eta**2 - a_val * eta**3)
                    denominator = (1 - eta)**3
                    return np.sqrt((2 * k_B * T / mass_particle) * (numerator / denominator))

                cs_SPT = theoretical_cs(0)           # a = 0
                cs_Henderson = theoretical_cs(0.125) # a = 0.125

                results_speed.append((L0, eta, c_eq19, c_s_eq_20, c_s_eq21, cs_SPT, cs_Henderson, c_eq19_paper_ex_freq, c_21_paper_example))

            # Output DataFrame
            speed_df = pd.DataFrame(
                results_speed,
                columns=["L0", "η", "c_s_19",'c_s_eq_20', "c_s_21", "c_s_SPT", "c_s_Henderson", "c_s_19_paper_example_freq",'c_s_21_paper_example_freq']
            )
            print(speed_df.to_string(index=False, float_format="{:.3f}".format))





            ######## PLOTS   ###########################             
            # # 3D Plot: Peak Frequency vs L0 and M

            L0_vals_arr = np.array(L0_vals)
            M_vals_arr = np.array(M_vals)
            freqs_arr = np.array(freqs)

            # Create grid
            L0_unique = np.unique(L0_vals_arr)
            M_unique = np.unique(M_vals_arr)
            L0_grid, M_grid = np.meshgrid(L0_unique, M_unique)

            # Fill Z grid with matching frequencies
            Z = np.empty_like(L0_grid, dtype=float)
            Z[:] = np.nan

            for i, L0 in enumerate(L0_unique):
                for j, M in enumerate(M_unique):
                    # Find index matching this (L0, M)
                    matches = [(k, f) for k, (l, m, f) in enumerate(zip(L0_vals, M_vals, freqs)) if l == L0 and m == M]
                    if matches:
                        Z[j, i] = matches[0][1]  # take first match

            # Plot
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(L0_grid, M_grid, Z, cmap=cm.viridis, edgecolor='k', linewidth=0.5, alpha=0.9)

            ax.set_xlabel("Box Length L0")
            ax.set_ylabel("Wall Mass Factor M")
            ax.set_zlabel("Peak Frequency [Hz]")
            ax.set_title("3D Surface of Peak Frequencies")
            fig.colorbar(surf, shrink=0.5, aspect=10, label="Frequency")

            plt.tight_layout()
            plt.show()



            # Marker setup
            unique_L0s = sorted(set(L0_vals))
            marker_styles = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'H', '8']
            L0_to_marker = {L0: marker_styles[i % len(marker_styles)] for i, L0 in enumerate(unique_L0s)}

            vmin = min(M_vals)
            vmax = max(M_vals)

            plt.figure(figsize=(12, 6))
            for i in range(len(freqs)):
                marker = L0_to_marker[L0_vals[i]]
                plt.scatter(freqs[i], amps[i], marker=marker, c=[M_vals[i]],
                            cmap='viridis', s=70, edgecolors='k', alpha=0.9,
                            vmin=vmin, vmax=vmax)
                plt.annotate(f"L0={L0_vals[i]}, M={M_vals[i]}\nη={eta_vals[i]:.3f}",
                            (freqs[i], amps[i]), fontsize=8, alpha=0.7)

            cbar = plt.colorbar()
            cbar.set_label("Wall Mass Factor")

            # Legend for L0 values
            legend_elements = [plt.Line2D([0], [0], marker=L0_to_marker[L0], color='w',
                                        markerfacecolor='gray', markeredgecolor='k',
                                        markersize=8, label=f"L0={L0}")
                            for L0 in unique_L0s]
            plt.legend(handles=legend_elements, title="Box Length (L0)", loc='upper right', fontsize=8)

            plt.title(plot_title)
            plt.xlabel("Peak Frequency [Hz]")
            plt.ylabel("Amplitude")
            plt.yscale("log")
            plt.xlim(0, max(freqs) * 1.05)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()





            #  === GET Speed of sound by SLOPE!!! of PLOT and interpolate Peak Frequency(y-axis) vs K/2π(L0 - 1) (x_axis) ===:

            vmin = min(M_vals)
            vmax = max(M_vals) 

            #### PLOT and interpolate for slPeak Frequency vs K/2π(L0 - 1) ===
            plt.figure(figsize=(12, 6))
            for i in range(len(freqs)):

                marker = L0_to_marker[L0_vals[i]]
                wall_mass_factor = M_vals[i]
                K_root = find_transcendental_root(wall_mass_factor=wall_mass_factor, N=50)
                K = K_root

                #plot_transcendental(wall_mass_factor, N)
                print(f"K root for wall_mass_factor={wall_mass_factor}: {K_root:.6f}")

                x_values = [K / (2 * np.pi * (L0 - 1)) for L0 in L0_vals]
                plt.scatter(x_values[i], freqs[i], label=None, marker=marker, c=[M_vals[i]],
                            cmap='viridis', s=70, edgecolors='k', alpha=0.9,
                            vmin=vmin, vmax=vmax)
                plt.annotate(f"L0={L0_vals[i]}, M={M_vals[i]}\nη={eta_vals[i]:.3f}",
                            (x_values[i], freqs[i]), fontsize=8, alpha=0.7)

            cbar = plt.colorbar()
            cbar.set_label("Wall Mass Factor")
            # Legend for L0 marker styles
            legend_elements = [plt.Line2D([0], [0], marker=L0_to_marker[L0], color='w',
                                        markerfacecolor='gray', markeredgecolor='k',
                                        markersize=8, label=f"L0={L0}")
                            for L0 in unique_L0s]
            plt.legend(handles=legend_elements, title="Box Length (L0)", loc='upper right', fontsize=8)

            plt.xlabel(r"$\frac{K}{2\pi(L_0 - 1)}$")
            plt.ylabel("Peak Frequency [Hz]")
            plt.title("Fundamental Frequency vs Theoretical K-term")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()





            #### 3rd plot: Peak Frequency vs Wall Mass Factor
            unique_L0s = sorted(set(L0_vals))
            unique_Ms = sorted(set(M_vals))

            colors = viridis(np.linspace(0, 1, len(unique_L0s)))
            L0_to_color = {L0: colors[i] for i, L0 in enumerate(unique_L0s)}

            # Helper function for x-jitter based on theoretical K-term
            def jitter_x(L0):
                return K / (2 * pi * (L0 - 1))


            # Plot41: Frequency vs Wall Mass
            plt.figure(figsize=(10, 6))
            for L0 in unique_L0s:
                xs = [M_vals[i] for i in range(len(M_vals)) if L0_vals[i] == L0]
                ys = [freqs[i] for i in range(len(freqs)) if L0_vals[i] == L0]
                labels = [f"K/(2π(L0-1))={jitter_x(L0):.3f}" for i in range(len(freqs)) if L0_vals[i] == L0]
                plt.plot(xs, ys, marker='o', label=f'L0={L0}', color=L0_to_color[L0])
                for x, y, label in zip(xs, ys, labels):
                    plt.annotate(label, (x, y), fontsize=8, alpha=0.7)

            plt.xlabel("Wall Mass Factor")
            plt.ylabel("Peak Frequency [Hz]")
            plt.title("Peak Frequency vs Wall Mass Factor (grouped by Box Length)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()



            # === 4th plot: Peak Amplitude vs Wall Mass Factor (grouped by L0) ===
            plt.figure(figsize=(10, 6))
            for L0 in unique_L0s:
                xs = [M_vals[i] for i in range(len(M_vals)) if L0_vals[i] == L0]
                ys = [amps[i] for i in range(len(amps)) if L0_vals[i] == L0]
                labels = [f"K/(2π(L0-1))={jitter_x(L0):.3f}" for i in range(len(amps)) if L0_vals[i] == L0]
                
                plt.plot(xs, ys, marker='o', linestyle='-', label=f"L0={L0}", color=L0_to_color[L0])
                for x, y, label in zip(xs, ys, labels):
                    plt.annotate(label, (x, y), fontsize=8, textcoords="offset points", xytext=(0, 5), ha='center')

            plt.xlabel("Wall Mass Factor")
            plt.ylabel("Peak Amplitude")
            plt.title("Peak Amplitude vs Wall Mass Factor (grouped by Box Length)")
            plt.yscale("log")
            plt.grid(True, alpha=0.3)
            plt.legend(title="L0")
            plt.tight_layout()
            plt.show()


            # 5th and 6th plots: Grouped by L0, x-axis = L0, color = Mass

            unique_Ms = sorted(set(M_vals))
            unique_L0s = sorted(set(L0_vals))
            colors_mass = plt.cm.viridis(np.linspace(0, 1, len(unique_Ms)))
            M_to_color = {M: colors_mass[i] for i, M in enumerate(unique_Ms)}

            # Plot 5: Peak Frequency vs Box Length (colored by Wall Mass)
            plt.figure(figsize=(10, 6))
            for M in unique_Ms:
                xs = [L0_vals[i] for i in range(len(freqs)) if M_vals[i] == M]
                ys = [freqs[i] for i in range(len(freqs)) if M_vals[i] == M]
                plt.plot(xs, ys, marker='o', linestyle='-', label=f"M={M}", color=M_to_color[M])
            plt.xlabel("Box Length L0")
            plt.ylabel("Peak Frequency [Hz]")
            plt.title("Peak Frequency vs Box Length (colored by Mass)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()


            # Plot 6: Peak Amplitude vs Box Length (colored by Wall Mass)
            plt.figure(figsize=(10, 6))
            for M in unique_Ms:
                xs = [L0_vals[i] for i in range(len(amps)) if M_vals[i] == M]
                ys = [amps[i] for i in range(len(amps)) if M_vals[i] == M]
                plt.plot(xs, ys, marker='o', linestyle='-', label=f"M={M}", color=M_to_color[M])
            plt.xlabel("Box Length L0")
            plt.ylabel("Peak Amplitude")
            plt.title("Peak Amplitude vs Box Length (colored by Mass)")
            plt.yscale("log")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()


            # 7th and 8th plots: Grouped by M, x-axis = Mass, color = L0
            colors_L0 = plt.cm.plasma(np.linspace(0, 1, len(unique_L0s)))
            L0_to_color = {L0: colors_L0[i] for i, L0 in enumerate(unique_L0s)}

            # Plot 7: Peak Frequency vs Mass (colored by L0)
            plt.figure(figsize=(10, 6))
            for L0 in unique_L0s:
                xs = [M_vals[i] for i in range(len(freqs)) if L0_vals[i] == L0]
                ys = [freqs[i] for i in range(len(freqs)) if L0_vals[i] == L0]
                plt.plot(xs, ys, marker='o', linestyle='-', label=f"L0={L0}", color=L0_to_color[L0])
            plt.xlabel("Wall Mass Factor")
            plt.ylabel("Peak Frequency [Hz]")
            plt.title("Peak Frequency vs Wall Mass Factor (colored by Box Length)")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()

            # Plot 8: Peak Amplitude vs Mass (colored by L0)
            plt.figure(figsize=(10, 6))
            for L0 in unique_L0s:
                xs = [M_vals[i] for i in range(len(amps)) if L0_vals[i] == L0]
                ys = [amps[i] for i in range(len(amps)) if L0_vals[i] == L0]
                plt.plot(xs, ys, marker='o', linestyle='-', label=f"L0={L0}", color=L0_to_color[L0])
            plt.xlabel("Wall Mass Factor")
            plt.ylabel("Peak Amplitude")
            plt.title("Peak Amplitude vs Wall Mass Factor (colored by Box Length)")
            plt.yscale("log")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()


        else:
            print("⚠️ No results collected.")

    ############# END BATCH PROCESSING #############
    else:
        print("🔬 Batch mode disabled — analyzing single file:")
        # --- Process a single file ---
        process_file("wall_position.csv")  # You can modify this to any filename
        #process_file(os.path.join(folder_path, single_plot_filename))  # You can modify this to any filename
        



if __name__ == "__main__":
    main()