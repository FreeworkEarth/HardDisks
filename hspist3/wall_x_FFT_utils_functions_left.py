
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
from scipy.signal import butter, filtfilt
from scipy.stats import linregress


from scipy.optimize import root_scalar
from matplotlib.cm import viridis
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#### SINGLE FILE FFT ANALYSIS ####
# --- Load the data ---
df = pd.read_csv("wall_position.csv")
print("Columns:", df.columns)

# Clean column names just in case
df.columns = df.columns.str.strip()

# Extract time and wall position
time = df['Time'].values
wall_x = df['Wall_X'].values

# --- Compress and print wall positions for inspection ---
compression_factor = 500
sampled_indices = np.arange(0, len(wall_x), compression_factor)
compressed_df = pd.DataFrame({
    "Time": time[sampled_indices],
    "Wall_X": wall_x[sampled_indices]
})
print("Compressed Wall Position Data:")
print(compressed_df.to_string(index=False))

# Calculate time step
dt = np.mean(np.diff(time))
print(f"Time step (dt): {dt:.6f} seconds")
sampling_rate = 1 / dt

# --- Center wall position ---
wall_x_centered = wall_x - np.mean(wall_x)





q

##################################################

##### PLOTS
# --- Plot: Amplitude Spectrum ---

# --- savgol smoothing Smooth Amplitude for Peak Finding Only ---
window_length_savgol = 20
polyorder_savgol = 10
smoothed_amplitude = savgol_filter(positive_amplitude, window_length=window_length_savgol, polyorder=polyorder_savgol)
smootthing_sigma_ampl, dominant_freq = estimate_sigma_from_dominant_freq(wall_x, dt)
smootthing_sigma_ampl, centroid = estimate_sigma_from_spectral_centroid(wall_x, dt, min_sigma=1.0, max_sigma=10.0, freq_floor=0.01, freq_ceiling=2.0)
print(f"Estimated smoothing sigma for amplitude: {smootthing_sigma_ampl:.2f}")


# --- Plot: Original Amplitude Spectrum with Smart Annotated Peaks ---
fig, ax = plt.subplots(figsize=(12, 5))
# Plot original amplitude spectrum and Gaussian smoothed curve
ax.plot(positive_freqs, positive_amplitude, label="Original Amplitude", color='dodgerblue')
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
ax.set_title("Wall X Amplitude Spectrum with Smart Peak Annotation")
ax.set_xlim(0, 2.5)
ax.set_ylim(auto=True)  # Let matplotlib handle it
ax.set_xscale('linear')
ax.set_yscale('log')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()


highest_peak = annotate_peaks_smart(
    positive_freqs,
    positive_amplitude,
    ax,
    label_prefix="Amp Peak",
    yscale='log',
    min_freq=0.01,         # Ignores below 0.5 Hz
    smooth_sigma=2,       # Smoothing strength
    prominence = None,  # <<< manually override,       # Tweak depending on your signal
    color='green',
    show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
    smoothed_color='red',         # ← customize color if you want)
)
if highest_peak:
    print(f"Highest Peak Frequency: {highest_peak[0]} Hz, Amplitude: {highest_peak[1]:.2e}")





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
ax.set_xlim(0, 2.5)
ax.set_ylim(auto=True)  # Let matplotlib handle it
ax.set_xscale('linear')
ax.set_yscale('log')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()


highest_peak_vel = annotate_peaks_smart(
    positive_freqs_vel,
    positive_amplitude_vel,
    ax,
    label_prefix="Velocity Peak",
    yscale='log',
    min_freq=0.001,
    smooth_sigma=2,
    prominence=None,
    color='green',
    show_smoothed_curve=True,      # ← enable internal smoothed FFT curve
    smoothed_color='red'          # ← customize color if you want
)
if highest_peak_vel:
    print(f"Highest Peak Vel Frequency: {highest_peak_vel[0]} Hz, Amplitude: {highest_peak_vel[1]:.2e}")





def annotate_peaks(x, y, ax, label_prefix="Peak", yscale='log', n_peaks=3, color='black'):
    peaks, _ = find_peaks(y, height=np.max(y) * 0.1)
    sorted_peaks = sorted(peaks, key=lambda i: y[i], reverse=True)[:n_peaks]

    legend_handles = []
    legend_labels = []

    for i, peak in enumerate(sorted_peaks):
        f = x[peak]
        a = y[peak]

        ax.axvline(f, color='red', linestyle='--', alpha=0.6)
        ax.axhline(a, color='blue', linestyle='--', alpha=0.6)
        point, = ax.plot(f, a, 'o', color=color)

        label = f"{label_prefix} {i+1}: {f:.2f} Hz, {a:.2e}" if yscale == 'log' else f"{label_prefix} {i+1}: {f:.2f} Hz, {a:.2f}"
        legend_handles.append(point)
        legend_labels.append(label)

    ax.legend(legend_handles, legend_labels, loc='upper right')



def annotate_peaks_smart1(freqs, amplitudes, ax, label_prefix="Peak",
                         yscale='linear', min_freq=0.0, smooth_sigma=2,
                         prominence=None, color='red', max_peaks=3,
                         show_smoothed_curve=True, smoothed_color='red'):
    """
    Smart peak detector and annotator with optional smoothed curve plot.
    - Detects peaks on smoothed data.
    - Shows peak markers on smoothed curve.
    - Can also optionally overlay smoothed FFT for visual clarity.
    """

    # Save full original arrays for plotting
    freqs_full = freqs
    amplitudes_full = amplitudes

    # Smooth entire spectrum first
    smoothed_full = gaussian_filter1d(amplitudes_full, sigma=smooth_sigma)

    # Apply min_freq mask for peak detection only
    valid_idx = (freqs_full > min_freq)
    freqs = freqs_full[valid_idx]
    smoothed_amp = smoothed_full[valid_idx]
    amplitudes = amplitudes_full[valid_idx]

    # Smooth FFT for stable peak detection
    smoothed_amp = gaussian_filter1d(amplitudes, sigma=smooth_sigma)

    # Optional: plot smoothed curve for visual reference
    if show_smoothed_curve:
        #ax.plot(freqs, smoothed_amp, label="Smoothed FFT Gaussian1D", color=smoothed_color, linestyle='--', alpha=0.6)
        ax.plot(freqs, smoothed_amp, label=f"Smoothed FFT Gaussian1D(σ={smooth_sigma})", color=smoothed_color, linestyle='--', alpha=0.6)    
    # Adaptive prominence
    if prominence is None:
        sharpness = np.max(np.abs(np.gradient(smoothed_amp))) / (np.mean(smoothed_amp) + 1e-12)
        base_prominence = 0.02 * sharpness
    else:
        base_prominence = prominence

    # Avoid low-garbage peaks
    floor = np.max(smoothed_amp) * 0.005 # 0.5% instead of 5%
    floor = np.median(smoothed_amp) + np.std(smoothed_amp) * 0.5 # more dynamically based on curve
    base_prominence = max(base_prominence, floor)

    # Search for a clean number of peaks
    for scale in [1.0, 1.5, 2.0, 3.0]:
        test_prom = base_prominence * scale
        peaks, properties = find_peaks(smoothed_amp, prominence=test_prom)
        if len(peaks) <= max_peaks and len(peaks) > 0:
            break
    else:
        print("⚠️ No clear peaks found after prominence tuning.")
        return

    print(f"✅ {len(peaks)} peak(s) detected with final prominence > {test_prom:.3f}")

    # Sort and process most prominent
    prominences = properties["prominences"]
    sorted_indices = np.argsort(prominences)[::-1][:max_peaks]
    selected_peaks = peaks[sorted_indices]

    print(f"Max smoothed amplitude: {np.max(smoothed_amp):.2e}")
    print(f"Base prominence estimate: {base_prominence:.2e}")


    for i, smoothed_idx in enumerate(selected_peaks):
        freq = freqs[smoothed_idx]
        amp_smoothed = smoothed_amp[smoothed_idx]

        # Use raw amplitude for annotation
        nearest_raw_idx = np.argmin(np.abs(freqs - freq))
        amp_raw = amplitudes[nearest_raw_idx]

        label = f"{label_prefix} {i+1}" if len(selected_peaks) > 1 else label_prefix
        label_text = f"{label} ({freq:.2f} Hz, {amp_raw:.2e})"

        # Vertical line and peak marker
        ax.axvline(freq, linestyle='--', color=color, alpha=0.6)
        ax.plot(freq, amp_smoothed, 'o', color=color, label=label_text)  # Use smoothed amp for marker

        ax.annotate(f"{freq:.2f} Hz\n{amp_raw:.2e}",
                    xy=(freq, amp_smoothed),
                    xytext=(5, 10),
                    textcoords='offset points',
                    fontsize=9, color=color)

        print(f"📍 Annotated {label} at {freq:.3f} Hz with amp {amp_raw:.2e}")

    if yscale == 'log':
        ax.set_yscale('log')

    # Clean up legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(), dict(zip(labels, handles)).keys(), fontsize=9)



def annotate_peaks_smart2(freqs, amplitudes, ax, label_prefix="Peak",
                         yscale='linear', min_freq=0.0, smooth_sigma=2,
                         prominence=None, color='red', max_peaks=3,
                         show_smoothed_curve=True, smoothed_color='red'):
    """
    Smart peak detector and annotator with optional smoothed curve plot.
    - Detects peaks on smoothed data.
    - Shows peak markers on smoothed curve.
    - Can also optionally overlay smoothed FFT for visual clarity.
    """

    # Smooth the full FFT for consistent visual overlay
    smoothed_full = gaussian_filter1d(amplitudes, sigma=smooth_sigma)

    # Plot full smoothed curve (starts at 0 Hz)
    if show_smoothed_curve:
        ax.plot(freqs, smoothed_full, label=f"Smoothed FFT Gaussian1D (σ={smooth_sigma})", color=smoothed_color, linestyle='--', alpha=0.6)

    # Apply min_freq mask for peak detection only
    valid_idx = (freqs > min_freq)
    freqs_valid = freqs[valid_idx]
    amps_valid = amplitudes[valid_idx]
    smoothed_valid = smoothed_full[valid_idx]

    # Adaptive prominence
    if prominence is None:
        sharpness = np.max(np.abs(np.gradient(smoothed_valid))) / (np.mean(smoothed_valid) + 1e-12)
        base_prominence = 0.02 * sharpness
    else:
        base_prominence = prominence

    # Avoid low-garbage peaks
    floor = np.max(smoothed_valid) * 0.05
    base_prominence = max(base_prominence, floor)


    # Search for a clean number of peaks
    for scale in [1.0, 1.5, 2.0, 3.0]:
        test_prom = base_prominence * scale
        peaks, properties = find_peaks(smoothed_valid, prominence=test_prom)
        if 0 < len(peaks) <= max_peaks:
            break
    else:
        print("⚠️ No clear peaks found after prominence tuning.")
        return

    print(f"✅ {len(peaks)} peak(s) detected with final prominence > {test_prom:.3f}")

    # Sort and process most prominent
    prominences = properties["prominences"]
    sorted_indices = np.argsort(prominences)[::-1][:max_peaks]
    selected_peaks = peaks[sorted_indices]

    for i, idx in enumerate(selected_peaks):
        freq = freqs_valid[idx]
        amp_smoothed = smoothed_valid[idx]

        # Use raw amplitude for annotation
        nearest_raw_idx = np.argmin(np.abs(freqs - freq))
        amp_raw = amplitudes[nearest_raw_idx]

        label = f"{label_prefix} {i+1}" if len(selected_peaks) > 1 else label_prefix
        label_text = f"{label} ({freq:.2f} Hz, {amp_raw:.2e})"

        # Vertical line and peak marker
        ax.axvline(freq, linestyle='--', color=color, alpha=0.6)
        ax.plot(freq, amp_smoothed, 'o', color=color, label=label_text)

        ax.annotate(f"{freq:.2f} Hz\n{amp_raw:.2e}",
                    xy=(freq, amp_smoothed),
                    xytext=(5, 10),
                    textcoords='offset points',
                    fontsize=9, color=color)

        print(f"📍 Annotated {label} at {freq:.3f} Hz with amp {amp_raw:.2e}")

    if yscale == 'log':
        ax.set_yscale('log')

    # Clean up legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(), dict(zip(labels, handles)).keys(), fontsize=9)









# Plot to visually compare
plt.figure(figsize=(12, 4))
plt.plot(positive_freqs, smoothed_amplitude, color='darkorange', label="Smoothed")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.title("Smoothed vs Original FFT Amplitude")
plt.legend()
plt.grid(True)
plt.xlim(0, 10)
plt.yscale('log')
plt.tight_layout()
plt.show()


# --- Plot: Amplitude Spectrum with Smoothed Peaks ---
# Call with smoothed_amplitude instead of original
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(positive_freqs, positive_amplitude, label="Amplitude Spectrum")
# Use smoothed amplitude for finding peaks but original for display
annotate_peaks_smart(
    positive_freqs,
    positive_amplitude,
    ax,
    label_prefix="Amp Peak",
    yscale='log',
    min_freq=0.5,         # Ignores below 0.5 Hz
    smooth_sigma=2,       # Smoothing strength
    prominence=0.03       # Tweak depending on your signal
)
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Amplitude")
ax.set_title("Amplitude Spectrum with Smoothed Peak Detection")
ax.set_xscale('linear')
ax.set_yscale('log')
ax.set_xlim(0, 10)
ax.set_ylim(1e-5, 3e1)
ax.grid(True)
plt.tight_layout()
plt.show()






# --- Plot: Windowed Velocity FFT ---
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(positive_freqs_vel, positive_amplitude_vel_windowed, label="Windowed FFT")
annotate_peaks(positive_freqs_vel, positive_amplitude_vel_windowed, ax, label_prefix="Win Vel Amp", yscale='log')
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Amplitude")
ax.set_title("FFT of Windowed Wall X Velocity")
ax.set_xscale('linear')
ax.set_yscale('log')
ax.set_xlim(0, 10)
ax.set_ylim(1e-5, 3*1e1)
ax.grid(True)
plt.tight_layout()
plt.show()





def annotate_spectrogram_peaks(freqs, spectrum, ax, n_peaks=3, color='white'):
    peaks, _ = find_peaks(spectrum, height=np.max(spectrum)*0.1)
    sorted_peaks = sorted(peaks, key=lambda i: spectrum[i], reverse=True)[:n_peaks]

    for i, peak in enumerate(sorted_peaks):
        f = freqs[peak]
        a = spectrum[peak]
        ax.axhline(f, color=color, linestyle='--', alpha=0.8, linewidth=1)
        ax.text(
            0.95, f, f"{f:.2f} Hz\n{a:.2e}", 
            color=color, ha='right', va='bottom', fontsize=9,
            transform=ax.get_yaxis_transform()
        )

# --- Spectrogram calculation for centered velocity ---
f_spec, t_spec, Sxx = spectrogram(velocity_centered, fs=sampling_rate, window='hann', nperseg=1024, noverlap=512)

# --- Spectrogram calculation for windowed velocity ---
f_spec_win, t_spec_win, Sxx_windowed = spectrogram(velocity_windowed, fs=sampling_rate, window='hann', nperseg=1024, noverlap=512)

# --- Spectrogram (Wall X Velocity, Centered) ---
fig, ax = plt.subplots(figsize=(12, 5))
img = ax.pcolormesh(t_spec, f_spec, 10*np.log10(Sxx), shading='gouraud', cmap='viridis')
ax.set_title("Spectrogram of Wall X Velocity (Centered)")
ax.set_ylabel("Frequency [Hz]")
ax.set_xlabel("Time [s]")
fig.colorbar(img, ax=ax, label='Power [dB]')

# Average across time to detect peaks
power_avg = np.mean(Sxx, axis=1)
annotate_spectrogram_peaks(f_spec, power_avg, ax, n_peaks=3, color='white')

plt.tight_layout()
plt.show()

# --- Spectrogram (Windowed Velocity) ---
fig, ax = plt.subplots(figsize=(12, 5))
img = ax.pcolormesh(t_spec, f_spec, 10*np.log10(Sxx_windowed), shading='gouraud', cmap='plasma')
ax.set_title("Spectrogram of Windowed Wall X Velocity")
ax.set_ylabel("Frequency [Hz]")
ax.set_xlabel("Time [s]")
fig.colorbar(img, ax=ax, label='Power [dB]')

power_avg_windowed = np.mean(Sxx_windowed, axis=1)
annotate_spectrogram_peaks(f_spec, power_avg_windowed, ax, n_peaks=3, color='white')

plt.tight_layout()
plt.show()

