import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from pathlib import Path
from scipy import stats

# Configuration constants
overscan_start = 6152
y_overscan_start = 1536
prescan = 8
running_window = 40
col_threshold = 6
cti_window = 100

def subtract_pedestal(data):
    row_medians = np.median(data[:, overscan_start:], axis=1)
    subtracted = data - row_medians[:, np.newaxis]
    noise_os = np.std(subtracted[:, overscan_start:])
    mean_os = np.mean(subtracted[:, overscan_start:])
    threshold = mean_os + 5 * noise_os
    valid_pixels = subtracted[:, overscan_start:][subtracted[:, overscan_start:] < threshold]
    print("Noise [ADU] : ", round(np.std(valid_pixels), 1))
    return subtracted

def compute_median_images(fz_files):
    all_images = [[] for _ in range(4)]
    all_noises = [[] for _ in range(4)]

    for f in fz_files:
        with fits.open(f) as file:
            for channel in range(1, 5):
                raw = np.array(file[channel].data)
                print(f"Loaded {f.name} channel {channel-1}")
                sub = subtract_pedestal(raw)
                all_images[channel - 1].append(sub)
                noise_os = np.std(sub[:, overscan_start:])
                mean_os = np.mean(sub[:, overscan_start:])
                threshold = mean_os + 5 * noise_os
                valid = sub[:, overscan_start:][sub[:, overscan_start:] < threshold]
                all_noises[channel - 1].append(np.std(valid))

    median_images = [np.median(np.stack(imgs), axis=0) for imgs in all_images]
    median_noises = [round(np.median(noises), 1) for noises in all_noises]
    return median_images, median_noises

def running_median_absolute_deviation(arr, window_size):
    medians, mads = [], []
    arr = np.where(np.ma.getmask(arr), np.nan, arr)
    for i in range(len(arr)):
        start = max(0, i - window_size // 2)
        end = min(len(arr), i + window_size // 2 + 1)
        medians.append(np.nanmedian(arr[start:end]))
        mads.append(stats.median_abs_deviation(arr[start:end], nan_policy='omit'))
    return np.ma.masked_array(medians), np.ma.masked_array(mads)

def find_col_traps(image, threshold=col_threshold):
    masked = np.ma.masked_array(image)
    masked[:, :prescan] = np.ma.masked
    masked[:, overscan_start - prescan:] = np.ma.masked
    col_medians = np.ma.median(masked, axis=0)
    mask = np.ma.getmaskarray(col_medians)
    running_medians, running_mads = running_median_absolute_deviation(col_medians, window_size=running_window)
    high = running_medians + threshold * running_mads
    low = running_medians - threshold * running_mads
    traps = ~mask & ((col_medians > high) | (col_medians < low))
    return high, low, traps, np.where(traps)[0]

def calculate_cti_x(col_medians):
    pixels = col_medians[overscan_start - cti_window:overscan_start]
    overscan_val = col_medians[overscan_start]
    overscan = col_medians[overscan_start + cti_window:]
    px_mean = np.mean(pixels)
    px_std = np.std(pixels)
    os_mean = np.mean(overscan)
    os_std = np.std(overscan)
    val = px_mean - os_mean
    unc_val = np.sqrt(px_std**2 + os_std**2)
    os_step = overscan_val - os_mean
    os_unc = np.sqrt(2) * os_std
    cti = os_step / val
    cti_unc = cti * np.sqrt((unc_val / val)**2 + (os_unc / os_step)**2)
    return cti * 100, cti_unc * 100

def calculate_cti_y(row_medians):
    pixels = row_medians[y_overscan_start - cti_window:y_overscan_start]
    overscan_val = row_medians[y_overscan_start + 1]
    overscan = row_medians[y_overscan_start + 1:]
    px_mean = np.mean(pixels)
    px_std = np.std(pixels)
    os_mean = np.mean(overscan)
    os_std = np.std(overscan)
    val = px_mean - os_mean
    unc_val = np.sqrt(px_std**2 + os_std**2)
    os_step = overscan_val - os_mean
    os_unc = np.sqrt(2) * os_std
    cti = os_step / val
    cti_unc = cti * np.sqrt((unc_val / val)**2 + (os_unc / os_step)**2)
    return cti * 100, cti_unc * 100

def save_image_detailed(images, ccd_image_name, log_file):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    axs = axs.flatten()

    for i, (image, channel, noise) in enumerate(images):
        ax = axs[i]
        im = ax.imshow(image, cmap='cividis', origin='lower', vmin=-1, vmax=10, aspect='auto')
        ax.set_title(f'{ccd_image_name} – Channel {channel - 1}', fontsize=14, fontweight="bold")
        ax.set_xlabel("Column Number", fontsize=12)
        ax.set_ylabel("Row Number", fontsize=12)

        divider = make_axes_locatable(ax)
        axtop = divider.append_axes("bottom", size="30%", pad=0.8, sharex=ax)
        axright = divider.append_axes("right", size="15%", pad=0.8, sharey=ax)
        cax = divider.append_axes("right", size="5%", pad=0.3)

        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label("Signal [ADU]", fontsize=12)

        col_medians = np.ma.median(image, axis=0)
        row_medians = np.ma.median(image, axis=1)
        col_medians[:prescan] = np.nan

        condition_high, condition_low, _, trap_locations = find_col_traps(image)

        axtop.step(np.arange(len(col_medians)), col_medians, where='mid', label='Median')
        axtop.step(np.arange(len(condition_high)), condition_high, color='red', linestyle='dashed', label='Median + 6 MAD')
        axtop.step(np.arange(len(condition_low)), condition_low, color='red', linestyle='dotted', label='Median - 6 MAD')
        axtop.set_ylabel("ADU", fontsize=10)
        axtop.set_xlabel("Column Number", fontsize=10)
        axtop.set_title("Column Median", fontsize=10)

        legend = axtop.legend(loc='lower left', bbox_to_anchor=(1.0, 0.5), borderaxespad=0., fontsize=8, frameon=True)

        axright.step(row_medians, np.arange(image.shape[0]), where='mid')
        axright.set_xlabel("ADU", fontsize=10)
        axright.set_title("Row Median", fontsize=10)

        cti_x, cti_x_err = calculate_cti_x(col_medians)
        cti_y, cti_y_err = calculate_cti_y(row_medians)

        log_file.write("-" * 72 + "\n")
        log_file.write(f"---------------------     Chhannel : {channel - 1}  -------------------------\n")
        log_file.write(f"Median Noise : {noise} ADU\n")
        log_file.write(f"Column Defects Location: {trap_locations.tolist()}\n")
        log_file.write(f"Number of Column Defects: {len(trap_locations)}\n")
        log_file.write(f"CTI X [%]: {cti_x:.2f} ± {cti_x_err:.2f}\n")
        log_file.write(f"CTI Y [%]: {cti_y:.2f} ± {cti_y_err:.2f}\n")

    fig.suptitle(f'Column Defect Maps – {ccd_image_name}', fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{ccd_image_name}.png", dpi=300)
    plt.close()

def main():
    args = sys.argv
    if len(args) < 3:
        print(f"Usage: {args[0]} <fz_file1> [fz_file2 ...] <ccd_image_name>")
        sys.exit(1)

    *fz_files, ccd_image_name = args[1:]
    fz_files = [Path(f) for f in fz_files if Path(f).suffix == '.fz']

    for f in fz_files:
        if not f.exists():
            print(f"[ERROR] File not found: {f}")
            sys.exit(1)

    median_images, median_noises = compute_median_images(fz_files)
    images = [(median_images[i], i + 1, median_noises[i]) for i in range(4)]

    log_file_path = Path(f"{ccd_image_name}.log")
    with log_file_path.open('w') as log_file:
        save_image_detailed(images, ccd_image_name, log_file)

    print(f"[DONE] Saved: {ccd_image_name}.png and {ccd_image_name}.log")
    plt.show()

if __name__ == "__main__":
    main()
