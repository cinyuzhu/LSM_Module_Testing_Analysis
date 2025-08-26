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

                # Noise
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
    print("Column Defects Location:", np.where(traps)[0])
    print("Number of Column Defects:", np.sum(traps))
    return high, low, traps

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
    print(f"CTI X [%]: {cti * 100:.2f} ± {cti_unc * 100:.2f}")
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
    print(f"CTI Y [%]: {cti * 100:.2f} ± {cti_unc * 100:.2f}")
    return cti * 100, cti_unc * 100

def show_image_detailed(image, channel, noise):
    print("########################################################################")
    print(f"#####################     channel : {channel - 1}     #########################")
    print(f"Median Noise : {noise} ADU")

    fig, ax = plt.subplots()
    im = ax.imshow(image, cmap='viridis', origin='lower', vmin=-1, vmax=80)
    ax.set_title(f'Combined Image [channel {channel - 1}]', fontweight="bold")
    ax.set_xlabel("Column Number")
    ax.set_ylabel("Row Number")
    divider = make_axes_locatable(ax)
    axtop = divider.append_axes("bottom", size="20%", pad=0.7, sharex=ax)
    axright = divider.append_axes("right", size="10%", pad=0.5, sharey=ax)

    high, low, _ = find_col_traps(image)
    col_medians = np.median(image, axis=0)
    row_medians = np.median(image, axis=1)
    col_medians[:prescan] = np.nan

    axtop.step(np.arange(len(col_medians)), col_medians, where='mid')
    axtop.step(np.arange(len(high)), high, color='red', where='mid')
    axtop.step(np.arange(len(low)), low, color='red', where='mid')

    axright.step(row_medians, np.arange(len(row_medians)), where='mid')
    axtop.set_ylabel("[ADU]")
    axright.set_xlabel("[ADU]")
    axtop.set_title("Column-wise median", fontweight="bold")
    axright.set_title("Row-wise median", fontweight="bold")

    calculate_cti_x(col_medians)
    calculate_cti_y(row_medians)

def save_image_detailed(images, ccd_image_name):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    axs = axs.flatten()

    for i, (image, channel, noise) in enumerate(images):
        ax = axs[i]
        im = ax.imshow(image, cmap='viridis', origin='lower', vmin=-1, vmax=80)
        ax.set_title(f'Combined Image [Channel {channel - 1}]', fontweight="bold")
        ax.set_xlabel("Column Number")
        ax.set_ylabel("Row Number")
        divider = make_axes_locatable(ax)
        axtop = divider.append_axes("top", size="20%", pad=0.7, sharex=ax)
        axright = divider.append_axes("right", size="10%", pad=0.5, sharey=ax)

        high, low, _ = find_col_traps(image)
        col_medians = np.median(image, axis=0)
        row_medians = np.median(image, axis=1)
        col_medians[:prescan] = np.nan

        axtop.step(np.arange(len(col_medians)), col_medians, where='mid')
        axtop.step(np.arange(len(high)), high, color='red', where='mid')
        axtop.step(np.arange(len(low)), low, color='red', where='mid')

        axright.step(row_medians, np.arange(len(row_medians)), where='mid')
        axtop.set_ylabel("[ADU]")
        axright.set_xlabel("[ADU]")
        axtop.set_title("Column-wise median", fontweight="bold")
        axright.set_title("Row-wise median", fontweight="bold")

        calculate_cti_x(col_medians)
        calculate_cti_y(row_medians)

    fig.suptitle(f'CCD {ccd_image_name}', fontweight="bold")
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
    save_image_detailed(images, ccd_image_name)

    log_file = Path(f"{ccd_image_name}.log")
    with log_file.open('w') as f:
        original_stdout = sys.stdout
        sys.stdout = f
        try:
            for img, ch, noise in images:
                show_image_detailed(img, ch, noise)
        finally:
            sys.stdout = original_stdout

    print(f"[DONE] Saved: {ccd_image_name}.png and {ccd_image_name}.log")
    plt.show()

if __name__ == "__main__":
    main()
