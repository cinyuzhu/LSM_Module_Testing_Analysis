import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from pathlib import Path
from scipy import stats

# Config
overscan_start = 6152
prescan = 8
e_transient_end = 100
running_window = 40
col_threshold = 6

def subtract_pedestal(data, log_lines, ext_index, filepath):
    log_lines.append(f"Data loaded {filepath} Extension {ext_index}")
    row_medians = np.median(data[:, overscan_start:], axis=1)
    subtracted = data - row_medians[:, np.newaxis]
    noise_os = np.std(subtracted[:, overscan_start:])
    mean_os = np.mean(subtracted[:, overscan_start:])
    threshold = mean_os + 5 * noise_os
    valid_pixels = subtracted[:, overscan_start:][subtracted[:, overscan_start:] < threshold]
    noise = round(np.std(valid_pixels), 1)
    log_lines.append(f"Noise [ADU] :  {noise}")
    return subtracted

def running_median_absolute_deviation(arr, window_size):
    medians, mads = [], []
    arr = np.where(np.ma.getmask(arr), np.nan, arr)
    for i in range(len(arr)):
        start = max(0, i - window_size // 2)
        end = min(len(arr), i + window_size // 2 + 1)
        medians.append(np.nanmedian(arr[start:end]))
        mads.append(stats.median_abs_deviation(arr[start:end], nan_policy='omit'))
    return np.array(medians), np.array(mads)

def find_col_traps(image, threshold=col_threshold):
    masked_image = np.ma.masked_array(image)
    masked_image[:, :prescan] = np.ma.masked
    masked_image[:, overscan_start - prescan:] = np.ma.masked

    col_medians = np.ma.median(masked_image, axis=0)

    combined_medians = np.empty_like(col_medians)
    combined_mads = np.empty_like(col_medians)

    run_medians, run_mads = running_median_absolute_deviation(
        col_medians[:e_transient_end], running_window)
    combined_medians[:e_transient_end] = run_medians
    combined_mads[:e_transient_end] = run_mads

    combined_medians[e_transient_end:overscan_start] = np.ma.median(col_medians[e_transient_end:])
    combined_mads[e_transient_end:overscan_start] = stats.median_abs_deviation(col_medians[e_transient_end:])

    condition_high = combined_medians + threshold * combined_mads
    condition_low = combined_medians - threshold * combined_mads

    traps = (col_medians[prescan:overscan_start - prescan] > condition_high[prescan:overscan_start - prescan]) | \
            (col_medians[prescan:overscan_start - prescan] < condition_low[prescan:overscan_start - prescan])

    trap_indices = np.where(traps)[0] + prescan
    return condition_high, condition_low, trap_indices

def process_fz_and_plot_all_channels(fz_file, ccd_name):
    log_lines = []
    with fits.open(fz_file) as hdul:
        fig, axs = plt.subplots(2, 2, figsize=(18, 12))
        axs = axs.flatten()

        for i in range(4):
            data = np.array(hdul[i + 1].data)
            image = subtract_pedestal(data, log_lines, i + 1, fz_file)

            ax = axs[i]
            im = ax.imshow(image, cmap='cividis', origin='lower', vmin=-1, vmax=50, aspect='auto')
            ax.set_title(f'{ccd_name} – Channel {i}', fontsize=14, fontweight="bold")
            ax.set_xlabel("Column Number", fontsize=12)
            ax.set_ylabel("Row Number", fontsize=12)
            ax.set_ylim(1, 18)

            divider = make_axes_locatable(ax)
            axtop = divider.append_axes("bottom", size="30%", pad=0.8, sharex=ax)
            axright = divider.append_axes("right", size="15%", pad=0.8, sharey=ax)
            cax = divider.append_axes("right", size="5%", pad=0.3)

            # Colorbar
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label("Signal [ADU]", fontsize=12)
            cbar.ax.tick_params(labelsize=10)

            # Column and row medians
            col_medians = np.ma.median(image, axis=0)
            col_medians[:prescan] = np.nan
            row_medians = np.ma.median(image, axis=1)

            condition_high, condition_low, trap_indices = find_col_traps(image)
            log_lines.append(f"Defect Locations:  {trap_indices.tolist()}")

            # Plot column median + thresholds
            axtop.step(np.arange(image.shape[1]), col_medians, where='mid', label='Median')
            axtop.step(np.arange(image.shape[1]), condition_high, color='red', linestyle='dashed', label='Median + 6 MAD')
            axtop.step(np.arange(image.shape[1]), condition_low, color='red', linestyle='dotted', label='Median - 6 MAD')
            axtop.set_ylabel("ADU", fontsize=10)
            axtop.set_xlabel("Column Number", fontsize=10)
            axtop.set_title("Column Median", fontsize=10)

            axtop.set_ylim(-10, 10)
            axtop.tick_params(labelsize=9)

            axtop.legend(
                loc='upper right',
                bbox_to_anchor=(1.02, 1.35),
                borderaxespad=0.,
                fontsize=8,
                frameon=True
            )
            axtop.grid(True, linestyle='--', alpha=0.3)

            # Row medians
            axright.step(row_medians, np.arange(image.shape[0]), where='mid')
            axright.set_xlabel("ADU", fontsize=10)
            axright.set_title("Row Median", fontsize=10)
            axright.tick_params(labelsize=9)

        fig.suptitle(f'Serial Register Defect Maps – {ccd_name}', fontsize=16, fontweight="bold")
        plt.tight_layout(rect=[0, 0, 0.96, 0.95])
        outname = f"{ccd_name}.png"
        plt.savefig(outname, dpi=300)
        plt.close()
        print(f"[SAVED] {outname}")

    # Write to log
    log_path = f"{ccd_name}.log"
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines))
    print(f"[LOG SAVED] {log_path}")

def main():
    args = sys.argv
    if len(args) != 3:
        print(f"Usage: {args[0]} input.fz CCD_LABEL")
        sys.exit(1)

    fz_file = Path(args[1])
    ccd_name = args[2]

    if not fz_file.exists():
        print(f"[ERROR] File not found: {fz_file}")
        sys.exit(1)

    process_fz_and_plot_all_channels(fz_file, ccd_name)
    print("[DONE] Combined plot and log generated.")

if __name__ == "__main__":
    main()
