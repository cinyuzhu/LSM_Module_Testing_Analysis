
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew


def passes_selection(cluster, selection):
    try:
        return eval(selection, {}, cluster)
    except:
        return False

def make_composite(listfile, selection, front=True, efact=0.0029):
    files = []
    with open(listfile, 'r') as f:
        for line in f:
            fname = line.strip()
            if fname:
                files.append(fname)

    nbins = 400
    bmin, bmax = -4, 4
    composite = np.zeros((nbins, nbins), dtype=np.float64)
    sum_skew_y = 0
    total_clusters = 0
    # Lists to store cluster-level variables for plotting
    energy_list = []
    npix_list = []
    dx_list = []
    dy_list = []
    rmsxy_list = []

    for fname in files:
        print(f"Processing {fname}")
        with uproot.open(fname) as f:
            tree = f["clustersRec"]
            data = tree.arrays([
                "pixels_x", "pixels_y", "pixels_E",
                "PosX", "PosY", "wSTD_X", "wSTD_Y",
                "Energy", "Npix", "wSTD_XY","has_seed"
            ])

            for i in range(len(data["pixels_x"])):
                for j in range(len(data["pixels_x"][i])):
                    cluster = {
                        "pixels_x": data["pixels_x"][i][j],
                        "pixels_y": data["pixels_y"][i][j],
                        "pixels_E": data["pixels_E"][i][j],
                        "x": data["PosX"][i][j],
                        "y": data["PosY"][i][j],
                        "dx": data["wSTD_X"][i][j],
                        "dy": data["wSTD_Y"][i][j],
                        "ene": data["Energy"][i][j],
                        "np": data["Npix"][i][j],
                        # "rmsxy": np.sqrt((data["wSTD_X"][i][j]**2 + data["wSTD_Y"][i][j]**2)/2),
                        "rmsxy": data["wSTD_XY"][i][j],
                        "valid": data["has_seed"][i][j] == 1,
                    }

                    if not passes_selection(cluster, selection):
                        continue

                    # Print basic info about the cluster
                    print(f"\nCluster {total_clusters + 1}:")
                    print(f"  Cluster Energy: {cluster['ene']}")
                    print(f"  N pixels: {len(data['pixels_E'][i][j])}")
                    print(f"  Sample pixel energies: {data['pixels_E'][i][j][:5]}")

                    x = np.array(cluster["pixels_x"])
                    y = np.array(cluster["pixels_y"])
                    e = np.array(cluster["pixels_E"])

                    if len(x) == 0:
                        continue

                    x_shift = np.round(x - cluster["x"]).astype(int) + nbins // 2
                    y_shift = np.round(y - cluster["y"]).astype(int) + nbins // 2

                    valid_mask = (x_shift >= 0) & (x_shift < nbins) & (y_shift >= 0) & (y_shift < nbins)
                    x_shift = x_shift[valid_mask]
                    y_shift = y_shift[valid_mask]
                    e = e[valid_mask]
                    # Print after masking
                    print(f"  Total charge after mask: {np.sum(e)}")
                   


                    # Border test
                    border_mask = ((x_shift < 50) | (x_shift > nbins - 50)) & ((y_shift < 50) | (y_shift > nbins - 50))
                    borq = np.sum(e[border_mask]) / 2500.0

                    print(f"  Border charge (borq): {borq}")

                    if np.abs(borq) < 25:
                        for xi, yi, ei in zip(x_shift, y_shift, e):
                            # if not (0 <= xi < nbins and 0 <= yi < nbins):
                            #     print(f"  WARNING: Shifted index out of bounds: xi={xi}, yi={yi}")
                            # else:
                            #     print(f"  Filling composite[{yi}, {xi}] with {ei / 2500.0}")
                            composite[yi, xi] += ei / 2500.0
                            energy_list.append(cluster["ene"])
                            npix_list.append(cluster["np"])
                            dx_list.append(cluster["dx"])
                            dy_list.append(cluster["dy"])
                            rmsxy_list.append(cluster["rmsxy"])
                        total_clusters += 1
                        sum_skew_y += cluster["y"] + cluster["dy"] if front else cluster["y"] - cluster["dy"]

    if total_clusters == 0:
        print("No clusters matched the selection.")
        return

    composite /= total_clusters

    print(f"Composite from {total_clusters} clusters")
    print(f"Raw sum of composite: {np.sum(composite)} keV")
    total_el = np.sum(composite)/ 3.8
    print(f"Total integral: {total_el:.4f} electrons")
    # print(f"Total integral: {np.sum(composite) * efact * 1000 / 3.8:.1f} electrons")
    print(f"Mean Skew Direction (Y): {sum_skew_y / total_clusters:.3f}")

    # Pixel center projections (central +/- 0.5)
    xaxis = np.linspace(bmin, bmax, nbins)
    yaxis = np.linspace(bmin, bmax, nbins)

    def get_range_idx(axis, lo, hi):
        return np.where((axis >= lo) & (axis < hi))[0]

    xc = get_range_idx(xaxis, -0.5, 0.5)
    yc = get_range_idx(yaxis, -0.5, 0.5)
    xup = get_range_idx(xaxis, -0.5, 0.5)
    yup = get_range_idx(yaxis, 0.5, 1.5)
    xdown = get_range_idx(xaxis, -0.5, 0.5)
    ydown = get_range_idx(yaxis, -1.5, -0.5)
    xleft = get_range_idx(xaxis, -1.5, -0.5)
    yleft = get_range_idx(yaxis, -0.5, 0.5)
    xright = get_range_idx(xaxis, 0.5, 1.5)
    yright = get_range_idx(yaxis, -0.5, 0.5)

    central = composite[np.ix_(yc, xc)].sum()
    above = composite[np.ix_(yup, xup)].sum()
    below = composite[np.ix_(ydown, xdown)].sum()
    left = composite[np.ix_(yleft, xleft)].sum()
    right = composite[np.ix_(yright, xright)].sum()

    print(f"Above fraction: {above / central:.3f}")
    print(f"Below fraction: {below / central:.3f}")
    print(f"Right fraction: {right / central:.3f}")
    print(f"Left fraction:  {left / central:.3f}")

    # Parameters
    nbins = composite.shape[0]
    bmin, bmax = -4, 4
    bin_edges = np.linspace(bmin, bmax, nbins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    # Projection limits: -4 to -3
    proj_range = (-4, -3)
    proj_idx = np.where((bin_centers >= proj_range[0]) & (bin_centers <= proj_range[1]))[0]

    # 1D Projections
    hbelow = np.sum(composite[proj_idx, :], axis=0)  # X projection (horizontal strip)
    hleft = np.sum(composite[:, proj_idx], axis=1)   # Y projection (vertical strip)

    # Normalize
    efact = 1  # Optional energy factor if needed
    hbelow /= hbelow.sum() if hbelow.sum() > 0 else 1
    hleft  /= hleft.sum() if hleft.sum() > 0 else 1

    # Stats
    print(f"\nProjection X:")
    print(f"  Integral: {np.sum(hbelow) * efact:.2f} electrons")
    print(f"  Mean:     {np.mean(hbelow):.6f}, RMS: {np.std(hbelow):.6f}, Skew: {skew(hbelow):.6f}")

    print(f"\nProjection Y:")
    print(f"  Integral: {np.sum(hleft) * efact:.2f} electrons")
    print(f"  Mean:     {np.mean(hleft):.6f}, RMS: {np.std(hleft):.6f}, Skew: {skew(hleft):.6f}")

    # ---- PLOTTING ----
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Composite image
    im = axs[0].imshow(composite, extent=[bmin, bmax, bmin, bmax], origin='lower', cmap='jet')
    axs[0].set_title("composite", fontsize=16, weight='bold')
    axs[0].axhspan(proj_range[0], proj_range[1], xmin=0, xmax=1, edgecolor='lime', facecolor='none', linestyle=':', linewidth=2)
    axs[0].axvspan(proj_range[0], proj_range[1], ymin=0, ymax=1, edgecolor='lime', facecolor='none', linestyle=':', linewidth=2)
    axs[0].text(-3.8, 3.2, "left", color='limegreen', rotation=90, fontsize=14)
    axs[0].text(-1.0, -3.8, "below", color='limegreen', fontsize=14)
    plt.colorbar(im, ax=axs[0], fraction=0.046, pad=0.04)

    # Projection plot
    axs[1].plot(bin_centers, hbelow, color='black', label='hbelow')
    axs[1].plot(bin_centers, hleft, color='red', label='hleft')
    axs[1].set_title("Projections")
    axs[1].legend()

    plt.tight_layout()
    plt.show()

    # Cluster-level variable distribution plots 
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    axs = axs.flatten()

    axs[0].hist(energy_list, bins=50, color='dodgerblue', alpha=0.7)
    axs[0].set_title("Cluster Energy [keV]")

    axs[1].hist(npix_list, bins=50, color='orange', alpha=0.7)
    axs[1].set_title("Number of Pixels per cluster")

    axs[2].hist(dx_list, bins=50, color='green', alpha=0.7)
    axs[2].set_title("Cluster sigma_x [pixel]")

    axs[3].hist(dy_list, bins=50, color='red', alpha=0.7)
    axs[3].set_title("Cluster sigma_y [pixel]")

    axs[4].hist(rmsxy_list, bins=50, color='purple', alpha=0.7)
    axs[4].set_title("Cluster sigma_xy [pixel]")

    # Hide empty subplot
    fig.delaxes(axs[5])

    plt.tight_layout()
    plt.show()

# === Example usage ===
if __name__ == "__main__":
    make_composite("test_list.txt", "valid", front=False)
    # make_composite("test_list.txt", "valid and ene > 4.5 and ene < 7 and rmsxy > 0.7", front=False)





