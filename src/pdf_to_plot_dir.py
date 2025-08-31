# usecase: wildcard pattern search, take the first 4 matches
# python plot_pdfs.py --files "trace_*.pdf" outputname
#
# compatible to the old version: take the first 4 files in a directory ()
# Example: python plot_pdfs.py --directory ./pdf_outputs outputname


import pdfplumber
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
import argparse
from natsort import natsorted
import glob
import os


def plot_pdf(pdf_paths, out_name, render_dpi, save_dpi):
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))
    axs = axs.flatten()

    for i, pdf_path in enumerate(pdf_paths[:4]):  # Only take first 4
        ax = axs[i]
        with pdfplumber.open(pdf_path) as pdf:
            page = pdf.pages[0]
            # ↑ Render the page to a raster at a higher DPI
            im = page.to_image(resolution=render_dpi)
            image_np = np.array(im.original)  # PIL Image -> numpy

        # Avoid Matplotlib smoothing the raster
        ax.imshow(image_np, cmap='gray', interpolation='nearest', aspect='equal')
        ax.axis('off')
        ax.set_title(Path(pdf_path).name, fontweight='bold', fontsize=10)

    plt.tight_layout()
    # save_dpi controls the *figure* DPI; render_dpi controls the source raster DPI
    plt.savefig(f'{out_name}.png', dpi=save_dpi, bbox_inches='tight', pad_inches=0.1)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Process and plot PDF files.")
    parser.add_argument('--files', type=str, help="Wildcard pattern for PDF files (e.g., 'trace_*.pdf').")
    parser.add_argument('--directory', type=str, help="Directory containing PDF files.")
    parser.add_argument('--render-dpi', type=int, default=500, help="DPI to rasterize PDF pages.")
    parser.add_argument('--save-dpi', type=int, default=500, help="DPI for the output PNG.")
    parser.add_argument('out_name', type=str, help="Name of the CCD for output image.")
    parser.add_argument('--cleanup', action='store_true', help="Delete the 4 PDF files after processing.")

    args = parser.parse_args()

    if not (args.files or args.directory):
        sys.exit("Provide --files <pattern> or --directory <dir>.")
    
    pdf_paths = []
    if args.files:
        # Expand wildcard matches
        pdf_paths = natsorted(glob.glob(args.files))
        if len(pdf_paths) < 4:
            sys.exit("The pattern must match at least 4 PDF files.")
    elif args.directory:
        pdf_dir = Path(args.directory)
        if not pdf_dir.exists() or not pdf_dir.is_dir():
            sys.exit(f"Directory {pdf_dir} does not exist or is not a directory.")
        pdf_paths = natsorted([str(p) for p in pdf_dir.glob("*.pdf")])
        if len(pdf_paths) < 4:
            sys.exit("The directory must contain at least 4 PDF files.")
    else:
        sys.exit("You must provide either a wildcard pattern with --files or a directory with PDFs.")

    # Take only first 4 matches
    pdf_paths = pdf_paths[:4]
    print("Using files:\n  " + "\n  ".join(pdf_paths))


    # Plot the PDFs
    plot_pdf(pdf_paths, args.out_name, args.render_dpi, args.save_dpi)

    if args.cleanup:
        for f in pdf_paths:
            try:
                os.remove(f)
                print(f"Deleted {f}")
            except Exception as e:
                print(f"Could not delete {f}: {e}")

if __name__ == "__main__":
    main()
