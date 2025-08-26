import matplotlib.pyplot as plt
from PIL import Image
import os

# Specify the folder containing the PNG files and the output file path
folder_path = "./lowtemp_image4_figs/"  # Replace with your folder path
output_file = "composite_figure.png"  # Output file name

# Define the file naming pattern
extensions = [1, 2, 3, 4]
types = [
    "Composite_cluster__back",
    "Composite_cluster__front",
    "ctifrontevents",
    "energyspectrum",
]

# Organize filenames into a 4x4 grid
file_grid = [
    [f"{file_type}_dm01_image4_ext{ext}.png" for file_type in types]
    for ext in extensions
]

# Load images and validate they exist
images = []
for row in file_grid:
    img_row = []
    for file_name in row:
        file_path = os.path.join(folder_path, file_name)
        if os.path.exists(file_path):
            img_row.append(Image.open(file_path))
        else:
            raise FileNotFoundError(f"File not found: {file_path}")
    images.append(img_row)

# Get image dimensions
image_width, image_height = images[0][0].size

# Create the composite figure
dpi = 300  # High resolution
figure_width = (image_width * 4) / dpi
figure_height = (image_height * 4) / dpi
fig, axes = plt.subplots(4, 4, figsize=(figure_width, figure_height), dpi=dpi)

# Add images to the grid
for i, row in enumerate(images):
    for j, img in enumerate(row):
        axes[i, j].imshow(img)
        axes[i, j].axis("off")  # Hide axes

# Adjust layout and save the figure
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
plt.close()

print(f"Composite figure saved as {output_file}")

