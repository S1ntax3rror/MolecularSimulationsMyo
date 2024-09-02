import matplotlib.pyplot as plt
from PIL import Image

# Load images
img1 = Image.open('CGenFFdynamics.png')
img2 = Image.open('CGenFFPocketIndex.png')
img3 = Image.open('MorseMDCMdynamics.png')
img4 = Image.open('MorseMDCMpocketIndex.png')

# Create a figure and a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(16, 9))

# Display images in the grid
axs[1, 0].imshow(img1)
axs[0, 0].axis('off')  # Hide axis

axs[0, 0].imshow(img2)
axs[0, 1].axis('off')  # Hide axis

axs[1, 1].imshow(img3)
axs[1, 0].axis('off')  # Hide axis

axs[0, 1].imshow(img4)
axs[1, 1].axis('off')  # Hide axis

# Adjust layout
plt.tight_layout()

# Save or show the result
plt.savefig('combined_image.png', dpi=300)
plt.show()