import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

if len(sys.argv) != 5:
    print("Usage: python script.py <file.dat> <title> <X_title> <Y_title>")
    sys.exit(0)

dat, title, xtitle, ytitle = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

print(f"Title: {title}")

# Load the data
data = pd.read_csv(dat, delimiter=' ', header=None)
n_values, values = data[0], data[1]

# Calculate min and max values
maxY, minY = values.max(), values.min()
maxX, minX = n_values.max(), n_values.min()

print(f"MaxY: {maxY}, Min: {minY}")
print(f"MaxX: {maxX}, Min: {minX}")

# Calculate scale
scaleY = (maxY - minY) / 10
scaleX = (maxX - minX) // 10

# Adjust max and min values
maxX += scaleX
minX -= scaleX
maxY += scaleY
minY -= scaleY

print(f"NewMaxY: {maxY}, NewMin: {minY}")
print(f"NewMaxX: {maxX}, NewMin: {minX}")

# Plot the main graph
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(n_values, values, color='red')
plt.xlabel(xtitle)
plt.ylabel(ytitle)
plt.title(title)
plt.axis([minX, maxX, minY, maxY])

# Calculate and plot density
kde = gaussian_kde(values)
dens_x = np.linspace(values.min(), values.max(), 1000)
dens_y = kde(dens_x)

# Normalize density
dens_y /= dens_y.sum()

# Save population data
population = dens_y * np.diff(np.concatenate(([dens_x[0]], dens_x)))

# Normalize the population values so their sum is 1
total_population = population.sum()
population /= total_population

# Create the population dataframe with the desired syntax
population_df = pd.DataFrame({'Index': np.arange(1, len(dens_x) + 1), 'Value': dens_x, 'Density': population})
population_df.to_csv(f"{dat}.population.csv", index=False, header=False, float_format='%.17g', quoting=1)

# Plot the density graph
plt.subplot(1, 2, 2)
plt.plot(population, dens_x, color='red')
plt.xlabel('Density')
plt.gca().invert_yaxis()

# Save the figure
plt.savefig(f"{dat}.png")
plt.show()

