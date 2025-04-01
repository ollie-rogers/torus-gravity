import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Step 1: Read the CSV file into a DataFrame (pandas will automatically treat the first row as header)


def load_data(file_path):
    # Read the CSV with the first row as the header
    data = pd.read_csv(file_path)

    # Convert columns to numeric (force errors to NaN if conversion fails)
    data["rho_1"] = pd.to_numeric(data["rho_1"], errors='coerce')
    data["theta_1"] = pd.to_numeric(data["theta_1"], errors='coerce')
    data["field_strength"] = pd.to_numeric(
        data["field_strength"], errors='coerce')

    # Drop any rows with NaN values that resulted from non-numeric data
    data = data.dropna()

    return data

# Step 2: Generate the contour plot


def create_contour_plot(data):
    # Extracting values
    rho = data['rho_1'].values  # rho_1
    theta = data['theta_1'].values  # theta_1
    field_strength = data['field_strength'].values  # field_strength

    # Convert polar (rho, theta) to Cartesian (x, y)
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)

    # Step 3: Create a grid for x and y (meshgrid)
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    # Define the resolution of the grid (increase for more detail)
    # Generate a 2D grid of x and y
    grid_x, grid_y = np.meshgrid(np.linspace(x_min, x_max, 200),
                                 np.linspace(y_min, y_max, 200))

    # Step 4: Interpolate field strength values onto the grid
    grid_field_strength = griddata(
        (x, y), field_strength, (grid_x, grid_y), method='cubic')

    # Step 5: Plot the contour
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(grid_x, grid_y, grid_field_strength,
                           levels=20, cmap='viridis')  # You can adjust levels
    plt.colorbar(contour)
    plt.title('Contour Plot of Field Strength (Cartesian)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


def main():
    # Path to your CSV file containing the data
    file_path = 'cross_section_data.csv'  # Update this to the actual file path

    # Load the data
    data = load_data(file_path)

    # Create the contour plot
    create_contour_plot(data)


if __name__ == "__main__":
    main()
