#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>

// Constants
const double MAJOR = 5e6;        // Major radius
const double MINOR = 1e6;        // Minor radius
const double GRAVITY_CONST = 6.67e-11; // Gravitational constant
const double MASS_EARTH = 5.972e24;  // Mass of Earth
const int NUM_RHO_STEPS = 10;
const int NUM_THETA_STEPS = 10;
const int NUM_PHI_STEPS = 10;

// Function to define the mass of a step
double define_mass_of_step() {
    int total_num_of_steps = NUM_RHO_STEPS * NUM_THETA_STEPS * NUM_PHI_STEPS;
    return MASS_EARTH / total_num_of_steps;
}

// Function to calculate distance between two points in toroidal coordinates
double distance(double rho_1, double rho_2, double theta_1, double theta_2, double phi_1, double phi_2) {
    double plane_d = std::sqrt(rho_1 * rho_1 + rho_2 * rho_2 - 2 * rho_1 * rho_2 * std::cos(theta_2 - theta_1));
    double chord_d = 2 * MAJOR * std::sin((phi_1 - phi_2) / 2);
    return std::sqrt(plane_d * plane_d + chord_d * chord_d);
}

// Function to calculate the field strength at a point
double field_strength_from_point(double d, double phi_2, double mass_of_step) {
    if (d == 0) {
        return 0;
    } else {
        return (GRAVITY_CONST * mass_of_step) / (d * d) * std::cos(phi_2 / 2);
    }
}

// Function to calculate field strength from the whole donut
double field_strength_from_whole_donut_at_given_point(double rho_1, double theta_1, double phi_1, double mass_of_step) {
    double rho_2 = 0;
    double theta_2 = 0;
    double phi_2 = 0;
    double total_potential = 0;
    double rho_step = MINOR / NUM_RHO_STEPS;
    double theta_step = 360.0 / NUM_THETA_STEPS;
    double phi_step = 180.0 / NUM_PHI_STEPS;

    while (phi_2 < 180) {
        theta_2 = 0;
        while (theta_2 < 360) {
            rho_2 = 0;
            while (rho_2 < MINOR) {
                total_potential += field_strength_from_point(
                    distance(rho_1, rho_2, theta_1, theta_2, phi_1, phi_2), phi_2, mass_of_step);
                rho_2 += rho_step;
            }
            theta_2 += theta_step;
        }
        phi_2 += phi_step;
    }

    return total_potential * 2;
}

// Function to convert toroidal coordinates to Cartesian coordinates
void toroidal_to_cartesian(double rho, double theta, double phi, double &x, double &y, double &z, double R = MAJOR) {
    x = (R + rho * std::cos(std::radians(theta))) * std::cos(std::radians(phi));
    y = (R + rho * std::cos(std::radians(theta))) * std::sin(std::radians(phi));
    z = rho * std::sin(std::radians(theta));
}

int main() {
    // Mass of a single step
    double step_mass = define_mass_of_step();
    
    double rho_1 = 0;
    double theta_1 = 0;
    double phi_1 = 0;
    double rho_step = MINOR / NUM_RHO_STEPS;
    double theta_step = 360.0 / NUM_THETA_STEPS;
    double phi_step = 360.0 / NUM_PHI_STEPS;
    
    std::vector<std::vector<double>> field_data;

    while (phi_1 < 360) {
        theta_1 = 0;
        while (theta_1 < 360) {
            rho_1 = 0;
            while (rho_1 < 1e7) {
                double field_strength = field_strength_from_whole_donut_at_given_point(rho_1, theta_1, phi_1, step_mass);
                field_data.push_back({rho_1, theta_1, phi_1, field_strength});
                rho_1 += rho_step * 10;
            }
            theta_1 += theta_step;
        }
        phi_1 += phi_step;
    }

    // Now, field_data holds the [rho, theta, phi, field_strength] for each point

    // Prepare data for plotting
    std::vector<double> x_data, y_data, z_data, field_strength_data;
    
    for (auto& entry : field_data) {
        double rho = entry[0], theta = entry[1], phi = entry[2], field_strength = entry[3];
        double x, y, z;
        toroidal_to_cartesian(rho, theta, phi, x, y, z);
        x_data.push_back(x);
        y_data.push_back(y);
        z_data.push_back(z);
        field_strength_data.push_back(field_strength);
    }

    // Save to CSV file to plot in external tool
    std::ofstream file("field_data.csv");
    for (size_t i = 0; i < x_data.size(); ++i) {
        file << x_data[i] << "," << y_data[i] << "," << z_data[i] << "," << field_strength_data[i] << "\n";
    }
    file.close();

    std::cout << "Data saved to field_data.csv." << std::endl;

    return 0;
}
