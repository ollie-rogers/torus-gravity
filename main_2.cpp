#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

const double GRAVITY_CONST = 6.67e-11;
const double MASS_EARTH = 5.972e24;


// steps
const int RHO_STEPS_CALC = 10;
const int THETA_STEPS_CALC = 10;
const int PHI_STEPS_CALC = 10;
const int RHO_STEPS_PLOT = 10;
const int THETA_STEPS_PLOT = 100;
const int PHI_STEPS_PLOT = 100;  // currently useless

// plot parameters
const double RHO_MAX_PLOT = 6e6;
const double MAJOR = 5e6;
const double MINOR = 1e6;

struct DataPoint {
    double rho_1;
    double theta_1;
    double field_strength;
};

// Function to calculate delta parameters for calculation
std::tuple<double, double, double> change_in_params_calc() {
    double delta_rho = MINOR / RHO_STEPS_CALC;
    double delta_theta = 2 * M_PI / (THETA_STEPS_CALC + 2);
    double delta_phi = M_PI / (PHI_STEPS_CALC + 2);
    return std::make_tuple(delta_rho, delta_theta, delta_phi);
}

// Function to calculate delta parameters for plotting
std::tuple<double, double> change_in_params_plot() {
    double delta_rho = RHO_MAX_PLOT / RHO_STEPS_PLOT;
    double delta_theta = (2 * M_PI) / THETA_STEPS_PLOT;
    return std::make_tuple(delta_rho, delta_theta);
}

// Function to find the mass of an element
double find_mass_element() {
    return MASS_EARTH / (RHO_STEPS_CALC * THETA_STEPS_CALC * PHI_STEPS_CALC);
}

// Function to calculate distance
double distance(double rho_1, double rho_2, double theta_1, double theta_2, double phi_2) {
    double plane_d = std::sqrt(rho_1 * rho_1 + rho_2 * rho_2 - 2 * rho_1 * rho_2 * std::cos(theta_2 - theta_1));
    double chord_d = 2 * MAJOR * std::sin((-1 * phi_2) / 2);
    return std::sqrt(plane_d * plane_d + chord_d * chord_d);
}

// Function to calculate field strength from a point
double field_strength_from_point(double d, double phi_2) {
    double mass_element = find_mass_element();
    if (d == 0) {
        return 0;
    } else {
        return (GRAVITY_CONST * mass_element) / (d * d) * std::cos(phi_2);
    }
}

// Function to calculate the total field strength from a torus
double field_strength_from_torus(double rho_1, double theta_1) {
    double mass_element = find_mass_element();
    auto [delta_rho, delta_theta, delta_phi] = change_in_params_calc();

    double rho_2 = 0;
    double theta_2 = delta_theta;
    double phi_2 = delta_phi;

    double total_field_strength = 0;

    while (phi_2 < (M_PI - delta_phi)) {
        theta_2 = delta_theta;

        while (theta_2 < (2 * M_PI - delta_theta)) {
            rho_2 = 0;

            while (rho_2 <= MINOR) {
                total_field_strength += field_strength_from_point(distance(rho_1, rho_2, theta_1, theta_2, phi_2), phi_2);
                rho_2 += delta_rho;
            }

            theta_2 += delta_theta;
        }

        phi_2 += delta_phi;
    }

    return 2 * total_field_strength;
}

// Function to create the cross section array
std::vector<DataPoint> create_cross_section_array() {
    std::vector<DataPoint> cross_section_array;

    auto [delta_rho_plot, delta_theta_plot] = change_in_params_plot();

    double rho_1 = 0;
    double theta_1 = 0;

    while (theta_1 < 2 * M_PI) {
        rho_1 = 0;

        while (rho_1 <= RHO_MAX_PLOT) {
            double field_strength = field_strength_from_torus(rho_1, theta_1);
            cross_section_array.push_back({rho_1, theta_1, field_strength});
            rho_1 += delta_rho_plot;
        }

        theta_1 += delta_theta_plot;
    }

    return cross_section_array;
}

// Function to save the data to a file
void save_data_to_file(const std::vector<DataPoint>& data) {
    std::ofstream outfile("cross_section_data.csv");

    // Writing headers
    outfile << "rho_1,theta_1,field_strength\n";

    // Writing data
    for (const auto& point : data) {
        outfile << point.rho_1 << "," << point.theta_1 << "," << point.field_strength << "\n";
    }

    outfile.close();
}

int main() {
    std::vector<DataPoint> data = create_cross_section_array();
    save_data_to_file(data);
    std::cout << "Data saved to cross_section_data.csv\n";
    return 0;
}
