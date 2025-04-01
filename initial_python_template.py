# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 16:45:18 2025

@author: tracy
"""

# notes:
# 1) is it phi_2/2 or phi_2 in field strength from point return

import numpy as np
import matplotlib.pyplot as plt


# defined constants
GRAVITY_CONST = 6.67*10**(-11)
MASS_EARTH = 5.972*10**24

# steps
RHO_STEPS_CALC = 10
THETA_STEPS_CALC = 10
PHI_STEPS_CALC = 10
RHO_STEPS_PLOT = 10
THETA_STEPS_PLOT = 100
PHI_STEPS_PLOT = 100  # currently useless

# plot parameters
RHO_MAX_PLOT = 3*10**6
MAJOR = 5*10**6
MINOR = 10**6


def change_in_params_calc():

    delta_rho = MINOR / RHO_STEPS_CALC
    delta_theta = 2*np.pi / (THETA_STEPS_CALC+2)
    delta_phi = np.pi / (PHI_STEPS_CALC+2)

    return delta_rho, delta_theta, delta_phi


def change_in_params_plot():

    delta_rho = RHO_MAX_PLOT / RHO_STEPS_PLOT
    delta_theta = (2*np.pi) / THETA_STEPS_PLOT

    return delta_rho, delta_theta


def find_mass_element():

    return MASS_EARTH / (RHO_STEPS_CALC*THETA_STEPS_CALC*PHI_STEPS_CALC)


def distance(rho_1, rho_2, theta_1, theta_2, phi_2):

    plane_d = np.sqrt(rho_1**2 + rho_2**2 - 2*rho_1 *
                      rho_2*np.cos(theta_2 - theta_1))

    chord_d = 2*MAJOR*np.sin((-1*phi_2)/2)

    return np.sqrt(plane_d**2 + chord_d**2)


def field_strength_from_point(d, phi_2):

    mass_element = find_mass_element()

    if d == 0:

        return 0

    else:

        return ((GRAVITY_CONST*mass_element)/(d**2))*np.abs(np.cos(phi_2))


def field_strength_from_torus(rho_1, theta_1):

    mass_element = find_mass_element()
    delta_rho, delta_theta, delta_phi = change_in_params_calc()

    rho_2 = 0
    theta_2 = delta_theta
    phi_2 = delta_phi

    total_field_strength = 0

    while phi_2 < (np.pi - delta_phi):

        theta_2 = delta_theta

        while theta_2 < (2*np.pi - delta_theta):

            rho_2 = 0

            while rho_2 <= MINOR:

                total_field_strength += field_strength_from_point(
                    distance(rho_1, rho_2, theta_1, theta_2, phi_2), phi_2)

                rho_2 += delta_rho

            theta_2 += delta_theta

        phi_2 += delta_phi

    return 2*total_field_strength


def create_cross_section_array():

    cross_section_array = []

    delta_rho_plot, delta_theta_plot = change_in_params_plot()

    rho_1 = 0
    theta_1 = 0

    while theta_1 < 2*np.pi:

        rho_1 = 0

        while rho_1 < RHO_MAX_PLOT:

            field_strength = field_strength_from_torus(rho_1, theta_1)

            cross_section_array.append([rho_1, theta_1, field_strength])

            rho_1 += delta_rho_plot

        theta_1 += delta_theta_plot

    return cross_section_array


def main():

    print(create_cross_section_array())

    return 0


main()
