# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 13:08:56 2023

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'font.size': 15,
    'figure.figsize': (8,6),
    'savefig.dpi': 300,
    'grid.alpha': 0.75,
    'figure.constrained_layout.use':True
})
plt.style.use('seaborn-dark')


def stable_optical_cavity_mode(mirror_radius1_m, mirror_radius2_m, cavity_length_m):
    """Calculates the stability parameter g1g2 for an optical cavity.
    
    Args:
    mirror_radius1_m (float): Radius of curvature of the first cavity mirror in meters.
    mirror_radius2_m (float): Radius of curvature of the second cavity mirror in meters.
    cavity_length_m (float): Length of the cavity in meters.
    
    Returns:
    float: The stability parameter g1g2.
    """
    g1 = 1 - cavity_length_m/mirror_radius1_m
    g2 = 1 - cavity_length_m/mirror_radius2_m
    
    return g1*g2

def plot_optical_cavity_stability(mirror_radius1_m, mirror_radius2_m, cavity_length_m):
    """Plots the stability parameter g1g2 for an optical cavity as a function of the second mirror radius.
    
    Args:
    mirror_radius1_m (float): Radius of curvature of the first cavity mirror in meters.
    mirror_radius2_m (array): Array of radii of curvature of the second cavity mirror in meters.
    cavity_length_m (float): Length of the cavity in meters.
    
    Returns:
    the plot
    """

    plt.plot(mirror_radius2_m, stable_optical_cavity_mode(mirror_radius1_m, mirror_radius2_m, cavity_length_m), '.-')
    plt.hlines(1, mirror_radius2_m[0], mirror_radius2_m[-1], 'r', linestyle='dashed')
    plt.hlines(0, mirror_radius2_m[0], mirror_radius2_m[-1], 'r', linestyle='dashed', label='region of stability\n(existence)')
    plt.legend()
    plt.ylim([-0.5,1.5])
    plt.xlabel('Radius of curvature sample [m]')
    plt.ylabel('$g_1 g_2$')
    plt.savefig('g1g2.png')

if __name__ == '__main__':
    
    mirror_radius1_m = 1
    mirror_radius2_m = np.linspace(-5,5,1000)
    cavity_length_m = 0.1866
    tolerance = 0.0025

    plot_optical_cavity_stability(mirror_radius1_m, mirror_radius2_m, cavity_length_m)
    
    # Find the closest values of g1*g2 to 0 and 1
    g1g2 = stable_optical_cavity_mode(mirror_radius1_m, mirror_radius2_m, cavity_length_m)
    closest_to_0 = mirror_radius2_m[np.argmin(np.abs(g1g2))]
    closest_to_1 = mirror_radius2_m[np.argmin(np.abs(g1g2 - 1))]

    print(f"Closest value of g1*g2 to 0 is {stable_optical_cavity_mode(mirror_radius1_m, closest_to_0, cavity_length_m):.4f} at R = {closest_to_0:.4f} m")
    print(f"Closest value of g1*g2 to 1 is {stable_optical_cavity_mode(mirror_radius1_m, closest_to_1, cavity_length_m):.4f} at R = {closest_to_1:.4f} m")
    
