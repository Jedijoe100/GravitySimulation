"""
Created by Joseph Kent
27/04/2020
This file holds the main calculation functions
There is a list bellow
calculate_collisions
calculate_distance
"""
import numpy as np


def calculate_collisions(planetessimals, graviational_distance, gravitational_constant, chemical_properties, density_regions):
    """Calculates the collision if both objects clip each other"""
    for particle in planetessimals:
        particle.set_acceleration_to_zero()
    for i in range(len(planetessimals)):
        for m in range(i + 1, len(planetessimals)):
            distance = calculate_distance(planetessimals[i].displacement_vector, planetessimals[m].displacement_vector)
            if distance <= planetessimals[i].radius + planetessimals[m].radius:
                planetessimals[i].collison(planetessimals[m], chemical_properties)
            else:
                planetessimals[i].calculate_acceleration(planetessimals[m], graviational_distance, distance, gravitational_constant)
                """for density_curve in density_regions:
                    chemical, amount = density_curve.get_density(calculate_distance(planetessimals[i].displacement_vector, np.array([0.0, 0.0, 0.0])), planetessimals[i].mass)
                    if amount > 10 ** -100:
                        planetessimals[i].composition[chemical] = planetessimals[i].composition.get(chemical, 0) + amount
                        planetessimals[i].velocity_vector = planetessimals[i].momentum() / (planetessimals[i].mass + amount)"""


def calculate_distance(vector1, vector2):
    """calculates the distance between the two vectors"""
    return np.linalg.norm(vector1 - vector2)

def sphere_volume(radius):
    return 4/3 * np.pi * radius
