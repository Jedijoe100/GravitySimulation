"""
Joseph Kent
Date: 26/04/2020
Description: The main file for a project that generates a solar system with planetessimals
Coordinate system used will be (radius, longitude, latitude)

"""
import numpy as np
import time
from object_file import Planetessimal, OptimisableValues, Paths, DensityCurve
import calculations as cal
import image_processing as vis

NUMBER_OF_INITIAL_PLANETESSIMALS = 100
NUMBER_OF_SIMULATIONS = 1000000
COMPOSITION_CLOUD = {"water": (0.9, 0.6), "rock": (3.0, 0.4), "hydrogen": (0.09, 0)}
MASS_CONSTRAINTS = [1.0, 10.0]
DISPLACEMENT_CONSTRAINTS = 1000
STORAGE_ADDRESS = "sub_files/"
CONSTANT_VALUES = OptimisableValues(
    gravity=1,
    acceleration_error=0.0000001,
    spin=0.1,
    maximum_distance=10,
    simulations_before_calculation=1
)
DENSITY_REGIONS = [
    DensityCurve("hydrogen", 0.1, 0, 10000, 1000),
    DensityCurve("rock", 0.1, 2000, 6000, 0),
    DensityCurve("water", 0.1, 4000, 10000, 0)
]


class Simulation:
    def __init__(self, constant_values, seed=1, optimising=True):
        np.random.seed(seed)
        self.constant_values = constant_values
        self.optimising = optimising
        self.planetessimals = []
        self.objects_cleaned = 0
        self.gravitational_distance = (self.constant_values.gravity / self.constant_values.acceleration_error) ** 0.5

    def run(self):
        self.initialize()
        start = time.process_time()
        self.accretion_disk()
        time_finished = time.process_time() - start
        if not self.optimising:
            vis.display_final(NUMBER_OF_INITIAL_PLANETESSIMALS, time_finished, self.planetessimals,
                              self.objects_cleaned)
        return 1 / (self.objects_cleaned ** 2 * time_finished + time_finished)

    def initialize(self):
        self.planetessimals.append(Planetessimal(id_1=0, mass=1000, displacement_vector=np.array([0.0, 0.0, 0.0]),
                                                 velocity_vector=np.array([0.0, 0.0, 0.0]), composition="hydrogen"))
        for i in range(1, NUMBER_OF_INITIAL_PLANETESSIMALS):
            mass = np.random.randint(low=MASS_CONSTRAINTS[0], high=MASS_CONSTRAINTS[1])
            displacement_vector = np.array(
                [np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS),
                 np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS),
                 np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS)])
            angle = np.random.uniform(0, 2 * np.pi)
            velocity_vector = self.constant_values.spin * np.random.uniform(0, 1) * np.array(
                [np.cos(angle), np.sin(angle), 0])
            composition_element = np.random.choice([key for key in COMPOSITION_CLOUD.keys()], 1,
                                                   p=[substance[1] for substance in COMPOSITION_CLOUD.values()])
            self.planetessimals.append(Planetessimal(id_1=i, mass=mass, displacement_vector=displacement_vector,
                                                     velocity_vector=velocity_vector,
                                                     composition=str(composition_element[0])))

    def accretion_disk(self):
        """simulates the accretion disk and returns a list of the objects with their
        positions, velocities, accelerations, mass and composition"""
        paths = Paths([particle.id for particle in self.planetessimals])
        for i in range(NUMBER_OF_SIMULATIONS):
            paths = self.calculate_displacements(paths)
            self.shift_center()
            self.collision_calculations(i)
            if i % 50000 == 0 and not self.optimising:
                print("{} number of simulations have occured or {}%".format(i, 100 * i / NUMBER_OF_SIMULATIONS))
                paths.display(
                    "percentage_{}_for_{}_particles_and_{}_simulations".format(100 * i / NUMBER_OF_SIMULATIONS,
                                                                               NUMBER_OF_INITIAL_PLANETESSIMALS,
                                                                               NUMBER_OF_SIMULATIONS),
                    STORAGE_ADDRESS)
                paths = Paths([particle.id for particle in self.planetessimals])

    def shift_center(self):
        center = self.planetessimals[0].displacement_vector
        for space_object in self.planetessimals:
            space_object.shift_center(center)

    def calculate_displacements(self, paths):
        for space_object in self.planetessimals:
            if not self.optimising:
                paths.add(space_object.id, list(space_object.displacement_vector), space_object.radius)
            space_object.calculate_displacement()
        return paths

    def collision_calculations(self, i):
        if i % self.constant_values.simulations_before_calculation == 0:
            cal.calculate_collisions(self.planetessimals, self.gravitational_distance, self.constant_values.gravity,
                                     COMPOSITION_CLOUD, DENSITY_REGIONS)
            self.objects_cleaned += clean_list(self.planetessimals, self.constant_values)

def convert(composition, mass):
    return {composition.name: (composition.density, mass)}


def clean_list(planetessimals, constant_values):
    """returns the array with all objects that have mass 0 or have a distance larger than 10*that of the initial size"""
    objects_removed = 0
    for space_object in planetessimals:
        distance = cal.calculate_distance(space_object.displacement_vector, np.array([0, 0, 0]))
        if space_object.mass == 0 or DISPLACEMENT_CONSTRAINTS * constant_values.maximum_distance < distance:
            planetessimals.remove(space_object)
            if DISPLACEMENT_CONSTRAINTS * constant_values.maximum_distance < distance:
                objects_removed += 1
    return objects_removed

if __name__ == '__main__':
    # optimisation(16, 1000, CONSTANT_VALUES)
    seed = int(input("Please type a seed: "))
    simulation = Simulation(CONSTANT_VALUES, seed, False)
    simulation.run()