"""
Joseph Kent
Date: 26/04/2020
Description: The main file for a project that generates a solar system with planetessimals
Coordinate system used will be (radius, longitude, latitude)

"""
import numpy as np
import time
from multiprocessing.dummy import Pool as ThreadPool
from object_file import Planetessimal, OptimisableValues, Paths
import calculations as cal
import image_processing as vis

NUMBER_OF_INITIAL_PLANETESSIMALS = 100
NUMBER_OF_SIMULATIONS = 1000000
COMPOSITION_CLOUD = {"water": (0.9, 0.6), "rock": (3, 0.4), "hydrogen": (0.09, 0)}
MASS_CONSTRAINTS = [1, 10]
DISPLACEMENT_CONSTRAINTS = 1000
STORAGE_ADDRESS = "sub_files/"

CONSTANT_VALUES = OptimisableValues(1, 0.0000000001, 0.1, 10, 1)


def main():
    # optimisation(16, 1000, CONSTANT_VALUES)
    seed = int(input("Please type a seed: "))
    run(CONSTANT_VALUES, seed, False)


def run(constant_values, seed=1, optimising=True):
    np.random.seed(seed)
    planetessimals = initialize(constant_values)
    start = time.process_time()
    paths, objects_cleaned = accretion_disk(planetessimals, constant_values, not (optimising), NUMBER_OF_SIMULATIONS)
    time_finished = time.process_time() - start
    if not optimising:
        vis.display_final(NUMBER_OF_INITIAL_PLANETESSIMALS, time_finished, planetessimals, objects_cleaned)
    return 1 / (objects_cleaned ** 2 * time_finished + time_finished)


def accretion_disk(planetessimals, constant_values, show_graphs, length_simulation=1000000):
    """simulates the accretion disk and returns a list of the objects with their
    positions, velocities, accelerations, mass and composition"""
    paths = Paths([particle.id for particle in planetessimals])
    objects_cleaned = 0
    gravitational_distance = (constant_values.gravity / constant_values.acceleration_error) ** 0.5
    length_paths = len(planetessimals)
    for i in range(length_simulation):
        for space_object in planetessimals:
            if show_graphs:
                paths.add(list(space_object.displacement_vector))
            space_object.calculate_displacement()
        for space_object in planetessimals:
            space_object.shift_center(planetessimals[0].displacement_vector)
        if i % constant_values.simulations_before_calculation == 0:
            cal.calculate_collisions(planetessimals, gravitational_distance, constant_values.gravity, COMPOSITION_CLOUD)
            objects_cleaned += clean_list(planetessimals, constant_values)
        if i % 50000 == 0 and show_graphs:
            print("{} number of simulations have occured or {}%".format(i, 100 * i / length_simulation))
            vis.process_path(paths,
                             "percentage_{}_for_{}_particles_and_{}_simulations".format(100 * i / length_simulation,
                                                                                        NUMBER_OF_INITIAL_PLANETESSIMALS,
                                                                                        NUMBER_OF_SIMULATIONS),
                             constant_values.maximum_distance, DISPLACEMENT_CONSTRAINTS, STORAGE_ADDRESS)
            paths.clear()
            for j in range(length_paths):
                paths.append([])

    return paths, objects_cleaned


def initialize(constant_values):
    planetessimals = [Planetessimal(id_1=0, mass=1000, displacement_vector=np.array([0.0, 0.0, 0.0]),
                                        velocity_vector=np.array([0.0, 0.0, 0.0]), composition="hydrogen")]
    for i in range(1, NUMBER_OF_INITIAL_PLANETESSIMALS):
        mass = np.random.randint(low=MASS_CONSTRAINTS[0], high=MASS_CONSTRAINTS[1])
        displacement_vector = np.array(
            [np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS),
             np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS),
             np.random.uniform(-DISPLACEMENT_CONSTRAINTS, DISPLACEMENT_CONSTRAINTS)])
        angle = np.random.uniform(0, 2 * np.pi)
        velocity_vector = constant_values.spin * np.random.uniform(0, 1) * np.array([np.cos(angle), np.sin(angle), 0])
        composition_element = np.random.choice([key for key in COMPOSITION_CLOUD.keys()], 1, p=[substance[1] for substance in COMPOSITION_CLOUD.values()])
        planetessimals.append(Planetessimal(id_1=i, mass=mass, displacement_vector=displacement_vector,
                                                velocity_vector=velocity_vector, composition=str(composition_element[0])))
    return planetessimals


def convert(composition, mass):
    return {composition.name: (composition.density, mass)}


def clean_list(planetessimals, constant_values):  # needs more work
    """returns the array with all objects that have mass 0 or have a distance larger than 10*that of the initial size"""
    objects_removed = 0
    for space_object in planetessimals:
        distance = cal.calculate_distance(space_object.displacement_vector, np.array([0, 0, 0]))
        if space_object.mass == 0 or DISPLACEMENT_CONSTRAINTS * constant_values.maximum_distance < distance:
            planetessimals.remove(space_object)
            if DISPLACEMENT_CONSTRAINTS * constant_values.maximum_distance < distance:
                objects_removed += 1
    return objects_removed


def optimisation(pop_size, optimisation_number, constant_values):
    """an optimisation algorithm that optimises for the best program performance"""
    i = 0
    best_pop = []
    pop_values = []
    for pop in range(pop_size):
        pop_values.append(constant_values.return_variation())
    while i < optimisation_number:
        pool = ThreadPool()
        results = pool.map(run, pop_values)
        pool.close()
        pool.join()
        maximum_value = max(results)
        index_value = results.index(maximum_value)
        best_pop.clear()
        best_pop.append(pop_values[index_value])
        for pop in pop_values:
            pop = best_pop[0].return_variation()
        print("simulation {} with a maximum value {}".format(i, maximum_value))
        i += 1
    print(
        "The best simulation had gravity {}, spin {}, maximum distance {}, acceleration error {}, simulations before calculation".format(
            best_pop[0].gravity, best_pop[0].spin, best_pop[0].maximum_distance, best_pop[0].acceleration_error,
            best_pop[0].simulations_before_calculation))


main()
