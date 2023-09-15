"""
Processes all of the visual aspects of the system
Joseph Kent
30/04/2020
"""


def display_final(number_of_initial_particles, time_finished, planetessimals, objects_cleaned):
    """displays the final objects"""
    print("Total time ellapsed to do {} simulations: {}".format(number_of_initial_particles, time_finished))
    print("With an initial {} planetessimals, there are now {}, {} were removed because they were to far".format(
        number_of_initial_particles,
        len(planetessimals), objects_cleaned))
    for space_object in planetessimals:
        print(space_object)