"""
File that stores all of the classes
Joseph Kent
30/04/2020
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from calculations import sphere_volume


class Planetessimal:
    def __init__(self, id_1, mass, displacement_vector, velocity_vector, composition):
        self.id = id_1
        self.radius = mass ** 0.5
        self.composition = {composition: mass}
        self.mass = mass
        self.displacement_vector = displacement_vector
        self.velocity_vector = velocity_vector
        self.acceleration = np.array([0.0, 0.0, 0.0])

    def __str__(self):
        return "Object {} with a mass {} at ({}, {}, {})".format(self.id, self.mass,
                                                                 round(self.displacement_vector[0]),
                                                                 round(self.displacement_vector[1]),
                                                                 round(self.displacement_vector[2]))

    def __repr__(self):
        return "planetesssimal({}, {}, {}, {}, {})".format(self.id, self.mass, self.displacement_vector,
                                                           self.velocity_vector, self.composition)

    def collison(self, other_object, chemical_properties):
        self.velocity_vector = (self.momentum() + other_object.momentum()) / (
                self.mass + other_object.mass)
        self.mass += other_object.mass
        for key, value in other_object.composition.items():
            self.composition[key] = self.composition.get(key, 0) + value
        self.calculate_radius(chemical_properties)
        other_object.mass = 0

    def calculate_acceleration(self, object_2, gravitational_distance, distance, gravitational_constant):
        if distance < gravitational_distance * object_2.mass ** 0.5 or distance < gravitational_distance * self.mass ** 0.5:
            direction_vector = object_2.displacement_vector - self.displacement_vector
            acceleration = (direction_vector / (
                    direction_vector ** 2).sum() ** 0.5) * gravitational_constant / (
                                   distance ** 2)
            self.acceleration += acceleration * object_2.mass
            object_2.acceleration -= acceleration * self.mass

    def calculate_displacement(self):
        """calculates the displacement and velocity vectors"""
        self.velocity_vector += self.acceleration
        self.displacement_vector += self.velocity_vector

    def momentum(self):
        """returns the objects momentum"""
        return self.mass * self.velocity_vector

    def shift_center(self, new_center):
        """shifts the objects center by an amount"""
        self.displacement_vector -= new_center

    def set_acceleration_to_zero(self):
        """sets the acceleration to zero"""
        self.acceleration = np.array([0.0, 0.0, 0.0])

    def calculate_radius(self, chemical_properties):
        volume = 0
        for key, value in self.composition.items():
            volume += value / chemical_properties[key][0]
        self.radius = (volume * 3 / (4 * np.pi)) ** (1 / 3)


class DensityCurve:
    def __init__(self, resource, amount, minimum, maximum, threshold):
        self.resource = resource
        self.amount = amount
        self.maximum = maximum
        self.minimum = minimum
        self.threshold = threshold
        self.points = {}

    def get_density(self, distance, mass):
        if mass > self.threshold:
            first_half = np.e ** (-self.minimum + distance)
            second_half = np.e ** (-self.maximum + distance)
            amount = self.amount * (first_half / (first_half + 1) - second_half / (second_half + 1)) * (
                        self.amount - self.points.get(distance // 10, 0))
            self.points[distance // 10] = self.points.get(distance // 10, 0) + amount
            return self.resource, amount
        else:
            return self.resource, 0


class OptimisableValues:
    def __init__(self, gravity, acceleration_error, spin, maximum_distance, simulations_before_calculation):
        self.gravity = gravity
        self.acceleration_error = acceleration_error
        self.spin = spin
        self.maximum_distance = maximum_distance
        self.simulations_before_calculation = simulations_before_calculation

    def return_variation(self):
        weights = np.random.uniform(0.6, 1.4, 5)
        return OptimisableValues(self.gravity * weights[0], round(self.acceleration_error * weights[1] + 1),
                                 self.spin * weights[2], self.maximum_distance * weights[3],
                                 self.simulations_before_calculation * weights[4])


class Paths:
    def __init__(self, id_list):
        self.dictionary_of_paths = {}
        for id_value in id_list:
            self.dictionary_of_paths[id_value] = Path()

    def display(self, filename, storage_address):
        """shows an image of all of the objects paths"""
        ax = plt.axes(projection='3d')
        ax.grid(False)
        for path in list(self.dictionary_of_paths.values()):
            x, y, z, radius = path.get_lists()
            ax.scatter3D(x, y, z,
                         s=radius)
        # plt.savefig(storage_address + filename + ".png")
        plt.show()


    def add(self, index, displacement_vector, radius):
        """adds a value to the list of positions"""
        self.dictionary_of_paths[index].add_value(displacement_vector[0], displacement_vector[1],
                                                  displacement_vector[2], radius)

class Path:
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.radius = []

    def add_value(self, x, y, z, radius):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.radius.append(radius)

    def get_lists(self):
        return self.x, self.y, self.z, self.radius


class InnerPlanet:
    def __iter__(self, radius, composition, temperature):
        self.layers = []
        self.radius = radius
        self.composition = composition
        self.volume = sphere_volume(radius)
        for i in range(radius // 1):
            self.layers.append(Layer(composition, (sphere_volume(i + 1) - sphere_volume(i)) / self.volume))


class Layer:
    def __init__(self, composition, percentage):
        self.composition = {}
        for key, value in list(composition.items()):
            self.composition[key] = value * percentage


class Chemicals:
    def __init__(self, name, density, probability, triple_point, critical_point, gradient):
        self.triple_point = triple_point
        self.critical_point = critical_point
        self.gradient = gradient
        self.name = name
        self.density = density
        self.probability = probability

    def display(self):
        """displays a graph of the state of matter"""
        x = np.array(range(0, self.triple_point[0]+1))
        y = eval("({}*x**2)/({}**2)".format(self.triple_point[1],self.triple_point[0]))
        plt.plot(x, y)
        x = np.array(range(self.triple_point[0], self.critical_point[0]+1))
        y = eval("(({})/({}))*(x-{})**2 + {}".format(self.triple_point[1]-self.critical_point[1], self.critical_point[0]-self.triple_point[0]**2, self.triple_point[0], self.triple_point[1]))
        plt.plot(x, y)
        y = np.array(range(self.triple_point[1], 308))
        x = eval("{}*(y-{})**3+{}".format(1/self.gradient,self.triple_point[1],self.triple_point[0]))
        plt.plot(x, y)
        plt.show()

    def chemical_type(self, pressure, temperature):
        if temperature < self.triple_point[0]:
            if (self.triple_point[1]*temperature**2)/(self.triple_point[0]**2) >= pressure:
                return "gas"
            elif pressure < self.triple_point[1]:
                return "solid"
        elif temperature < self.gradient*(pressure*self.triple_point[1])**3+self.triple_point[0] and pressure > self.triple_point[1]:
            return "solid"
        elif temperature < self.critical_point[0] and pressure <= ((self.triple_point[1]-self.critical_point[1])/(self.critical_point[0]-self.triple_point[0]**2))*(temperature-self.triple_point[0])**2 +self.triple_point[1]:
            return "gas"
        elif temperature >= self.critical_point[0]:
            if pressure <= self.critical_point[0]:
                return "gas"
            else:
                return "super_critical"
        else:
            return "liquid"