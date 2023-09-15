import numpy as np
import matplotlib.pyplot as plt
from calculations import sphere_volume

MAX_DISPLAY = 20
TEST_X = [0.1, 0.1, 2, 6, 6, 4, 4]
TEST_Y = [2, 5, 2, 4, 15, 15, 30]


def main():
    table_of_chemicals = TableOfChemicals()
    table_of_chemicals.append(Chemicals("water", 1, 0.1, [273, 611], [647, 22064], 3001, (0, 0, 1)))
    table_of_chemicals.append(Chemicals("carbon_dioxide", 1, 0.1, [1, 3], [5, 2], 3001, (0.3, 0.3, 1)))
    table_of_chemicals.append(Chemicals("rock", 1, 0.7, [1, 3], [5, 2], 3001, (0.3, 0.3, 0.3)))
    table_of_chemicals.append(Chemicals("nitrogen", 1, 0.1, [1, 3], [5, 2], 3001, (0.3, 0.3, 0.3)))
    table_of_chemicals["water"].display()
    planet = InnerPlanet(100, table_of_chemicals.return_random_chemical(3, 4), 100)
    for i in range(len(TEST_X)):
        print(TEST_X[i], TEST_Y[i], table_of_chemicals["water"].chemical_type(TEST_X[i], TEST_Y[i]))
    planet.display(table_of_chemicals)


class TableOfChemicals:
    def __init__(self):
        self.chemicals = {}

    def __getitem__(self, item):
        return self.chemicals[item]

    def append(self, chemical):
        self.chemicals[chemical.name] = chemical

    def return_random_chemical(self, number, mass):
        chemical_name = []
        chemical_probabilities = []
        for value in self.chemicals.values():
            chemical_name.append(value.name)
            chemical_probabilities.append(value.probability)
        chemical_names = np.random.choice(chemical_name, number, p=chemical_probabilities)
        chemicals_returning = {}
        for name in chemical_names:
            chemicals_returning[name] = chemicals_returning.get(name, 0) + mass / number
        return chemicals_returning


class Chemicals:
    def __init__(self, name, density, probability, triple_point, critical_point, gradient, colour):
        self.triple_point = triple_point
        self.critical_point = critical_point
        self.gradient = gradient
        self.name = name
        self.density = density
        self.probability = probability
        self.colour = colour

    def display(self):
        """displays a graph of the state of matter"""
        x = np.array(range(0, self.triple_point[0] + 1))
        y = eval("({}*x**2)/({}**2)".format(self.triple_point[1], self.triple_point[0]))
        plt.plot(x, y)
        x = np.array(range(self.triple_point[0], self.critical_point[0] + 1))
        y = eval("(({})/({}))*(x-{})**2 + {}".format(self.triple_point[1] - self.critical_point[1],
                                                     self.critical_point[0] - self.triple_point[0] ** 2,
                                                     self.triple_point[0], self.triple_point[1]))
        plt.plot(x, y)
        y = np.array(range(self.triple_point[1], 30))
        x = eval("{}*(y-{})**3+{}".format(1 / self.gradient, self.triple_point[1], self.triple_point[0]))
        plt.plot(x, y)
        plt.scatter(TEST_X, TEST_Y)
        plt.show()

    def chemical_type(self, temperature, pressure):
        if temperature < self.triple_point[0] and (self.triple_point[1] * temperature ** 2) / (
                self.triple_point[0] ** 2) >= pressure:
            return "gas"
        elif temperature < self.triple_point[0] and pressure < self.triple_point[1]:
            return "solid"
        elif temperature < 1 / self.gradient * (pressure - self.triple_point[1]) ** 3 + self.triple_point[
            0] and pressure > self.triple_point[1]:
            return "solid"
        elif temperature < self.critical_point[0] and pressure <= ((self.triple_point[1] - self.critical_point[1]) / (
                self.critical_point[0] - self.triple_point[0] ** 2)) * (temperature - self.triple_point[0]) ** 2 + \
                self.triple_point[1]:
            return "gas"
        elif temperature >= self.critical_point[0] and pressure <= self.critical_point[0]:
            return "gas"
        elif temperature >= self.critical_point[0]:
            return "super_critical"
        else:
            return "liquid"


class InnerPlanet:
    def __init__(self, radius, composition, temperature):
        self.layers = []
        self.radius = radius
        self.composition = composition
        self.volume = sphere_volume(radius)
        for i in range(radius // 1):
            self.layers.append(Layer(composition, (sphere_volume(i + 1) - sphere_volume(i)) / self.volume))

    def display(self, table_of_chemicals):
        fig, ax = plt.subplots(1, 1, figsize=(2 * 6, 2 * 3.5), subplot_kw={'projection': 'polar'})
        for i, layer in enumerate(self.layers):
            composition = layer.composition
            red_total = 0
            blue_total = 0
            green_total = 0
            weights = 0
            for key, value in composition.items():
                red, green, blue = table_of_chemicals[key].colour
                red_total += red*value
                blue_total += blue*value
                green_total += green*value
                weights += value
            display_colour = ((red_total/weights), (green_total/weights), (blue_total/weights))
            print(display_colour)
            x = np.arange(0.0, 2 * np.pi, 0.01)
            y = i
            ax.fill_between(x, y, color=display_colour)
        plt.show()

    def evolve_structure(self):
        """evolves the structure"""


class Layer:
    def __init__(self, composition, percentage):
        self.composition = {}
        for key, value in list(composition.items()):
            self.composition[key] = value * percentage


main()
