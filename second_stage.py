DENSITY_CONSTANT = 1


class Layer:
    def __init__(self, heat, material, heat_transfer):
        self.temperature = heat
        self.material = material
        self.heat_transfer = heat_transfer

    def density(self, pressure):
        return DENSITY_CONSTANT * pressure/self.temperature

    def transfer_heat(self, above, below):
        amount_of_heat = (self.temperature - below.heat) * self.heat_transfer * below.heat_transfer


if __name__ == '__main__':
    for i in range(1000):
