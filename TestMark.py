from modsim import *
from matplotlib import pyplot as plt

# declare all units
m = UNITS.meter
s = UNITS.second
kg = UNITS.kilogram
degree = UNITS.degree
radian = UNITS.radian

hero_height = 2

condition = Condition(x = 2 * m,
                      y = (365.76 + hero_height) * m,
                      g = 9.8 * m/s**2,
                      mass = 145e-3 * kg,
                      diameter = 73e-3 * m,
                      rho = 1.2 * kg/m**3,
                      C_d = 0.3,
                      angle = -60 * degree,
                      velocity = 50 * m/s,
                      w = 80 * radian/s,
                      duration = 15 * s)

def make_system(condition):
    theta = np.deg2rad(condition.angle)
    vx, vy = pol2cart(theta, condition.velocity)
    init = State(x = condition.x, y = condition.y, vx = vx, vy = vy)
    area = np.pi * (condition.diameter / 2) ** 2
    ts = linspace(0, condition.duration, 101)

    return System(init = init, g = condition.g, mass = condition.mass, area = area, rho = condition.rho, C_d = condition.C_d, ts = ts)

make_system(condition)
