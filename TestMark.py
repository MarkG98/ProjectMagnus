from modsim import *
from matplotlib import pyplot as plt

# declare all units
m = UNITS.meter
s = UNITS.second
kg = UNITS.kilogram
degree = UNITS.degree
radian = UNITS.radian

hero_height = 1.5

condition = Condition(x = 2 * m,
                      y = (365.76 + hero_height) * m,
                      g = 9.8 * m/s**2,
                      mass = 145e-3 * kg,
                      diameter = 73e-3 * m,
                      rho = 1.2 * kg/m**3,
                      C_d = 0.3,
                      angle = -60 * degree,
                      velocity = 70/2 * m/s,
                      w = 40/2 * radian/s,
                      duration = 15 * s)

def make_system(condition):
    theta = np.deg2rad(condition.angle)
    vx, vy = pol2cart(theta, condition.velocity)
    init = State(x = condition.x, y = condition.y, vx = vx, vy = vy)
    area = np.pi * (condition.diameter / 2) ** 2
    ts = linspace(0, condition.duration, 101)

    return System(init = init, g = condition.g, mass = condition.mass, area = area, rho = condition.rho, C_d = condition.C_d, ts = ts)

baseballSystem = make_system(condition)

def slope_func(state, t, system):
    x, y, vx, vy = state
    unpack(system)

    a_grav = Vector(0, -g, 0) # acceleration due to gravity
    v_ball = Vector(vx, vy) # ball velocity vector

    # creates unit vector in direction of the magus force
    x, y = pol2cart(v_ball.angle - (pi / 2) * radian, 1)
    magnus_direction = Vector(x, y, 0)

    w_vector = Vector(0, 0, condition.w)

    # calculates acceleration due to the magnus force
    #f_magnus = (pi ** 2) * ((condition.diameter / 2) ** 3) * rho * v_ball.mag * condition.w * magnus_direction.hat()
    f_magnus = ((pi ** 2) * ((condition.diameter / 2) ** 3) * rho) * w_vector.cross(v_ball)
    a_magnus = f_magnus / mass

    # calculates acceleration due to drag
    f_drag = -rho * v_ball.mag * v_ball * C_d * (area / 2)
    a_drag = f_drag / mass

    a_total = a_grav + a_magnus + a_drag

    return v_ball.x, v_ball.y, a_total.x, a_total.y

#slope_func(baseballSystem.init, 0, baseballSystem)

run_odeint(baseballSystem, slope_func)

xs = baseballSystem.results.x
ys = baseballSystem.results.y

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Ball flight Path')
plt.grid(True)
plt.axvline(x = 2, color='r')
plt.plot(xs, ys)
plt.show()
