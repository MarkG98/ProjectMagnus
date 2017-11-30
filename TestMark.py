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
                      angle = -10 * degree,
                      velocity = 70 * m/s,
                      w = 61.4674989621 * radian/s,
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
    v_ball = Vector(vx, vy,0) # ball velocity vector

    # creates unit vector in direction of the magus force
    x, y = pol2cart(v_ball.angle - (pi / 2) * radian, 1)
    magnus_direction = Vector(x, y, 0)

    w_vector = Vector(0, 0, - condition.w)

    # calculates acceleration due to the magnus force
    #f_magnus = (pi ** 2) * ((condition.diameter / 2) ** 3) * rho * v_ball.mag * condition.w * magnus_direction.hat()
    f_magnus = ((pi ** 2) * ((condition.diameter / 2) ** 3) * rho) * w_vector.cross(v_ball)
    a_magnus = f_magnus / mass

    # calculates acceleration due to drag
    f_drag = -rho * v_ball.mag * v_ball * C_d * (area / 2)
    a_drag = f_drag / mass

    a_total = a_grav + a_magnus + a_drag

    return v_ball.x, v_ball.y, a_total.x, a_total.y

run_odeint(baseballSystem, slope_func)
xs = baseballSystem.results.x
ys = baseballSystem.results.y

def error_func(w):
  condition.set(w = w * radian/s)
  baseballSystem = make_system(condition)
  run_odeint(baseballSystem, slope_func)
  T = interpolate(baseballSystem.results.y)
  return T(0)

def sweep_func():
  for angle in linspace(-10,-80,2):
    condition.set(angle=angle * degree)
    for velocity in linspace(30,150,2):
      condition.set(velocity=velocity* m/s)
      for w in linspace(20,120,2):
        condition.set(w=w* radian/s)
        baseballSystem = make_system(condition)
        run_odeint(baseballSystem,slope_func)
        xs = baseballSystem.results.x
        ys = baseballSystem.results.y
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.title('Ball flight Path')
        plt.xlim([0,70])
        plt.ylim([0,400])
        plt.grid(True)
        plt.axvline(x = 2, color='r')
        plt.plot(xs, ys)
        plt.show()

def interpolate_range(results,y_value):
  ys = results
  descent = ys
  T = interp_inverse(descent)
  t_land = T(y_value)
  return t_land

def error_func(w):
  condition.set(w = w * radian/s)
  baseballSystem = make_system(condition)
  run_odeint(baseballSystem, slope_func)
  T = interpolate(baseballSystem.results.x)
  t_land = interpolate_range(baseballSystem.results.y, 0)
  return T(t_land)


#value = error_func(w = 67, interpolate_range = interpolate_range)
solution = fsolve(error_func, 67)
w_ideal = solution[0]
print(w_ideal)

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Ball flight Path')
plt.xlim([0,70])
plt.ylim([0,400])
plt.grid(True)
plt.axvline(x = 2, color='r')
plt.plot(xs, ys)
plt.show()
