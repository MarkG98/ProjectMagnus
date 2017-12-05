from modsim import *
from matplotlib import pyplot as plt

# declare all units
m = UNITS.meter
s = UNITS.second
kg = UNITS.kilogram
degree = UNITS.degree
radian = UNITS.radian

hero_height = 1.5

"""condition = Condition(x = 2 * m,
                      y = (365.76 + hero_height) * m,
                      g = 9.8 * m/s**2,
                      mass = 145e-3 * kg,
                      diameter = 73e-3 * m,
                      rho = 1.2 * kg/m**3,
                      C_d = 0.3,
                      angle = 7 * degree,
                      velocity = 95 * m/s,
                      w = 70.97764981402517 * radian/s,
                      duration = 15 * s)"""
condition = Condition(x = 2,
                      y = (365.76 + hero_height),
                      g = 9.8,
                      mass = 145e-3,
                      diameter = 73e-3,
                      rho = 1.2,
                      C_d = 0.3,
                      angle = 7,
                      velocity = 95,
                      w = 70.97764981402517,
                      duration = 15)

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

def sweep_func():
  for angle in linspace(-10,-80,2):
    """condition.set(angle=angle * degree)"""
    condition.set(angle=angle)
    for velocity in linspace(30,150,2):
      """condition.set(velocity=velocity* m/s)"""
      condition.set(velocity=velocity)
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
  """condition.set(w = w * radian/s)"""
  condition.set(w = w)
  baseballSystem = make_system(condition)
  run_odeint(baseballSystem, slope_func)
  X = interpolate(baseballSystem.results.x)
  t_land = interpolate_range(baseballSystem.results.y, 0.78)
  return X(t_land)

def heightAtDoor(system):
    xs = system.results.x
    ys = system.results.y

    t_maxRange = xs.idxmax()

    finalDescent = xs.loc[t_maxRange:]
    T = interp_inverse(finalDescent)
    t_throughDoor = T(2)

    Y = interpolate(ys, kind='cubic')
    return Y(t_throughDoor)


#value = error_func(w = 67, interpolate_range = interpolate_range)
#solution = fsolve(error_func, 67)
#w_ideal = solution[0]
#print(w_ideal)

angle_array = linspace(6.5, 7.5, 4)

def angleSweep(system, angle_array):
    for angle in angle_array:
        condition.set(angle = angle * degree)
        solution = fsolve(error_func, 60)
        spinForButt = solution[0]

        condition.set(w = spinForButt * radian/s)
        system = make_system(condition)
        run_odeint(system, slope_func)

        height = heightAtDoor(system)

        print("angle = ", angle, "w = ", condition.w, "height at door = ", height, '\n')

#angleSweep(baseballSystem, angle_array)
velocity_array = linspace(94,96,11)
def velocitySweep(system,velocity_array, angle_array):
    for angle in angle_array:
        heights = []
        condition.set(angle=angle)
        for velocity in velocity_array:
            """condition.set(velocity = velocity * m/s)"""
            condition.set(velocity = velocity)
            solution = fsolve(error_func, 60)
            spinForButt = solution[0]

            """condition.set(w = spinForButt * radian/s)"""
            condition.set(w = spinForButt)
            system = make_system(condition)
            run_odeint(system, slope_func)

            height = heightAtDoor(system)
            heights.append(height)
            print("angle =", angle, "velocity = ", velocity, "w = ", condition.w, "height at door = ", height, '\n')

        plt.plot(velocity_array, heights)
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Height From Ground (m)')
        plt.title('Linear Velocity Optimization')
        #plt.xlim([0,70])
        #plt.ylim([0,400])
        plt.grid(True)
        plt.axhline(y = 2, color='r')


    plt.show()

#velocitySweep(baseballSystem,velocity_array,angle_array)

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Ball flight Path')
plt.xlim([0,250])
plt.ylim([0,400])
plt.grid(True)
plt.axvline(x = 2, color='r')
plt.plot(xs, ys)
plt.show()
