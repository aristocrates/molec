"""
Module for Langevin integration methods
"""
import numpy as np
# if __name__ == "__main__":
#     import matplotlib
#     matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ctypes
import struct
import random
import time

lib = ctypes.cdll.LoadLibrary("./molec.so")

class LennardJones:
    def __init__(self, centers, m, rmin, epsilon, T, alpha, h,
                 rc = "default", box_width = 30.):
        nparticles = len(centers)
        centers_array = b""
        for p in centers:
            centers_array += struct.pack("dd", *p)

        if rc == "default":
            constructor = lib.LennardJones_new
            constructor.argtypes = (ctypes.c_int,     # nparticles
                                    ctypes.c_char_p,  # centers
                                    ctypes.c_double,   # m
                                    ctypes.c_double,   # rmin
                                    ctypes.c_double,   # epsilon
                                    ctypes.c_double,   # T
                                    ctypes.c_double,   # alpha
                                    ctypes.c_double,   # h
                                    ctypes.c_double)   # box_width
            self.obj = constructor(nparticles, centers_array, m, rmin,
                                   epsilon, T, alpha, h, box_width)
        else:
            constructor = lib.LennardJones_new_full
            constructor.argtypes = (ctypes.c_int,     # nparticles
                                    ctypes.c_char_p,  # centers
                                    ctypes.c_double,   # m
                                    ctypes.c_double,   # rmin
                                    ctypes.c_double,   # epsilon
                                    ctypes.c_double,   # T
                                    ctypes.c_double,   # alpha
                                    ctypes.c_double,   # h
                                    ctypes.c_double,   # rc
                                    ctypes.c_double)   # box_width
            self.obj = constructor(nparticles, centers_array, m, rmin,
                                   epsilon, T, alpha, h, rc, box_width)

    def get_box_width(self):
        ans_call = lib.LennardJones_get_box_width
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_T(self):
        ans_call = lib.LennardJones_get_T
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_alpha(self):
        ans_call = lib.LennardJones_get_alpha
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_h(self):
        ans_call = lib.LennardJones_get_h
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_m(self):
        ans_call = lib.LennardJones_get_m
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_rmin(self):
        ans_call = lib.LennardJones_get_rmin
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_rc(self):
        ans_call = lib.LennardJones_get_rc
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_epsilon(self):
        ans_call = lib.LennardJones_get_epsilon
        ans_call.restype = ctypes.c_double
        return ans_call(self.obj)

    def get_nparticles(self):
        ans_call = lib.LennardJones_get_nparticles
        ans_call.restype = ctypes.c_int
        return ans_call(self.obj)

    def get_centers(self):
        ans_call = lib.LennardJones_get_centers
        ans_call.restype = ctypes.POINTER(ctypes.c_double)
        return ans_call(self.obj)

    def get_center_data(self, plottable = False):
        """
        Returns the coordinates of each center point in a usable
        numpy array format
        """
        # gotten_centers = np.ctypeslib.as_array(gotten_centers, (
        return list_to_centers_tuple(np.ctypeslib.as_array(self.get_centers(), (2 * self.get_nparticles(),)), plot = plottable)

    def step(self):
        return lib.LennardJones_step(self.obj)

    def param_change(self, T, alpha, h, box_width, m = "nochange",
                     rmin = "nochange", rc = "nochange", epsilon = "nochange"):
        if all(np.array([ m, rmin, rc, epsilon]) == "nochange"):
            param_call = lib.LennardJones_param_change
            param_call.argtypes = (ctypes.c_void_p, # sys
                                   ctypes.c_double,  # T
                                   ctypes.c_double,  # alpha
                                   ctypes.c_double,  # h
                                   ctypes.c_double)  # box_width
            param_call(self.obj, T, alpha, h, box_width)
        else:
            if m == "nochange":
                m = self.get_m()
            if rmin == "nochange":
                rmin = self.get_rmin()
            if rc == "nochange":
                rc = self.get_rc()
            if epsilon == "nochange":
                epsilon = self.get_epsilon()
            param_call = lib.LennardJones_param_change_full
            param_call.argtypes = (ctypes.c_void_p, # sys
                                   ctypes.c_double,  # T
                                   ctypes.c_double,  # alpha
                                   ctypes.c_double,  # h
                                   ctypes.c_double,  # box_width
                                   ctypes.c_double,  # m
                                   ctypes.c_double,  # rmin
                                   ctypes.c_double,  # rc
                                   ctypes.c_double)  # epsilon
            param_call(self.obj, T, alpha, h, box_width, m, rmin, rc, epsilon)

    def deconstruct(self):
        """
        Frees memory associated with the simulation
        """
        lib.LennardJones_delete(self.obj)
        self.obj = None

def scatter_particles(centers, xlim = None, ylim = None, show=True):
    centers_x = [c[0] for c in centers]
    centers_y = [c[1] for c in centers]
    plt.scatter(centers_x, centers_y)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if show:
        plt.show()

def dist(p1, p2):
    assert(len(p1) == len(p2) == 2)
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
        
def not_too_close_random(init_min, width, centers_so_far, threshold = 1,
                         max_tries = 100, verbose = True):
    ans = (0, 0)
    try_again = True
    iteration = 0
    while try_again and iteration < max_tries:
        ans = (init_min + random.uniform(0, 1) * width, init_min + random.uniform(0, 1) * width)
        try_again = False
        for c in centers_so_far:
            if dist(ans, c) < threshold:
                try_again = True
                break
        iteration += 1
        if verbose:
            print("iteration: %s" % str(iteration))
    return ans

def get_initial_data(nelements, box_width, init_min, init_width, iterations,
                     steps_per_update = 100, step_size = 1e-4,
                     start_temp = 0.01, cool_iter = 50,
                     decay_factor = 0.25, verbose = False):
    """
    Cools a system down into a lattice

    By default uses a small step size to obtain sensible initial conditions
    without blowing up
    """
    centers = []
    for i in range(nelements):
        centers.append(not_too_close_random(init_min, init_width, centers,
                                            verbose = False))
    simulation = LennardJones(centers, m = 1, rmin = 1., epsilon = 0.8,
                              T = start_temp, alpha = 0.2, h = step_size,
                              rc = "default", box_width = box_width)
    try:
        k = 0
        for i in range(iterations):
            # every 50 cool the temperature down
            k += 1
            if k >= cool_iter:
                T = simulation.get_T()
                alpha = simulation.get_alpha()
                h = simulation.get_h()
                box_width = simulation.get_box_width()
                simulation.param_change(T * decay_factor, alpha, h, box_width)
                if verbose:
                    print("New temperature: %s" % str(simulation.get_T()))
                k = 0
            for j in range(steps_per_update):
                simulation.step()
            current_centers = list_to_centers_tuple(np.ctypeslib.as_array(simulation.get_centers(), (2 * nelements,)), plot = False)
    finally:
        simulation.deconstruct()

    return_centers = []
    for cent in current_centers:
        if 1 < cent[0] < box_width and 1 < cent[1] < box_width:
            print(cent)
            print(cent[0])
            print(cent[1])
            # for some unknown reason, appending the numpy array itself
            # causes the values to become garbage, possibly because something
            # is overwriting something else
            return_centers.append([cent[0], cent[1]])
    return return_centers

def list_to_centers_tuple(centers_list, plot = False):
    """
    Assumes centers_list is a numpy array
    """
    assert(len(centers_list) % 2 == 0)
    if not plot:
        return np.reshape(centers_list, (int(len(centers_list) / 2), 2))
    else:
        return np.transpose(np.reshape(centers_list, (int(len(centers_list) / 2), 2)))

if __name__ == "__main__":
    box_width = 30.
    init_min = 8.
    init_width = 14.
    centers = [(10, 10)]
    n = 100
    #n = 2
    #centers = [(10, 10), (10, 11)]
    centers = []
    for i in range(n):
        centers.append(not_too_close_random(init_min, init_width, centers))
    print(centers)
    scatter_particles(centers)
    simulation = LennardJones(centers, m = 1, rmin = 1., epsilon = 0.8,
                              T = 0.01, alpha = 0.2, h = 0.001,
                              rc = "default", box_width = box_width)
    try:
        # print out every parameter as a sanity check
        print("nparticles: " + str(simulation.get_nparticles()))
        print("mass: " + str(simulation.get_m()))
        print("rmin: " + str(simulation.get_rmin()))
        print("epsilon: " + str(simulation.get_epsilon()))
        print("T: " + str(simulation.get_T()))
        print("alpha: " + str(simulation.get_alpha()))
        print("h: " + str(simulation.get_h()))
        print("rc: " + str(simulation.get_rc()))
        print("box_width: " + str(simulation.get_box_width()))
        gotten_centers = simulation.get_centers()
        gotten_centers = np.ctypeslib.as_array(gotten_centers, (2 * n,))
        print(gotten_centers)
        simulation.step()
        gotten_centers2 = np.ctypeslib.as_array(simulation.get_centers(),
                                                (2 * n,))
        print(list_to_centers_tuple(gotten_centers2))
        plt.axis([0, 30, 0, 30])
        plt.ion()
        k = 0
        for i in range(1000):
            # every 50 cool the temperature down
            k += 1
            if k >= 50:
                T = simulation.get_T()
                alpha = simulation.get_alpha()
                h = simulation.get_h()
                box_width = simulation.get_box_width()
                simulation.param_change(T * 0.25, alpha, h, box_width)
                print("New temperature: %s" % str(simulation.get_T()))
                k = 0
            for j in range(10):
                simulation.step()
            current_centers = list_to_centers_tuple(np.ctypeslib.as_array(simulation.get_centers(), (2 * n,)), plot = True)
            #print(current_centers)
            plt.clf()
            plt.axis([0, 30, 0, 30])
            plt.scatter(*current_centers)
            plt.pause(0.001)
        while True:
            plt.pause(0.1)
    finally:
        simulation.deconstruct()
