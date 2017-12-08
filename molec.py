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

lib = ctypes.cdll.LoadLibrary("./molec.so")

class LennardJones:
    def __init__(self, centers, m, rmin, epsilon, T, alpha, h,
                 rc = "default", box_width = 30.):
        nparticles = len(centers)
        centers_array = b""
        for p in centers:
            centers_array += struct.pack("ff", *p)

        if rc == "default":
            constructor = lib.LennardJones_new
            constructor.argtypes = (ctypes.c_int,     # nparticles
                                    ctypes.c_char_p,  # centers
                                    ctypes.c_float,   # m
                                    ctypes.c_float,   # rmin
                                    ctypes.c_float,   # epsilon
                                    ctypes.c_float,   # T
                                    ctypes.c_float,   # alpha
                                    ctypes.c_float,   # h
                                    ctypes.c_float)   # box_width
            self.obj = constructor(nparticles, centers_array, m, rmin,
                                   epsilon, T, alpha, h, box_width)
        else:
            constructor = lib.LennardJones_new_full
            constructor.argtypes = (ctypes.c_int,     # nparticles
                                    ctypes.c_char_p,  # centers
                                    ctypes.c_float,   # m
                                    ctypes.c_float,   # rmin
                                    ctypes.c_float,   # epsilon
                                    ctypes.c_float,   # T
                                    ctypes.c_float,   # alpha
                                    ctypes.c_float,   # h
                                    ctypes.c_float,   # rc
                                    ctypes.c_float)   # box_width
            self.obj = constructor(nparticles, centers_array, m, rmin,
                                   epsilon, T, alpha, h, rc, box_width)

    def get_box_width(self):
        ans_call = lib.LennardJones_get_box_width
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_T(self):
        ans_call = lib.LennardJones_get_T
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_alpha(self):
        ans_call = lib.LennardJones_get_alpha
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_h(self):
        ans_call = lib.LennardJones_get_h
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_m(self):
        ans_call = lib.LennardJones_get_m
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_rmin(self):
        ans_call = lib.LennardJones_get_rmin
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_rc(self):
        ans_call = lib.LennardJones_get_rc
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_epsilon(self):
        ans_call = lib.LennardJones_get_epsilon
        ans_call.restype = ctypes.c_float
        return ans_call(self.obj)

    def get_nparticles(self):
        ans_call = lib.LennardJones_get_nparticles
        ans_call.restype = ctypes.c_int
        return ans_call(self.obj)

    def step(self):
        return lib.LennardJones_step(self.obj)

    def param_change(self, T, alpha, h, box_width, m = "nochange",
                     rmin = "nochange", rc = "nochange", epsilon = "nochange"):
        if all(np.array([ m, rmin, rc, epsilon]) == "nochange"):
            param_call = lib.LennardJones_param_change
            param_call.argtypes = (ctypes.c_void_p, # sys
                                   ctypes.c_float,  # T
                                   ctypes.c_float,  # alpha
                                   ctypes.c_float,  # h
                                   ctypes.c_float)  # box_width
            param_call(self.obj, T, alpha, h)
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
                                   ctypes.c_float,  # T
                                   ctypes.c_float,  # alpha
                                   ctypes.c_float,  # h
                                   ctypes.c_float,  # box_width
                                   ctypes.c_float,  # m
                                   ctypes.c_float,  # rmin
                                   ctypes.c_float,  # rc
                                   ctypes.c_float)  # epsilon
            param_call(self.obj, T, alpha, h, box_width, m, rmin, rc, epsilon)

    def deconstruct(self):
        """
        Frees memory associated with the simulation
        """
        lib.LennardJones_delete(self.obj)
        self.obj = None

def scatter_particles(centers, show=True):
    centers_x = [c[0] for c in centers]
    centers_y = [c[1] for c in centers]
    plt.scatter(centers_x, centers_y)
    if show:
        plt.show()

def dist(p1, p2):
    assert(len(p1) == len(p2) == 2)
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
        
def not_too_close_random(box_width, centers_so_far, threshold = 1,
                         max_tries = 100, verbose = True):
    ans = (0, 0)
    try_again = True
    iteration = 0
    while try_again and iteration < max_tries:
        ans = (random.uniform(0, 1) * box_width, random.uniform(0, 1) * 30)
        try_again = False
        for c in centers_so_far:
            if dist(ans, c) < threshold:
                try_again = True
                break
        iteration += 1
        if verbose:
            print("iteration: %s" % str(iteration))
    return ans

if __name__ == "__main__":
    box_width = 30.
    n = 50
    centers = []
    for i in range(n):
        centers.append(not_too_close_random(box_width, centers))
    print(centers)
    scatter_particles(centers)
    simulation = LennardJones(centers, m = 1, rmin = 1, epsilon = 1,
                              T = 100, alpha = 0.2, h = 0.001,
                              rc = "default", box_width = 30.)
    print(simulation.get_nparticles())
    simulation.deconstruct()
