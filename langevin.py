"""
Module for Langevin integration methods
"""
import numpy as np

class Langevin:
    def __init__(self, TL, alpha, h):
        self.TL = TL
        self.alpha = alpha
        self.h = h
