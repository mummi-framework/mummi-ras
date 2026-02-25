"""
Simple interpolator class for testing purposes.
"""

import random
import numpy as np
from mummi_ras.ml.ls_point import LSPoint

class NaiveInterpolator:

    def __init__(self, source, target):
        """
        Constructor for the naive interpolator.
        """

        self.source = source
        self.target = target

        return 


    def draw_k_points(self, k):
        """
        Draw k points via a naive point to point interpolation.
        """

        new_points = [self.select_point() for i in range(k)]

        return new_points
    
    def select_point(self):
        """
        Select a single new point.
        """
        # choose one point in A and B each
        pnt_a = random.choice(self.source)
        pnt_b = random.choice(self.target)

        _vec = pnt_b.get_coordinates() - pnt_a.get_coordinates()

        lmbd = np.random.uniform()

        new_coord = lmbd * pnt_b.get_coordinates() + (1 - lmbd) * pnt_a.get_coordinates() 

        return new_coord



