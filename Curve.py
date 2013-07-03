#!/usr/bin/env python
"""
Curve.py: a representation of a three dimensional curve intended to simplify
reparameterisation and distance calculations
"""
import math
import numpy

class Curve(object):
    def __init__(self, points_ndarray):
        """
        C.__init__(numpy.ndarray)

        Takes an numpy array with each entry defining a point in three
        dimensional space and creates a new curve object to represent it using
        linear interpolation between the points.
        """
        self.set_points(points_ndarray)

        

    def val_at_arclength(self,t):
        """
        C.val_at_arclength(float) -> [float,float,float]

        Give the point corresponding to the given arclength for the curve.

        To get this, linear interpolation on the given points is used.
        """
        if t==0:
            return self.points[0]
        if t==self.arc_length:
            return self.points[-1]
        if not(0 < t < self.arc_length):
            return \
                ValueError("Input needs to be between 0 and the arclength of the curve")

        # find which two points this arclength t lies between
        top_limit = len(self.distances)-1
        bottom_limit = 0
        while top_limit != bottom_limit:
            middle = int((top_limit-bottom_limit)//2)
            top_dist, bottom_dist = self.distances[middle]

            if bottom_dist <= t < top_dist:
                top_limit = bottom_limit = middle
            elif t < bottom_dist:
                top_limit = middle
            elif t >= top_dist:
                bottom_limit = middle+1
        
        a,b = self.points[top_limit],self.points[top_limit+1]
        # contruct a function to linearly interpolate between these points
        dist_a = self.distances[top_limit][0]
        dist_b = self.distances[top_limit][1]
        f = lambda t: (b-a)/(dist_b-dist_a)*(t-dist_a) + a
        # return the interpolated value
        return f(t)

    def gen_num_points(self, total, loose_detail=False):
        """
        C.change_num_points(int)

        Using the arclength parameterisation of this curve, create a new curve
        with the new number of points.

        Note: to avoid losing information, this will default to only allowing larger
        number of points to be chosen. If this what you desire, set
        loose_detail to True.
        """
        if total < self.points.shape[0] and not loose_detail:
            raise ValueError(("Total ({0}) needs to be larger than number of points"+
            " present ({1}). If you want to loose detail, invoke this function with"+
            "loose_detail=True").format(total,self.points.shape[0]))
        return numpy.array([self.val_at_arclength(t) for t in 
                                            numpy.linspace(0,self.arc_length,total)])
    
    def set_points(self,points_ndarray):
        """
        C.set_points(numpy.ndarray)

        Change the points this curve represents
        """
        if not(isinstance(points_ndarray,numpy.ndarray)):
            raise TypeError("Points must be in a numpy.ndarray")
        if not(points_ndarray.shape[1] == 3 or points_ndarray.shape[0] >0):
            raise TypeError("points_ndarray must have more than one point "
                            "and each point must be an array of length three.")

        self.points = points_ndarray
        euc_length = lambda a,b: pow(sum(pow(a-b,2)),0.5)
        self.distances = [(0,euc_length(self.points[0],self.points[1]))]
        for i in range(1,len(self.points)-1):
            self.distances.append((self.distances[i-1][1],
                    self.distances[i-1][1]+euc_length(self.points[i],self.points[i+1])))
        self.arc_length = float(self.distances[-1][1])


    def __str__(self):
        """
        s.__str___()

        Give the string representation of the points this curve currently is
        storing.
        """
        return str(self.points)

if __name__=="__main__":
    c = Curve(numpy.array([[1,1,1],[2,2,2],[3,3,3]]))
    print(c)
    print str(c.gen_num_points(20))
