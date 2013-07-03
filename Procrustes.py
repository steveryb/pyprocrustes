#!/usr/bin/env python
"""
Procrustes.py: An implementation of Partial Ordinary Procrustes Analysis in
Python 2 for three dimensional curves using numpy.

Procrustes was a bandit who, in Greek mythology, had a rather esoteric modus
operandi: he would make his victims fit inside his bed by either stretching
them if they were too short or cutting off pieces of them if they were too
long (who said statistical analysis couldn't be morbid!).

Partial Ordinary Procrustes Analysis (which we'll be referring to as POPA) 
takes in two curves: a reference curve and a curve
to compare to the reference curve. It then compares the shape of these two
curves by removing the translational and scaling componenents of the two shapes
and then removes the rotational component from the curve being compared to get
the curves optimally "superimposed". Then, we apply some distance function to
score their closeness.
"""
import math
import numpy
class Procrustes(object):
    def __init__(self,reference_curve):
        """
        __init__(numpy.ndarray) 
        Procrustes accepts a reference curve representing a three
        dimensional curve which it will use to compare to all other curves
        given to it to analyse

        reference_curve must be a numpy nd of three tuples each representing a
        point in three dimensional space.
        """
        # Store the reference curve scaled and translated
        if not(isinstance(reference_curve,numpy.ndarray)):
            raise TypeError("Reference curve must be a numpy.ndarray")
        self.reference_curve = self.scale(self.translate(reference_curve))

    @staticmethod
    def translate(curve):
        """
        P.translate(numpy.ndarray) -> numpy.ndarray
        
        Take a three dimensional curve defined by a numpy ndarray and translate
        that curve so that the mean of each of spatial componenets lies at the
        origin.
        """
        try:
            return curve - sum(curve)/curve.shape[0]
        except AttributeError:
            raise TypeError("Scale accepts only numpy arrays")
    @staticmethod
    def scale(curve):
        """
        P.scale(numpy.ndarray) -> numpy.ndarray
        
        Take a three dimensional curve defined by a numpy ndarray and translate
        that curve so that the root mean square deviation from the points to the
        origin is 1.
        """
        return curve/Procrustes.rmsd(curve)

    @staticmethod
    def rmsd(curve):
        """
        P.rmsd(numpy.ndarray) -> float

        Give the root mean square deviation of a three dimensional curve's
        points from the origin
        """
        sum_of_squares = sum(sum(pow(curve,2)))
        return pow(sum_of_squares/curve.shape[0],0.5)

    @staticmethod
    def rotate(ref_curve, curve):
        """
        P.rotate(numpy.ndarray,numpy.ndarray) -> numpy.ndarray

        Rotate the curve around the reference curve in order to minimise the
        root mean squared deviation of the one curve from the other. This uses
        Kabsch's algorithm to calculate the optimum rotation matrix to this.
        """
        # step 1) calculate co_variance matrix
        co_variance_matrix = numpy.transpose(ref_curve).dot(curve)
        # step 2) use single value decomposition to represent our co_variance
        # matrix as the product of matrices v s w.
        v,s,w = numpy.linalg.svd(co_variance_matrix)
        # decide whether we need to adjust sign of rotation matrix to insure a
        # right-handed co-ordinate system
        sign = lambda x: 1 if x > 0 else 0 if x == 0 else -1
        v_t,w_t = [numpy.transpose(x) for x in (v,w)]
        d = sign(numpy.linalg.det(w_t.dot(v_t)))
        # finally, calculate the optimum rotation matrix
        rotation_matrix = \
                        w_t.dot(numpy.array([[1,0,0],[0,1,0],[0,0,d]])).dot(v_t)
        # now apply that rotation matrix to every point in the curve
        return numpy.array([rotation_matrix.dot(point) for point in
        curve])
    @staticmethod
    def superposition(ref_curve,curve):
        """
        P.superposition(numpy.ndarray,numpy.ndarray) -> 
                                                     (numpy.ndarray,numpy.ndarray)

        Convenience function to perform a procrustes superposition on two curves,
        using the reference curve as the template for the rotation of the other
        curve. The output will be a tuple with the refence curve and the other
        curve superposed respectively.
        """
        s_ref_curve,s_curve = [Procrustes.scale(Procrustes.translate(curve)) 
                                                for curve in [ref_curve,curve]]
        return (s_ref_curve, Procrustes.rotate(s_ref_curve,s_curve))
    
    @staticmethod
    def distance(ref_curve,curve):
        pass
if __name__ == "__main__":
    a = Procrustes(numpy.array([[1,2,3],[4,5,6],[7,8,9]]))

    print(Procrustes.rotate(a.reference_curve,numpy.array([[1,2,3],[4,5,6],[7,8,9]])))
    
    
