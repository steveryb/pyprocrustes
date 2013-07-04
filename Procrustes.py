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
from Curve import Curve

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


def rmsd(curve):
    """
    P.rmsd(numpy.ndarray) -> float

    Give the root mean square deviation of a three dimensional curve's
    points from the origin
    """
    sum_of_squares = sum(sum(pow(curve,2)))
    return pow(sum_of_squares/curve.shape[0],0.5)


def scale(curve):
    """
    P.scale(numpy.ndarray) -> numpy.ndarray
    
    Take a three dimensional curve defined by a numpy ndarray and translate
    that curve so that the root mean square deviation from the points to the
    origin is 1.
    """
    return curve/rmsd(curve)



def rotate(ref_curve, curve):
    """
    P.rotate(numpy.ndarray,numpy.ndarray) -> (numpy.ndarray,numpy.ndarray)

    Rotate the curve around the reference curve in order to minimise the
    root mean squared deviation of the one curve from the other. This uses
    Kabsch's algorithm to calculate the optimum rotation matrix to this.
    This expects both curves to have the same number of points.

    Return the reference curve and the rotated curve respectively.
    """
    if(ref_curve.shape[0] != curve.shape[0]):
        return TypeError("Curves have different numbers of points")
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


def superposition(ref_curve,curve):
    """
    P.superposition(numpy.ndarray,numpy.ndarray) -> 
                                                 (numpy.ndarray,numpy.ndarray)

    Convenience function to perform a procrustes superposition on two curves,
    using the reference curve as the template for the rotation of the other
    curve. The output will be a tuple with the refence curve and the other
    curve superposed respectively.

    Note that if the curves don't have the same number of points, this will
    take the curve with the smaller number of points and reparameterise the
    linear interpolation of its points using arclength.
    """
    # firstly, we need to make sure the ref_curve and the other curve have
    # the same number of points.
    num_points = max([c.shape[1] for c in [ref_curve,curve]])
    if ref_curve.shape[0] > curve.shape[0]:
        curve = Curve(curve).gen_num_points(ref_curve.shape[0])
    else:
        ref_curve = Curve(ref_curve).gen_num_points(curve.shape[0])
    s_ref_curve,s_curve = [scale(translate(curve)) 
                                            for curve in
                                            [ref_curve,curve]]
    return (s_ref_curve, rotate(s_ref_curve,s_curve))


def min_distance(ref_curve,curve):
    """
    P.min_distance(numpy.ndarray,numpy.ndarray) -> (int)

    For each point in a curve, find the smallest distance from it to the
    reference curve.

    Returns an array of tuples, each corresponding to the original
    points on the curve giving the minimum distance to the reference curve.
    """
    euc_length = lambda a,b: pow(sum(pow(a-b,2)),0.5)
    ref_curve_c = Curve(ref_curve)
    print(str(ref_curve_c),curve)
    return\
        tuple((euc_length(ref_curve_c.find_nearest_point(point),point) for point in curve))

def procrustes_distance(ref_curve,curve):
    """
    P.procrustes_distance(numpy.ndarray,numpy.ndarray) -> (int)

    Given a curve, use procrustes analysis to (as closely as possible)
    superimpose it on top of a reference curve. Then, find the minimum
    distance between each of the points on it and the reference curve.
    Finally return the sum of these distances to give a total distance.

    Uses superposition and min_distance, see those methods for caveats and
    details
    """
    super_imposed = superposition(ref_curve,curve)
    distances = min_distance(*super_imposed)
    print("dist",distances)
    return pow(sum(distances),0.5)

if __name__ == "__main__":
    print(procrustes_distance(numpy.array([[1,1,1],[2,2,2],[3,3,3],[4,4,4]]),numpy.array([[1,1,1],[3,3,3]])))
