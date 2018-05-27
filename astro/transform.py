"""
    transform.py

    Sources:
    Vallado, David - Fundamentals of Astrodynamics and Applications
"""

import numpy as np
from math import sin, cos, sqrt


class Transform(object):
    """
    Does all transformations, including between coordinate systems
    TODO: add unit tests
    """
    def __init__(self):
        """Initialize attributes
        """
        self._mu = 132712440018.0  # sun gravitational parameter, km^3/s^2

    def set_gravitational_parameter(self, mu):
        """
        Set the gravitational parameter of interest

        Args:
            mu (float) - standard gravitational parameter, km^3/s^2

        Returns:
            None
        """
        self._mu = mu

    def Rotx(self, angle):
        """
        Returns the transformation corresponding to the rotation about the x-axis
        by an angle.

        Args:
            angle (float) - angle in radians

        Returns:
            R (numpy matrix) - rotation matrix
        """
        cosangle = cos(angle)
        sinangle = sin(angle)

        # Return rotation matrix
        R = np.array([[1, 0, 0], \
                      [0, cosangle, sinangle], \
                      [0, -sinangle, cosangle]])
        return R

    def Roty(self, angle):
        """
        Returns the transformation corresponding to the rotation about the y-axis
        by an angle.

        Args:
            angle (float) - angle in radians

        Returns:
            R (numpy matrix) - rotation matrix
        """
        cosangle = cos(angle)
        sinangle = sin(angle)

        # Return rotation matrix:
        R = np.array([[cosangle, 0, -sinangle], \
                      [0, 1, 0], \
                      [sinangle, 0, cosangle]])
        return R

    def Rotz(self, angle):
        """
        Returns the transformation corresponding to the rotation about the z-axis
        by an angle.

        Args:
            angle (float) - angle in radians

        Returns:
            R (numpy matrix) - rotation matrix
        """
        cosangle = cos(angle)
        sinangle = sin(angle)

        # Return rotation matrix:
        R = np.array([[cosangle, sinangle, 0], \
                      [-sinangle, cosangle, 0], \
                      [0, 0, 1]])
        return R

    def classical_to_cartesian(self, a, e, i, node, w, TA):
        """ Convert classical elements to cartesian elements

        Args:
            a (float) = semi-major axis (km)
            e (float) = eccentricity
            i (float) = inclination (0-pi rad)
            node (float) = right ascension of ascending node (0-2pi rad)
            w (float) = argument of periapsis (0-2pi rad)
            TA (float) = true anomaly (0-2pi rad)

        Returns:
            r (numpy array) - cartesian position (km)
            v (numpy array) - cartesian velocity (km/s)
        """
        TOL = 1.0e-10  # tolerance for approximation

        # Convert
        p = a * (1 - e ** 2)  # semi-latus rectum, km

        # Determine what type of orbit is involved and set up angles for special cases
        if e < TOL:
            # argument of periapsis and RAAN are 0 for circular orbits
            w = 0.0
            # Circular equatorial
            if (i < TOL) or (abs(i - np.pi) < TOL):
                node = 0.0
        # Elliptical equatorial
        elif (i < TOL) or (abs(i - np.pi) < TOL):
            node = 0.0

        # Position and velocity vectors in perifocal coordinate system
        r_pqw = [0.0] * 3
        v_pqw = [0.0] * 3
        temp = p / (1.0 + e * cos(TA))
        r_pqw[0] = temp * cos(TA)
        r_pqw[1] = temp * sin(TA)

        if abs(p) < 0.0001:
            p = 0.0001

        v_pqw[0] = -sin(TA) * sqrt(self._mu) / sqrt(p)
        v_pqw[1] = (e + cos(TA)) * sqrt(self._mu) / sqrt(p)

        # Transform to cartesian
        tempvec = self.Rotz(-w).dot(r_pqw)
        tempvec = self.Rotx(-i).dot(tempvec)
        r = self.Rotz(-node).dot(tempvec)

        tempvec = self.Rotz(-w).dot(v_pqw)
        tempvec = self.Rotx(-i).dot(tempvec)
        v = self.Rotz(-node).dot(tempvec)

        return r, v