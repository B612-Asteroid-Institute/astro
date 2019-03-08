"""
    transform.py

    Sources:
    Vallado, David - Fundamentals of Astrodynamics and Applications
"""

import numpy as np
from math import sin, cos, acos, sqrt

TOL = 1.0e-10
TWO_PI = 2 * np.pi

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

    def cartesian_to_classical(self, r, v):
        """ Convert cartesian elements to classical elements

        Args:
            r (list or array) - cartesian position (km)
            v (list or array) - cartesian velocity (km/s)

        Returns:
            a (float) = semi-major axis (km)
            e (float) = eccentricity
            i (float) = inclination (0-pi rad)
            node (float) = right ascension of ascending node (0-2pi rad)
            w (float) = argument of periapsis (0-2pi rad)
            TA (float) = true anomaly (0-2pi rad)
        """

        # Ensure r, v are numpy arrays
        r = np.array(r)
        v = np.array(v)

        # Get magnitudes of r, v
        mag_r = np.linalg.norm(r)
        mag_v = np.linalg.norm(v)

        # Find h, n, and e vectors
        hbar = np.cross(r, v)  # angular momentum vector
        mag_h = np.linalg.norm(hbar)

        nbar = np.array([0.0] * 3)  # line of nodes
        ebar = np.array([0.0] * 3)  # eccentricity vector

        if mag_h > TOL:

            # Define line of nodes
            nbar[0] = -hbar[1]
            nbar[1] = hbar[0]
            mag_n = np.linalg.norm(nbar)
            c1 = mag_v ** 2 - self._mu / mag_r

            # Get eccentricity
            for i in range(3):
                ebar[i] = (c1 * r[i] - (r.dot(v)) * v[i]) / self._mu

            e = np.linalg.norm(ebar)

            # Get semi-major axis
            a = 0.5 * mag_v ** 2 - self._mu / mag_r
            if abs(a) > TOL:
                a = -self._mu / (2 * a)
            else:
                a = float('nan')

            # Get inclination
            i = acos(hbar[2] / mag_h)

            # Get right ascension of ascending node
            if mag_n > TOL:
                temp = nbar[0] / mag_n
                if abs(temp) > 1.:
                    temp = np.sign(temp)
                node = acos(temp)
                if nbar[2] <= 0.:
                    node = TWO_PI - node
            else:
                node = float('nan')

            # Get argument of periapsis and true anomaly
            if e > TOL:
                w = acos((nbar.dot(ebar)) / (mag_n * e))
                TA = acos((ebar.dot(r)) / (e * mag_r))
                if ebar[2] < 0.:
                    w = TWO_PI - w
                if r.dot(v) < 0.:
                    TA = TWO_PI - TA
            else:
                w = float('nan')
                TA = float('nan')
        else:
            a = float('nan')
            e = float('nan')
            i = float('nan')
            node = float('nan')
            w = float('nan')
            TA = float('nan')

        return a, e, i, node, w, TA