"""
Demo to convert orbital elements to cartesian and back
"""

from astro import Transform
import numpy as np

d2r = np.pi/180.
r2d = 180./np.pi

# Test orbital parameters
a = 1.52854101  # SMA, AU
e = 0.51809438
i = 4.42575336 * d2r
node = 304.98756885 * d2r
w = 281.25948665 * d2r
TA = 39.59664935 * d2r

# Convert SMA to km
AU = 149597870.0     # conversion from AU to km
a = a * AU

# Instantiate class and convert to cartesian
trans = Transform()
r, v = trans.classical_to_cartesian(a, e, i, node, w, TA)

print("Position: %s km" % r)
print("Velocity: %s km/s" % v)

# Convert back to orbital elements
a_new, e_new, i_new, node_new, w_new, TA_new = trans.cartesian_to_classical(r, v)
print("SMA: %s AU" % (a_new/AU))
print("Ecc: %s" % e_new)
print("Inc: %s deg" % (i_new * r2d))
print("RAAN: %s deg" % (node_new * r2d))
print("Arg of per: %s deg" % (w * r2d))
print("True anomaly: %s deg" % (TA * r2d))