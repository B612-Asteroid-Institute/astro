"""
Demo
"""

from astro import Transform
import numpy as np

# Test orbital parameters
a = 1.52854101  # SMA, AU
e = 0.51809438
i = 4.42575336 * np.pi/180
node = 304.98756885 * np.pi/180
w = 281.25948665 * np.pi/180
TA = 39.59664935 * np.pi/180

# Convert SMA to km
AU = 149597870.0     # conversion from AU to km
a = a * AU

# Instantiate class and convert to cartesian
trans = Transform()
r, v = trans.classical_to_cartesian(a, e, i, node, w, TA)

print("Position: %s km" % r)
print("Velocity: %s km/s" % v)