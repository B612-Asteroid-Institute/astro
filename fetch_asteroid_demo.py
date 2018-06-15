"""
Demo to fetch asteroid info from JPL Horizons database
"""

from astro import AsteroidDB
from datetime import datetime

str_format = '%Y-%m-%dT%H:%M:%S.%f'

obj_id = '2017 DR109'
epoch_str = '2018-03-23T00:00:00.000'
end_time_str = '2018-03-30T00:00:00.000'
epoch = datetime.strptime(epoch_str, str_format)
end_time = datetime.strptime(end_time_str, str_format)

# Instantiate class
asteroid = AsteroidDB()
#asteroid.set_step_size(30, 'day')
rv = asteroid.get_cartesian_states(obj_id, epoch, end_time)

for state in rv:
    print(state)
