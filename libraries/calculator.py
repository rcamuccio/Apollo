from config import Configuration
from libraries.utils import Utils
import numpy as np

class Calculator:

	@staticmethod
	def angular_distance(ra1, dec1, ra2, dec2):
		''' This function will determine the angular distance between two points.

		:parameter ra1 - The right ascension of point 1
		:parameter dec1 - The declination of point 1
		:parameter ra2 - The right ascension of point 2
		:parameter dec2 - The declination of point 2

		:return ang_dist_deg - The angular distance of the two points in degrees

		'''

		# convert to radians
		ra1_rad = np.deg2rad(ra1)
		dec1_rad = np.deg2rad(dec1)
		ra2_rad = np.deg2rad(ra2)
		dec2_rad = np.deg2rad(dec2)

		# angular distance in radians
		ang_dist_rad = np.arccos((np.sin(dec1_rad) * np.sin(dec2_rad)) + (np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra1_rad - ra2_rad)))

		# angular distance in degrees
		ang_dist_deg = np.rad2deg(ang_dist_rad)

		return ang_dist_deg