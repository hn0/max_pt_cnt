"""
	Collection of useful geographic static functions

	Author:  Hrvoje Novosel<hrvojedotnovosel@gmail.com>
	Created: Thu Aug 18 23:43:06 CEST 2016
"""

import math

class geoclass:

	earth_radius = 6371000

	@staticmethod
	def haversine_distance(lon1, lat1, lon2, lat2):
		"""
		"""
		lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
		
		dlon = lon2 - lon1 
		dlat = lat2 - lat1 
		a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
		c = 2 * math.asin(math.sqrt(a)) 
		return c * geoclass.earth_radius
