from scipy import spatial
import numpy as np

def find(array, value):

	"""

	Find the closest point to "value" in the array "array"
	"array" should be an array with size N by 2 (column 0 = lon, column 1 = lat)
	Carine G. R. Costa, 2017

	"""

	array[spatial.KDTree(array).query(value)[1]]
	distance,index = spatial.KDTree(array).query(value)

	return distance, index
