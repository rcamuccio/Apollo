import numpy as np

class Observatory:

	@staticmethod
	def get_area_primary(name):

		diameter_primary = Observatory.get_diameter_primary(name)

		area_primary = np.pi * (diameter_primary / 2) ** 2

		return area_primary

	@staticmethod
	def get_area_secondary(name):

		diameter_secondary = Observatory.get_diameter_secondary(name)

		area_secondary = np.pi * (diameter_secondary / 2) ** 2

		return area_secondary

	@staticmethod
	def get_bandpass(name):

		if name == 'ctmo':
			bandpass = [65, 149, 133, 149, 280]
		elif name == 'macon':
			bandpass = [147, 141, 147, 147]
		elif name == 'oafa':
			bandpass = [50]

		return bandpass

	@staticmethod
	def get_central_wavelength(name):

		if name == 'ctmo':
			central_wavelength = [352.5, 475.5, 628.5, 769.5, 960.0]
		elif name == 'macon':
			central_wavelength = [473.5, 638.5, 775.5, 922.5]
		elif name == 'oafa':
			central_wavelength = [530.0]

		return central_wavelength

	@staticmethod
	def get_collecting_area(name):

		telescope = Observatory.get_telescope(name)
		area_primary = Observatory.get_area_primary(name)
		area_secondary = Observatory.get_area_secondary(name)

		if telescope == 'reflector':
			collecting_area = area_primary - area_secondary
		elif telescope == 'refractor':
			collecting_area = area_primary

		return collecting_area

	@staticmethod
	def get_declination_limit(name):

		if name == 'ctmo':
			declination_limit = -34.00
		elif name == 'macon':
			declination_limit = 33.84
		elif name == 'oafa':
			declination_limit = 26.66

		return declination_limit

	@staticmethod
	def get_diameter_primary(name):

		if name == 'ctmo':
			diameter_primary = 0.432
		elif name == 'macon':
			diameter_primary = 0.610
		elif name == 'oafa':
			diameter_primary = 0.508

		return diameter_primary

	@staticmethod
	def get_diameter_secondary(name):

		if name == 'ctmo':
			diameter_secondary = 0.190
		elif name == 'macon':
			diameter_secondary = 0.280
		elif name == 'oafa':
			diameter_secondary = 0.0

		return diameter_secondary

	@staticmethod
	def get_etendue(name):

		collecting_area = Observatory.get_collecting_area(name)
		fov = Observatory.get_fov(name)

		etendue = collecting_area * fov

		return etendue

	@staticmethod
	def get_elevation(name):

		if name == 'ctmo':
			elevation = 11.5
		elif name == 'macon':
			elevation = 4650
		elif name == 'oafa':
			elevation = 2420

		return elevation

	@staticmethod
	def get_field_separation(name):

		if name == 'ctmo':
			field_separation = None
		elif name == 'macon':
			field_separation = 1.19
		elif name == 'oafa':
			field_separation = 1.2648678231857113

		return field_separation

	@staticmethod
	def get_field_size(name):

		pixel_scale = Observatory.get_pixel_scale(name)
		pixel_number = Observatory.get_pixel_number(name)

		field_size = pixel_scale * pixel_number / 3600.

		return field_size

	@staticmethod
	def get_focal_length(name):

		if name == 'ctmo':
			focal_length = 2.939
		elif name == 'macon':
			focal_length = 3.974
		elif name == 'oafa':
			focal_length = 3.7

		return focal_length

	@staticmethod
	def get_fov(name):

		field_size = Observatory.get_field_size(name)

		fov = field_size ** 2

		return fov

	@staticmethod
	def get_gain(name):

		if name == 'ctmo':
			gain = 1.385
		elif name == 'macon':
			gain = 2.18
		elif name == 'oafa':
			gain = 2.18

		return gain

	@staticmethod
	def get_latitude(name):

		if name == 'ctmo':
			latitude = 25.995789
		elif name == 'macon':
			latitude = -24.62055556
		elif name == 'oafa':
			latitude = -31.8023

		return latitude

	@staticmethod
	def get_longitude(name):

		if name == 'ctmo':
			longitude = -97.568956
		elif name == 'macon':
			longitude = -67.32833333
		elif name == 'oafa':
			longitude = -69.3265

		return longitude

	@staticmethod
	def get_pixel_number(name):

		if name == 'ctmo':
			pixel_number = 4096
		elif name == 'macon':
			pixel_number = 10560
		elif name == 'oafa':
			pixel_number = 10560

		return pixel_number

	@staticmethod
	def get_pixel_scale(name):

		if name == 'ctmo':
			pixel_scale = 0.6305
		elif name == 'macon':
			pixel_scale = 0.468
		elif name == 'oafa':
			pixel_scale = 0.4959

		return pixel_scale

	@staticmethod
	def get_qe_atmosphere(name):
		
		if name == 'ctmo':
			qe_atmosphere = [0.8, 0.8, 0.9, 0.9, 0.9]
		elif name == 'macon':
			qe_atmosphere = [0.8, 0.9, 0.9, 0.9]
		elif name == 'oafa':
			qe_atmosphere = [0.9]

		return qe_atmosphere

	@staticmethod
	def get_qe_ccd(name):

		if name == 'ctmo':
			qe_ccd = 0.9
		elif name == 'macon':
			qe_ccd = 0.9
		elif name == 'oafa':
			qe_ccd = 0.9

		return qe_ccd

	@staticmethod
	def get_qe_filter(name):

		if name == 'ctmo':
			qe_filter = 0.9
		elif name == 'macon':
			qe_filter = 0.9
		elif name == 'oafa':
			qe_filter = 1.0

		return qe_filter

	@staticmethod
	def get_qe_primary(name):

		if name == 'ctmo':
			qe_primary = 0.96
		elif name == 'macon':
			qe_primary = 0.96
		elif name == 'oafa':
			qe_primary = 0.92

		return qe_primary

	@staticmethod
	def get_qe_secondary(name):

		if name == 'ctmo':
			qe_secondary = 0.96
		elif name == 'macon':
			qe_secondary = 0.96
		elif name == 'oafa':
			qe_secondary = 1.0

		return qe_secondary

	@staticmethod
	def get_readout_noise(name):

		if name == 'ctmo':
			readout_noise = 15.8
		elif name == 'macon':
			readout_noise = 5.0
		elif name == 'oafa':
			readout_noise = 5.0

		return readout_noise

	@staticmethod
	def get_seeing(name):

		if name == 'ctmo':
			seeing = 4.0
		elif name == 'macon':
			seeing = 0.93
		elif name == 'oafa':
			seeing = 2.5

		return seeing

	@staticmethod
	def get_sky(name):

		if name == 'ctmo':
			sky = [19.81, 19.81, 19.81, 19.81, 19.81]
		elif name == 'macon':
			sky = [22.1, 21.1, 20.1, 18.7]
		elif name == 'oafa':
			sky = [21.96]

		return sky

	@staticmethod
	def get_telescope(name):

		if name == 'ctmo':
			telescope = 'reflector'
		elif name == 'macon':
			telescope = 'reflector'
		elif name == 'oafa':
			telescope = 'refractor'

		return telescope

	@staticmethod
	def get_total_throughput(name):

		qe_atmosphere = Observatory.get_qe_atmosphere(name)
		qe_ccd = Observatory.get_qe_ccd(name)
		qe_filter = Observatory.get_qe_filter(name)
		qe_primary = Observatory.get_qe_primary(name)
		qe_secondary = Observatory.get_qe_secondary(name)

		total_throughput = qe_atmosphere * qe_ccd * qe_filter * qe_primary * qe_secondary

		return total_throughput