from config import Configuration
from libraries.calculator import Calculator
from libraries.observatory import Observatory
from libraries.utils import Utils
from astroplan import Observer
from astropy import units as u
from astropy.coordinates import EarthLocation, get_sun, SkyCoord
from astropy.time import Time
from datetime import timedelta
from ligo.skymap.postprocess import crossmatch
import numpy as np
import os
import pandas as pd
import logging
logging.getLogger('healpy').setLevel(logging.WARNING)

class Priority:

	@staticmethod
	def field_generator(field_sep):
		''' This function will generate survey fields based on the provided separation of each field.

		:parameter field_sep - The desired field separation required between each field center [deg]

		:return field_list - The field list is both saved to file and returned as a Pandas DataFrame

		'''

		declination_limit = Observatory.get_declination_limit(Configuration.OBSERVATORY)
		field_size = Observatory.get_field_size(Configuration.OBSERVATORY)
		latitude = Observatory.get_latitude(Configuration.OBSERVATORY)

		survey_fields_path = Configuration.ANALYSIS_SURVEY_DIRECTORY + Configuration.OBSERVATORY + '_fields.cat'

		# set up variables necessary for geometry
		deg_to_rad = np.pi / 180.

		gal_idx = 0
		tot_idx = 0

		if not os.path.isfile(survey_fields_path):
			Utils.log('Now generating survey fields...', 'info')

			if latitude < 0e0:
				# set up the field data frames and start with the first field
				c = SkyCoord(ra=0e0, dec=-90e0, frame='icrs', unit='deg')
				field_list = pd.DataFrame(data=[['00.000', 0e0, -90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', 300., 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# set up the number of declination strips
				field_number = int(np.ceil((90 + declination_limit) / field_sep))

				# separation of the fields in declination
				declination_strips = -90 + np.arange(0, field_number) * field_sep

			else:
				# set up the field data frames and start with the first field
				c = SkyCoord(ra=0e0, dec=90e0, frame='icrs', unit='deg')
				field_list = pd.DataFrame(data=[['00.000', 0e0, 90e0, c.galactic.l.to_value(), c.galactic.b.to_value(), 'main_survey', 300., 1, 0, 0, 0, 0]], columns=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

				# set up the number of declination strips
				field_number = int(np.ceil((90 - declination_limit) / field_sep))

				# separation of the fields in declination
				declination_strips = 90 - np.arange(0, field_number) * field_sep

			# now loop through and generate the fields
			eo = 0
			for idx in range(1, field_number):
				nfr = np.ceil(360. * np.cos(declination_strips[idx] * deg_to_rad) / field_sep)
				ra_sep = 360. / nfr

				if eo == 1:
					ra_off = 0.5 * ra_sep
					eo = 0
				else:
					ra_off = 0.
					eo = 1

				for idy in range(0, int(nfr)):
					# set up the field name with the first hex field
					if len(hex(idx).split('x')[1]) == 1:
						field_1 = '0' + hex(idx).split('x')[1]
					else:
						field_1 = hex(idx).split('x')[1]

					# set up the second hex field
					if len(hex(idy).split('x')[1]) == 1:
						field_2 = '00' + hex(idy).split('x')[1]
					elif len(hex(idy).split('x')[1]) == 2:
						field_2 = '0' + hex(idy).split('x')[1]
					else:
						field_2 = hex(idy).split('x')[1]

					# get the galactic coordinates of the field
					c_idy = SkyCoord(ra=idy * ra_sep + ra_off, dec=declination_strips[idx], frame='icrs', unit='deg')
					if (c_idy.galactic.b.to_value() < Configuration.GALACTIC_PLANE) and (c_idy.galactic.b.to_value() > -1 * Configuration.GALACTIC_PLANE):
						moon_phase = 1
						gal_idx += 1
					else:
						moon_phase = 0

					# set up the series for appending
					field = pd.Series(data=[field_1 + '.' + field_2, idy * ra_sep + ra_off, declination_strips[idx], c_idy.galactic.l.to_value(), c_idy.galactic.b.to_value(), 'main_survey', Configuration.EXPOSURE_TIME, 1, 0, 0, 0, moon_phase], index=['field_id', 'ra', 'dec', 'l', 'b', 'program', 'exposure_time', 'cadence', 'ephemeris', 'period', 'observations', 'moon_phase'])

					# append the series
					field_list.loc[tot_idx] = field
					tot_idx += 1

			field_list.to_csv(survey_fields_path, sep=' ', header=True, index=False, float_format='%.3f')

			print('Total fields:', tot_idx)
			print('Galactic fields:', gal_idx)
			print('Number of declination strips:', field_number)
			print('Declination strips:', declination_strips)

		else:
			# if the file exists already, then just read the field list in
			Utils.log('Reading extant survey fields.', 'info')
			field_list = pd.read_csv(survey_fields_path, header=0, sep=' ')

		return field_list

	@staticmethod
	def find_field(survey_fields, ra, dec, moon_phase=2, ang_extend=0, program='main_survey'):
		''' This function will select the appropriate survey field based on the RA/Dec of the target you are interested in observing.

		:parameter survey_fields - The set of survey fields
		:parameter ra - The right ascension of the target [deg]
		:parameter dec - The declination of the target [deg]
		:parameter moon_phase - Update this to be 0 if you want dark time, 1 for bright time, or 2 for any time
		:parameter ang_extend - If this is an extended object, provide its angular extension to grab more than one survey field (assumes circular extension).
		:parameter program - If the input is something other than 'main_survey', the program is overwritten

		:return field - The survey field name to be observed for the specific science

		'''

		# get the angular distance between the given position and the survey fields
		ang_dist = survey_fields.apply(lambda x: Calculator.angular_distance(float(x.ra), float(x.dec), ra, dec), axis=1)

		# find the closest survey fields for our objects
		if ang_extend != 0:

			# if an angular extent is given, find the fields within that extent
			if len(survey_fields[ang_dist < ang_extend]) > 0:

				# all fields selected where the field center is within the angular extent
				fields = survey_fields[ang_dist < ang_extend].copy() # copies the new fields
				ang_dists = ang_dist[ang_dist < ang_extend] # gets their angular distances

				# updates the program name to NOT be the main survey
				if program != 'main_survey':
					fields['program'] = program # updates the program name
					fields['moon_phase'] = moon_phase # update when to observe during the lunar cycle

			else:
				fields = 'No field found.'

		else:
			# if no angular extent is provided, then just get the closest field within the survey radius
			fov = Observatory.get_fov(Configuration.OBSERVATORY)

			if np.min(ang_dist) <= fov / 2:

				# pull out the field row
				fields = survey_fields.loc[np.argmin(ang_dist)].copy() # copies the new fields
				ang_dists = np.argmin(ang_dist) # gets their angular distances

				if program != 'main_survey':
					fields['program'] = program # updates the program name
					fields['moon_phase'] = moon_phase # update when to observe during the lunar cycle

			else:
				fields = 'No field found.'

		return fields, ang_dists

	@staticmethod
	def generate_targets(detection_time=None):

		declination_limit = Observatory.get_declination_limit(Configuration.OBSERVATORY)
		latitude = Observatory.get_latitude(Configuration.OBSERVATORY)

		if Configuration.CATALOG == 'glade24':

			Utils.log('Reading GLADE 2.4 catalog.', 'info')

			catalog = Configuration.CATALOGS_DIRECTORY + 'GLADE_2.4.txt'
			delimiter = ' '
			usecols = [1, 6, 7, 8]
			header = None
			names = ['GWGC', 'ra', 'dec', 'dist']
			low_memory = False

		elif Configuration.CATALOG == 'glade+':

			Utils.log('Reading GLADE+ catalog.', 'info')

			catalog = Configuration.CATALOGS_DIRECTORY + 'GLADE+.txt'
			delimiter = ' '
			usecols = [2, 8, 9, 32]
			header = None
			names = ['GWGC', 'ra', 'dec', 'dist']
			low_memory = False

		df = pd.read_csv(catalog, delimiter=delimiter, usecols=usecols, header=header, names=names, low_memory=low_memory)

		Utils.log('Slicing catalog.', 'info')

		# drop NaNs from catalog
		df = df.dropna()

		# slice catalog by declination
		lim_dec = df.dec > declination_limit
		new_df = df[lim_dec]

		# slice catalog by right ascension
		right_now = detection_time or Time.now()
		right_now = Time(right_now)

		sun = get_sun(right_now)

		horizon = -15 * u.degree
		min_height = 30 * u.degree

		alpha_obs_min = (sun.ra - horizon + min_height).degree
		alpha_obs_max = (sun.ra + horizon - min_height).degree

		circum = 90.0 - abs(new_df.dec) < abs(latitude)
		alfa_min = new_df.ra > float(alpha_obs_min)
		alfa_max = new_df.ra <= float(alpha_obs_max)

		case_1 = (alfa_min & alfa_max) | circum
		case_2 = (alfa_min | alfa_max) | circum

		if alpha_obs_max > alpha_obs_min:
			final_df = new_df[case_1]
		else:
			final_df = new_df[case_2]

		Utils.log('Catalog slicing complete.', 'info')

		return final_df
	
	@staticmethod
	def night_field_selector(survey_fields, commissioning_fields, field_select='airmass'):
		''' This function will select fields for tonight based on the current UTC, moon, alert, etc.

		:parameter survey_fields - The full set of survey fields with current observing history
		:parameter commissioning_fields - The full set of commissioning fields to choose from
		:parameter field_select - This should be 'airmass' or 'observations' depending on how you want to prioritize final field selections

		:return obs_list - The list of fields to observe tonight (in order)

		'''

		# merge the lists together to make a giant list selection
		#toros_fields = pd.concat([survey_fields, commissioning_fields]).reset_index(drop=True)
		toros_fields = commissioning_fields.copy()

		# set up the observatory location parameters
		longitude = Observatory.get_longitude(Configuration.OBSERVATORY)
		latitude = Observatory.get_latitude(Configuration.OBSERVATORY)
		elevation = Observatory.get_elevation(Configuration.OBSERVATORY)
		location = EarthLocation.from_geodetic(longitude*u.deg, latitude*u.deg, elevation*u.m)

		# for the observing position to be at the toros location
		toros = Observer(location=location, name='TOROS', timezone='America/Argentina/Salta')
		field_c = SkyCoord(ra=toros_fields.ra.to_list(), dec=toros_fields.dec.to_list(), frame='icrs', unit='deg')

		# get the time for the evening hours
		jd_start = toros.tonight(horizon=-6*u.deg)[0]
		jd_end = toros.tonight(horizon=-6*u.deg)[1]

		# get the moon information
		if toros.moon_illumination(jd_start) >= 0.5:
			moon_bright = 1
		else:
			moon_bright = 0

		# set up initial time-keeping variable
		time_of_exp = jd_start

		# set up exposure count time-keeper variable
		n_exposure = 0

		# opens the observation plan for the night
		time_of_exp.format = 'jd'
		UTC = -3
		date_today = Time((time_of_exp.value + UTC / 24.), format='jd').iso.split(' ')[0]
		f = open(Configuration.ANALYSIS_DIRECTORY + 'nightly_observing_list_' + date_today + '.txt', 'w')

		# keep looping through simulated exposures until you get to the end of the night
		while time_of_exp < jd_end:

			# get the alt/az for the fields and the moon
			moon = toros.moon_altaz(time=time_of_exp)
			targets_exp = toros.altaz(time_of_exp, field_c)

			# make sure that the moon is above the horizon before settling on dark time
			if moon.alt.value > 0:
				watch_out_for_moon = 1
			else:
				watch_out_for_moon = 0

			# update the hour angles for the expected exposure time
			toros_fields['field_ha_exp'] = toros.target_hour_angle(time_of_exp, field_c).value
			toros_fields['airmass_exp'] = targets_exp.secz.value
			toros_fields['altitude_exp'] = targets_exp.alt.value
			toros_fields['azimuth_exp'] = targets_exp.az.value
			toros_fields['horizon_exp'] = toros.target_is_up(time_of_exp, field_c)

			# get the possible fields based on observing cadences
			possible_fields = toros_fields[((toros_fields.airmass_exp >= 1) & (toros_fields.airmass_exp <= 3)) &
											((toros_fields.moon_phase == moon_bright) |
											(toros_fields.moon_phase == 2))].sort_values(by=['dec'], ascending=True)

			'''
			# get the possible fields based on observing cadences
			possible_fields = toros_fields[((toros_fields.airmass_exp >= 1) & (toros_fields.airmass_exp <= 3)) &
											~(toros_fields.field_ha_exp.between(3, 21)) &
											((toros_fields.moon_phase == moon_bright) |
											(toros_fields.moon_phase == 2))].sort_values(by=['dec'], ascending=True)
			'''

			# look for alerts, then science fields, then survey fields, then high cadence fields
			for prio in Configuration.FIELD_PRIORITY:

				# get all the fields in your priority list
				field_chk = possible_fields[possible_fields['program'] == prio]

				#if not field_chk.empty:
				#	print(f'time_of_exp: {time_of_exp}\n {field_chk[['field_id', 'airmass_exp']]}\n')

				# you have to have a field in the priority list
				if len(field_chk) >= 1:

					# get the alt and az of the positions of the fields, assuming the moon is up
					if watch_out_for_moon == 1:
						# the moon is above the horizon, so watch out
						field_chk['moon_dist'] = field_chk.apply(lambda x: Calculator.angular_distance(x.azimuth_exp, x.altitude_exp, moon.az.value, moon.alt.value), axis=1)
					else:
						# the moon is below the horizon, so don't worry
						field_chk['moon_dist'] = 999

					# select fields based on the distance to the moon
					field_chk = field_chk[field_chk['moon_dist'] > Configuration.MOON_DISTANCE]
					if len(field_chk[field_chk['moon_dist'] > Configuration.MOON_DISTANCE]) == 0:
						Utils.log('Oops! The Moon is in your field, observe something else.', 'info')
						break

					# select the field with the lowest airmass to observe
					if field_select == 'airmass':
						field_chk = field_chk[field_chk.airmass_exp == field_chk.airmass_exp.min()]
					else:
						field_chk = field_chk[field_chk.observations == field_chk.observations.min()]

					# get the index of the lowest airmass field
					field = field_chk.index[0]

					# add an observation count to it
					toros_fields.loc[field, 'observations'] += 1

					# increase your time step by the time of the observation
					next_exp = time_of_exp + timedelta(seconds=int(toros_fields.loc[field, 'exposure_time'] + Configuration.OVERHEAD + Configuration.READ_TIME))

					# update time_of_exp to be the next exposure time
					time_of_exp = next_exp 

					# set up your initial observing list
					if n_exposure == 0:
						# if first exposure, then just copy the data
						obs_list = pd.DataFrame([field_chk.iloc[0]]).reset_index(drop=True)
					else:
						# if it is n+1 exposures, concatenate
						obs_list = pd.concat([obs_list, field_chk], axis=0).reset_index(drop=True)

					# for readability get the local time
					obs_time = Time((time_of_exp.value + UTC / 24.), format='jd')

					# write out necessary information to log file
					if n_exposure == 0:
						header = 'date time field_id ra dec exp_time airmass moon_dist program ha alt\n'
						f.write(header)

					line = (obs_time.iso + ' ' + str(field_chk.field_id.values[0]) + ' ' + \
							str(field_chk.ra.values[0]) + ' ' + str(field_chk.dec.values[0]) + ' ' + \
							str(int(field_chk.exposure_time.values[0])) + ' ' + \
							str(np.around(field_chk.airmass_exp.values[0], decimals=2)) + ' ' + \
							str(np.around(field_chk.moon_dist.values[0], decimals=1)) + ' ' + \
							field_chk.program.values[0] + ' ' + \
							str(np.around(field_chk.field_ha_exp.values[0], decimals=2)) + ' ' + \
							str(np.around(field_chk.altitude_exp.values[0], decimals=2)) + '\n')

					f.write(line)

					# update the exposure count
					n_exposure += 1
					break

				else:
					# increase your time step by the time of the observation
					next_exp = time_of_exp + timedelta(seconds=int(Configuration.OVERHEAD + Configuration.READ_TIME))
					time_of_exp = next_exp # update time_of_exp to be the next exposure time

		f.close()

		return obs_list

	@staticmethod
	def sort_field_skymap(survey_fields, skymap, event_name):

		# copy the data frame
		selected_fields = survey_fields.copy().reset_index(drop=True)

		# get the coordinate values for each survey field in the appropriate units
		field_coords = SkyCoord(selected_fields.ra, selected_fields.dec, unit=u.deg)

		# crossmatch the survey fields with the sky map probabilities
		cross_match = crossmatch(skymap, field_coords)

		# move through each survey field and get the integrated area within the field area
		selected_fields['prob'] = cross_match.probdensity

		# sort all survey fields based on the probability strip
		selected_fields = selected_fields.sort_values(by='prob', ascending=False).copy().reset_index(drop=True)
		selected_fields['program'] = 'lvc_alert_' + event_name

		return selected_fields

	@staticmethod
	def sort_galaxy_skymap(catalog, skymap, event_name):

		# copy the catalog
		selected_galaxies = catalog.copy().reset_index(drop=True)

		# get the coordinate values for each galaxy in the appropriate units
		galaxy_coords = SkyCoord(selected_galaxies.ra, selected_galaxies.dec, unit=u.deg)

		# now crossmatch the galaxies with the sky map probabilities
		cross_match = crossmatch(skymap, galaxy_coords)

		# move through each galaxy and get the integrated area within the galaxy area
		selected_galaxies['prob'] = cross_match.probdensity

		# sort all galaxies based on the probability strip
		selected_galaxies = selected_galaxies.sort_values(by='prob', ascending=False).copy().reset_index(drop=True)

		# limit to top N galaxies
		final_galaxies = selected_galaxies.head(Configuration.GALAXY_LIMIT)

		return final_galaxies