from config import Configuration
from libraries.alerts import Alerts
from libraries.observatory import Observatory
from libraries.priority import Priority
from libraries.utils import Utils
from astropy import units as u
from astropy.table import QTable, Table
from base64 import b64decode
from io import BytesIO
import astropy_healpix as ah
import json
import numpy as np
import os

class Processor:

	@staticmethod
	def decode_classic_notice(value, verbose=False):

		contents = value.decode('utf-8').split('\n')

		keys = []
		items = []

		for content in contents:
			if len(content.split(': ')) > 1:
				key = content.split(': ')[0].strip()
				keys.append(key)

				item = content.split(': ')[1].strip()
				items.append(item)

				if verbose:
					print(key, item)

		params = dict(zip(keys, items))

		if verbose:
			print('----------')

		return params

	@staticmethod
	def process_igwn_gwalert(value):

		record = json.loads(value)

		superevent_id = record['superevent_id']
		time_created = record['time_created']
		url = record['urls']['gracedb']
		alert_type = record['alert_type']

		Utils.log('Event ' + str(superevent_id) + ' generated at time ' + str(time_created) + '.', 'info')
		Utils.log('GraceDB URL: ' + str(url) + '.', 'info')

		event_tag = superevent_id[0]
		full_type = event_tag + '-' + alert_type
		is_mock = False

		if (event_tag == 'M') or (event_tag == 'T'):
			Utils.log('This is a mock event!', 'info')
			is_mock = True
		elif event_tag == 'S':
			Utils.log('This is a true event!', 'info')
		else:
			Utils.log('This is an unknown event (event_tag == ' + str(event_tag) + ').', 'info')

		field_path = Configuration.ANALYSIS_GW_FIELDS_DIRECTORY + superevent_id + '-skymap-fields.txt'
		galaxy_path = Configuration.ANALYSIS_GW_GALAXIES_DIRECTORY + superevent_id + '-skymap-galaxies.txt'
		json_path = Configuration.ALERTS_GW_JSON_DIRECTORY + superevent_id + '.json'
		moc_path = Configuration.ALERTS_GW_MOC_DIRECTORY + superevent_id + '-90percent.moc.fits'
		param_path = Configuration.ALERTS_GW_PARAM_DIRECTORY + superevent_id + '.txt'
		skymap_path = Configuration.ALERTS_GW_SKYMAP_DIRECTORY + superevent_id + '-skymap.fits'

		# filter alert type
		if alert_type == 'RETRACTION':
			Utils.log('Alert for ' + superevent_id + ' was retracted!', 'info')

			# send team alert
			if is_mock:
				Utils.log('Passing alert due to being a mock event.', 'info')
			else:
				Utils.log('Alerting team of retraction.', 'info')
				Alerts.alert_team(alert_channel='igwn_gwalert', alert_type=full_type, event_name=superevent_id)

			# delete fields
			if os.path.isfile(field_path):
				Utils.log('Deleting fields for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + field_path)
			else:
				Utils.log('No fields for ' + superevent_id + ' to delete.', 'info')

			# delete galaxies
			if os.path.isfile(galaxy_path):
				Utils.log('Deleting galaxies for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + galaxy_path)
			else:
				Utils.log('No galaxies for ' + superevent_id + ' to delete.', 'info')

			# delete json
			if os.path.isfile(json_path):
				Utils.log('Deleting JSON for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + json_path)
			else:
				Utils.log('No JSON for ' + superevent_id + ' to delete.', 'info')

			# delete params
			if os.path.isfile(param_path):
				Utils.log('Deleting parameters for retracted alert ' + superevent_id + '.', 'info')
				os.system('rm ' + param_path)
			else:
				Utils.log('No parameters for ' + superevent_id + ' to delete.', 'info')

			# delete skymap
			if os.path.isfile(skymap_path):
				Utils.log('Deleting skymap for ' + superevent_id + '.', 'info')
				os.system('rm ' + skymap_path)
			else:
				Utils.log('No skymap for ' + superevent_id + ' to delete.', 'info')

			# delete moc
			if os.path.isfile(moc_path):
				Utils.log('Deleting 90 percent credible region for ' + superevent_id + '.', 'info')
				os.system('rm ' + moc_path)
			else:
				Utils.log('No 90 percent credible region for ' + superevent_id + ' to delete.', 'info')

			print()

		elif (alert_type == 'PRELIMINARY') or (alert_type == 'INITIAL') or (alert_type == 'UPDATE'):
			Utils.log('This is a ' + alert_type + ' alert for ' + superevent_id + '!', 'info')

			# send team alert
			if is_mock:
				Utils.log('Passing alert due to being a mock event.', 'info')
			else:
				Utils.log('Alerting team of retraction.', 'info')
				Alerts.alert_team(alert_channel='igwn_gwalert', alert_type=full_type, event_name=superevent_id)

			# save json
			if not os.path.isfile(json_path):
				Utils.log('Saving JSON for ' + superevent_id + '.', 'info')
				with open(json_path, 'w') as outfile:
					json.dump(record, outfile)
			else:
				Utils.log('JSON already saved for ' + superevent_id + '.', 'info')

			# read alert parameters and skymap
			Utils.log('Saving event parameters.', 'info')

			if not os.path.isfile(param_path):
				with open(param_path, 'w') as param_file:

					time = record['event']['time'] # string
					far = record['event']['far'] # float
					significant = record['event']['significant'] # bool

					instrument_string = ''
					instruments = record['event']['instruments'] # list
					for item in instruments:
						if item == instruments[-1]:
							instrument_string += item
						else:
							instrument_string += item + ', '

					pipeline = record['event']['pipeline'] # string
					search = record['event']['search'] # string

					param_file.write('Time: ' + str(time) + '\n')
					param_file.write('FAR: ' + str(far) + '\n')
					param_file.write('Significant: ' + str(significant) + '\n')
					param_file.write('Instruments: ' + instrument_string + '\n')
					param_file.write('Pipeline: ' + pipeline + '\n')
					param_file.write('Search: ' + search + '\n')
			
					group = record['event']['group'] # string

					param_file.write('Group: ' + group + '\n')

					if group == 'Burst':

						Utils.log('This is a burst event!', 'info')

						duration = record['event']['duration'] # 
						central_frequency = record['event']['central_frequency'] # 

						param_file.write('Duration: ' + str(duration) + '\n')
						param_file.write('Central Frequency: ' + str(central_frequency) + '\n')

					elif group == 'CBC':

						Utils.log('This is a merger event!', 'info')

						has_ns = record['event']['properties']['HasNS'] # float
						has_remnant = record['event']['properties']['HasRemnant'] # float
						has_mass_gap = record['event']['properties']['HasMassGap'] # float

						bns_class = record['event']['classification']['BNS'] # float
						nsbh_class = record['event']['classification']['NSBH'] # float
						bbh_class = record['event']['classification']['BBH'] # float
						terra_class = record['event']['classification']['Terrestrial'] # float

						param_file.write('HasNS: ' + str(has_ns) + '\n')
						param_file.write('HasRemnant: ' + str(has_remnant) + '\n')
						param_file.write('HasMassGap: ' + str(has_mass_gap) + '\n')
						param_file.write('BNS Class: ' + str(bns_class) + '\n')
						param_file.write('NSBH Class: ' + str(nsbh_class) + '\n')
						param_file.write('BBH Class: ' + str(bbh_class) + '\n')
						param_file.write('Terra Class: ' + str(terra_class) + '\n')

			else:
				Utils.log('Parameters already saved for ' + superevent_id + '.', 'info')

			Utils.log('Reading skymap.', 'info')
			skymap_str = record.get('event', {}).pop('skymap')

			if skymap_str:

				skymap_bytes = b64decode(skymap_str)
				buffer = BytesIO(skymap_bytes)
				skymap = QTable.read(buffer)

				hdr_object = skymap.meta['OBJECT'] # unique identifier for this event
				hdr_date = skymap.meta['DATE-OBS'] # UTC of observation
				hdr_mjd = skymap.meta['MJD-OBS'] # MJD of observation

				try:
					hdr_distmean = skymap.meta['DISTMEAN'] # posterior mean distance [Mpc]
				except KeyError:
					Utils.log('Could not read DISTMEAN.', 'info')

				try:
					hdr_diststd = skymap.meta['DISTSTD'] # posterior std distance [Mpc]
				except KeyError:
					Utils.log('Could not read DISTSTD.', 'info')
					
				# most probable sky location
				i = np.argmax(skymap['PROBDENSITY'])
				uniq = skymap[i]['UNIQ']

				level, ipix = ah.uniq_to_level_ipix(uniq)
				nside = ah.level_to_nside(level)

				ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
				ra_max = ra.deg
				dec_max = dec.deg

				# sort the pixels of the skymap by descending probability density
				skymap.sort('PROBDENSITY', reverse=True)

				# find the area of each pixel
				level, ipix	= ah.uniq_to_level_ipix(skymap['UNIQ'])
				pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

				# calculate the probability within each pixel: pixel area times probability density
				prob = pixel_area * skymap['PROBDENSITY']

				# calculate the cumulative sum of the probability
				cumprob = np.cumsum(prob)

				# find the pixel for which the probability sums to 0.9
				i = cumprob.searchsorted(0.9)

				# calculate the area of the 90 percent credible region: sum of the areas of the pixels up to that one
				area_90 = pixel_area[:i].sum()
				area_90.to_value(u.deg**2)

				# save skymap to FITS
				if not os.path.isfile(skymap_path):
					Utils.log('Saving skymap for event ' + superevent_id + '.', 'info')
					skymap.write(skymap_path, overwrite=True)

				else:
					Utils.log('Skymap already saved for ' + superevent_id + '.', 'info')

				# save 90 percent credible region to FITS
				if not os.path.isfile(moc_path):
					Utils.log('Saving 90 percent credible region for event ' + superevent_id + '.', 'info')

					# keep only the pixels that are within the 90 percent credible region
					skymap90 = skymap[:i]

					# sort the pixels by their UNIQ pixel index
					skymap90.sort('UNIQ')

					# delete all columns except for the UNIQ column
					skymap90 = skymap90['UNIQ',]
					skymap90.write(moc_path, overwrite=True)

				else:
					Utils.log('MOC 90 percent credible region already saved for ' + superevent_id + '.', 'info')

			# check to see if the field list from the skymap exists
			if not os.path.isfile(field_path):
				Utils.log('Generating field list from skymap for ' + superevent_id + '.', 'info')

				# pull in the survey fields
				sep = Observatory.get_field_separation(Configuration.OBSERVATORY)
				survey_fields = Priority.field_generator(sep)

				# generate a ranked list of fields within the skymap
				fields_prio = Priority.sort_field_skymap(survey_fields, skymap, superevent_id)
				fields_prio.to_csv(field_path, sep=' ', index=False)

			else:
				Utils.log('Field list already generated from skymap for ' + superevent_id + '.', 'info')

			# check to see if the galaxy list from the skymap exists
			if not os.path.isfile(galaxy_path):
				Utils.log('Generating galaxy list from skymap for ' + superevent_id + '.', 'info')

				# pull in the survey fields
				sep = Observatory.get_field_separation(Configuration.OBSERVATORY)
				survey_fields = Priority.field_generator(sep)

				# pull in the galaxy catalog
				catalog = Priority.generate_targets(detection_time=None)

				# generated a ranked list of galaxies within the skymap
				galaxies_prio = Priority.sort_galaxy_skymap(catalog, skymap, superevent_id)
				galaxies_prio.to_csv(galaxy_path, sep=' ', index=False)

			else:
				Utils.log('Galaxy list already generated from skymap for ' + superevent_id + '.', 'info')

			print()

	@staticmethod
	def process_lvk_nu_track_search(value):

		record = json.loads(value)

		alert_type = record['type']
		reference = record['reference']
		ref_id = record['ref_id']
		alert_datetime = record['alert_datetime']
		trigger_time = record['trigger_time']
		observation_start = record['observation_start']
		observation_stop = record['observation_stop']
		observation_livetime = record['observation_livetime']
		pval_generic = record['pval_generic']
		pval_bayesian = record['pval_bayesian']
		n_events_coincident = record['n_events_coincident']

		try:
			coincident_events = record['coincident_events']
		except KeyError:
			coincident_events = None

		try:
			most_probable_direction = record['most_probable_direction']
		except KeyError:
			most_probable_direction = None

		try:
			flux_sensitivity = record['neutrino_flux_sensitivity_range']['flux_sensitivity']
		except KeyError:
			flux_sensitivity = None

		try:
			sensitive_energy_range = record['neutrino_flux_sensitivity_range']['sensitive_energy_range']
		except KeyError:
			sensitive_energy_range = None

	@staticmethod
	def process_swift_bat_guano(value):

		record = json.loads(value)

		mission = record['mission']
		instrument = record['instrument']
		messenger = record['messenger']
		record_number = record['record_number']
		alert_datetime = record['alert_datetime']
		alert_tense = record['alert_tense']
		alert_type = record['alert_type']
		trigger_time = record['trigger_time']
		follow_up_event = record['follow_up_event']
		follow_up_type = record['follow_up_type']
		data_archive_page = record['data_archive_page']
		alert_id = record['id']

		Utils.log('Event BAT-GUANO/' + follow_up_event + ' generated at time ' + str(alert_datetime) + '.', 'info')
		Utils.log('Data Archive Page: ' + str(data_archive_page) + '.', 'info')

		# 1 - guano.example.json
		if record_number == 1:
			Utils.log('This is a BAT-GUANO (standard) alert for event ' + follow_up_event + '.', 'info')

			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']
		
		# 2 - guano.loc_map.example.json
		if record_number == 2:
			Utils.log('This is a BAT-GUANO (loc_map) alert for event ' + follow_up_event + '.', 'info')

			healpix_file = record['healpix_file']
			systematic_included = record['systematic_included']
			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']

		# 3 - guano.loc_arc_min.example.json
		elif record_number == 3:
			Utils.log('This is a BAT-GUANO (loc_arc_min) alert for event ' + follow_up_event + '.', 'info')

			ra = record['ra']
			dec = record['dec']
			ra_dec_error = record['ra_dec_error']
			containment_probability = record['containment_probability']
			systematic_included = record['systematic_included']
			image_snr = record['image_snr']
			image_duration = record['image_duration']
			image_energy_range = record['image_energy_range']
			rate_snr = record['rate_snr']
			rate_duration = record['rate_duration']
			rate_energy_range = record['rate_energy_range']
			classification = record['classification']
			far = record['far']

		# 4 - guano.retraction
		elif record_number == 4:
			Utils.log('BAT-GUANO alert for event ' + follow_up_event + ' was retracted!', 'info')

		else:
			print('Could not read BAT-GUANO notice.')
			print(topic)
			print(value)

	@staticmethod
	def process_fermi_gbm_alert(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_signif = params['TRIGGER_SIGNIF']
		trigger_dur = params['TRIGGER_DUR']
		e_range = params['E_RANGE']
		algorithm = params['ALGORITHM']
		detectors = params['DETECTORS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_gbm_fin_pos(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		e_range = params['E_RANGE']
		loc_algorithm = params['LOC_ALGORITHM']
		lc_url = params['LC_URL']
		loc_url = params['LOC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']		

	@staticmethod
	def process_fermi_gbm_flt_pos(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten = params['GRB_INTEN']
		data_signif = params['DATA_SIGNIF']
		integ_time = params['INTEG_TIME']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		hard_ratio = params['HARD_RATIO']
		loc_algorithm = params['LOC_ALGORITHM']
		most_likely = params['MOST_LIKELY']
		secmost_likely = params['2nd_MOST_LIKELY']
		try:
			detectors = params['DETECTORS']
		except KeyError:
			detectors = None
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_gbm_gnd_pos(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		data_signif = params['DATA_SIGNIF']
		data_interval = params['DATA_INTERVAL']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		e_range = params['E_RANGE']
		loc_algorithm = params['LOC_ALGORITHM']
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_gbm_pos_test(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten = params['GRB_INTEN']
		data_signif = params['DATA_SIGNIF']
		integ_time = params['INTEG_TIME']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		data_time_scale = params['DATA_TIME_SCALE']
		hard_ratio = params['HARD_RATIO']
		loc_algorithm = params['LOC_ALGORITHM']
		most_likely = params['MOST_LIKELY']
		secmost_likely = params['2nd_MOST_LIKELY']
		detectors = params['DETECTORS']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_gbm_subthresh(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trans_num = params['TRANS_NUM']
		full_id_num = params['FULL_ID_NUM']
		trans_ra = params['TRANS_RA']
		trans_dec = params['TRANS_DEC']
		trans_error = params['TRANS_ERROR']
		trans_duration = params['TRANS_DURATION']
		trans_date = params['TRANS_DATE']
		trans_time = params['TRANS_TIME']
		earth_angle = params['EARTH_ANGLE']
		spectral_class = params['SPECTRAL_CLASS']
		type_class = params['TYPE_CLASS']
		reliability = params['RELIABILITY']
		healpix_url = params['HEALPIX_URL']
		map_url = params['MAP_URL']
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_lat_monitor(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		source_obj = params['SOURCE_OBJ']
		ref_num = params['REF_NUM']
		ra = params['RA']
		dec = params['DEC']
		curr_flux = params['CURR_FLUX']
		base_flux = params['BASE_FLUX']
		significance = params['SIGNIFICANCE']
		time_scale = params['TIME_SCALE']
		energy_band = params['ENERGY_BAND']
		outburst_date = params['OUTBURST_DATE']
		outburst_time = params['OUTBURST_TIME']
		soln_status = params['SOLN_STATUS']
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_lat_pos_test(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		record_num = params['RECORD_NUM']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten_tot = params['GRB_INTEN_TOT']
		grb_inten1 = params['GRB_INTEN1']
		grb_inten2 = params['GRB_INTEN2']
		grb_inten3 = params['GRB_INTEN3']
		grb_inten4 = params['GRB_INTEN4']
		integ_dur = params['INTEG_DUR']
		first_photon = params['FIRST_PHOTON']
		last_photon = params['LAST_PHOTON']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		soln_status = params['SOLN_STATUS']
		burst_id = params['BURST_ID']
		temp_test_stat = params['TEMP_TEST_STAT']
		image_test_stat = params['IMAGE_TEST_STAT']
		loc_quality = params['LOC_QUALITY']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_fermi_pointdir(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		curr_point_ra = params['CURR_POINT_RA']
		curr_point_dec = params['CURR_POINT_DEC']
		curr_date = params['CURR_DATE']
		curr_time = params['CURR_TIME']
		delta_time = params['DELTA_TIME']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		future_ra_dec = params['FUTURE_RA_DEC']
		comments = params['COMMENTS']

	@staticmethod
	def	process_icecube_astrotrack_bronze(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		stream = params['STREAM']
		run_num = params['RUN_NUM']
		event_num = params['EVENT_NUM']
		src_ra = params['SRC_RA']
		src_dec = params['SRC_DEC']
		src_error = params['SRC_ERROR']
		src_error50 = params['SRC_ERROR50']
		discovery_date = params['DISCOVERY_DATE']
		discovery_time = params['DISCOVERY_TIME']
		revision = params['REVISION']
		energy = params['ENERGY']
		signalness = params['SIGNALNESS']
		far = params['FAR']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']
	
	@staticmethod
	def process_swift_actual_pointdir(value):

		params = Processor.decode_classic_notice(value)
	
		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		curr_point_ra = params['CURR_POINT_RA']
		curr_point_dec = params['CURR_POINT_DEC']
		curr_point_roll = params['CURR_POINT_ROLL']
		slew_time = params['SLEW_TIME']
		slew_date = params['SLEW_DATE']
		tgt_num = params['TGT_NUM']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_bat_grb_lc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_index = params['TRIGGER_INDEX']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		delta_time = params['DELTA_TIME']
		trigger_dur = params['TRIGGER_DUR']
		soln_status = params['SOLN_STATUS']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_bat_grb_pos_ack(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten = params['GRB_INTEN']
		trigger_dur = params['TRIGGER_DUR']
		trigger_index = params['TRIGGER_INDEX']
		bkg_inten = params['BKG_INTEN']
		bkg_time = params['BKG_TIME']
		bkg_dur = params['BKG_DUR']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		soln_status = params['SOLN_STATUS']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		merit_params = params['MERIT_PARAMS']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_bat_grb_pos_test(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten = params['GRB_INTEN']
		trigger_dur = params['TRIGGER_DUR']
		trigger_index = params['TRIGGER_INDEX']
		bkg_inten = params['BKG_INTEN']
		bkg_time = params['BKG_TIME']
		bkg_dur = params['BKG_DUR']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		grb_phi = params['GRB_PHI']
		grb_theta = params['GRB_THETA']
		soln_status = params['SOLN_STATUS']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		merit_params = params['MERIT_PARAMS']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_bat_ql_pos(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		rate_signif = params['RATE_SIGNIF']
		soln_status = params['SOLN_STATUS']
		flags = params['FLAGS']
		merit = params['MERIT']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_bat_scaledmap(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		map_date = params['MAP_DATE']
		map_time = params['MAP_TIME']
		trigger_dur = params['TRIGGER_DUR']
		soln_status = params['SOLN_STATUS']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		map_url = params['MAP_URL']
		sun_postn = params['SUN_POSTN']
		moon_postn = params['MOON_POSTN']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_fom_obs(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_index = params['TRIGGER_INDEX']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		flags = params['FLAGS']
		merit = params['MERIT']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_pointdir(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		next_point_ra = params['NEXT_POINT_RA']
		next_point_dec = params['NEXT_POINT_DEC']
		next_point_roll = params['NEXT_POINT_ROLL']
		slew_time = params['SLEW_TIME']
		slew_date = params['SLEW_DATE']
		obs_time = params['OBS_TIME']
		tgt_name = params['TGT_NAME']
		tgt_num = params['TGT_NUM']
		merit = params['MERIT']
		inst_modes = params['INST_MODES']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_sc_slew(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_index = params['TRIGGER_INDEX']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		slew_query = params['SLEW_QUERY']
		wait_time = params['WAIT_TIME']
		obs_time = params['OBS_TIME']
		inst_modes = params['INST_MODES']
		merit = params['MERIT']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_too_fom(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_index = params['TRIGGER_INDEX']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		flags = params['FLAGS']
		merit = params['MERIT']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_too_sc_slew(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_date = params['GRB_DATE']
		grb_time = params['GRB_TIME']
		trigger_index = params['TRIGGER_INDEX']
		rate_signif = params['RATE_SIGNIF']
		image_signif = params['IMAGE_SIGNIF']
		slew_query = params['SLEW_QUERY']
		wait_time = params['WAIT_TIME']
		obs_time = params['OBS_TIME']
		inst_modes = params['INST_MODES']
		merit = params['MERIT']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_uvot_dburst(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		roll = params['ROLL']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		ffilter = params['FILTER']
		exposure_id = params['EXPOSURE_ID']
		x_offset = params['X_OFFSET']
		y_offset = params['Y_OFFSET']
		width = params['WIDTH']
		height = params['HEIGHT']
		x_grb_pos = params['X_GRB_POS']
		y_grb_pos = params['Y_GRB_POS']
		binning_index = params['BINNING_INDEX']
		im_url = params['IM_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_uvot_dburst_proc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		roll = params['ROLL']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		ffilter = params['FILTER']
		exposure_id = params['EXPOSURE_ID']
		x_offset = params['X_OFFSET']
		y_offset = params['Y_OFFSET']
		width = params['WIDTH']
		height = params['HEIGHT']
		x_grb_pos = params['X_GRB_POS']
		y_grb_pos = params['Y_GRB_POS']
		binning_index = params['BINNING_INDEX']
		im_url = params['IM_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_uvot_fchart(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		point_roll = params['POINT_ROLL']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		ffilter = params['FILTER']
		bkg_mean = params['BKG_MEAN']
		n_stars = params['N_STARS']
		x_offset = params['X_OFFSET']
		y_offset = params['Y_OFFSET']
		x_max = params['X_MAX']
		y_max = params['Y_MAX']
		det_thresh = params['DET_THRESH']
		photo_thresh = params['PHOTO_THRESH']
		sl_url = params['SL_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_uvot_fchart_proc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		point_roll = params['POINT_ROLL']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		ffilter = params['FILTER']
		bkg_mean = params['BKG_MEAN']
		n_stars = params['N_STARS']
		x_offset = params['X_OFFSET']
		y_offset = params['Y_OFFSET']
		x_max = params['X_MAX']
		y_max = params['Y_MAX']
		det_thresh = params['DET_THRESH']
		photo_thresh = params['PHOTO_THRESH']
		sl_url = params['SL_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_centroid(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		counts = params['COUNTS']
		std_dev = params['STD_DEV']
		ph2_iter = params['PH2_ITER']
		error_code = params['ERROR_CODE']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_lc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		lc_start_date = params['LC_START_DATE']
		lc_start_time = params['LC_START_TIME']
		lc_stop_date = params['LC_STOP_DATE']
		lc_stop_time = params['LC_STOP_TIME']
		lc_live_time = params['LC_LIVE_TIME']
		delta_time = params['DELTA_TIME']
		n_bins = params['N_BINS']
		term_cond = params['TERM_COND']
		lc_url = params['LC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_position(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		grb_ra = params['GRB_RA']
		grb_dec = params['GRB_DEC']
		grb_error = params['GRB_ERROR']
		grb_inten = params['GRB_INTEN']
		grb_signif = params['GRB_SIGNIF']
		img_start_date = params['IMG_START_DATE']
		img_start_time = params['IMG_START_TIME']
		tam = params['TAM[0-3]']
		amplifier = params['AMPLIFIER']
		waveform = params['WAVEFORM']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_spectrum(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		spec_start_date = params['SPEC_START_DATE']
		spec_start_time = params['SPEC_START_TIME']
		spec_stop_date = params['SPEC_STOP_DATE']
		spec_stop_time = params['SPEC_STOP_TIME']
		livetime = params['LIVETIME']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		mode = params['MODE']
		waveform = params['WAVEFORM']
		lrpd_bias = params['LRPD_BIAS']
		term_cond = params['TERM_COND']
		spec_url = params['SPEC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_spectrum_proc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		spec_start_date = params['SPEC_START_DATE']
		spec_start_time = params['SPEC_START_TIME']
		spec_stop_date = params['SPEC_STOP_DATE']
		spec_stop_time = params['SPEC_STOP_TIME']
		livetime = params['LIVETIME']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		mode = params['MODE']
		waveform = params['WAVEFORM']
		lrpd_bias = params['LRPD_BIAS']
		term_cond = params['TERM_COND']
		spec_url = params['SPEC_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_sper(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		sper_start_date = params['SPER_START_DATE']
		sper_start_time = params['SPER_START_TIME']
		sper_stop_date = params['SPER_STOP_DATE']
		sper_stop_time = params['SPER_STOP_TIME']
		tot_expo_time = params['TOT_EXPO_TIME']
		n_pkts = params['N_PKTS']
		n_evts = params['N_EVTS']
		sper_url = params['SPER_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_sper_proc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		sper_start_date = params['SPER_START_DATE']
		sper_start_time = params['SPER_START_TIME']
		sper_stop_date = params['SPER_STOP_DATE']
		sper_stop_time = params['SPER_STOP_TIME']
		tot_expo_time = params['TOT_EXPO_TIME']
		n_pkts = params['N_PKTS']
		n_evts = params['N_EVTS']
		sper_url = params['SPER_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_threshpix(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		tp_start_date = params['TP_START_DATE']
		tp_start_time = params['TP_START_TIME']
		tp_expo_time = params['TP_EXPO_TIME']
		tp_url = params['TP_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']

	@staticmethod
	def process_swift_xrt_threshpix_proc(value):

		params = Processor.decode_classic_notice(value)

		title = params['TITLE']
		notice_date = params['NOTICE_DATE']
		notice_type = params['NOTICE_TYPE']
		trigger_num = params['TRIGGER_NUM']
		point_ra = params['POINT_RA']
		point_dec = params['POINT_DEC']
		tp_start_date = params['TP_START_DATE']
		tp_start_time = params['TP_START_TIME']
		tp_expo_time = params['TP_EXPO_TIME']
		tp_url = params['TP_URL']
		sun_postn = params['SUN_POSTN']
		sun_dist = params['SUN_DIST']
		moon_postn = params['MOON_POSTN']
		moon_dist = params['MOON_DIST']
		moon_illum = params['MOON_ILLUM']
		gal_coords = params['GAL_COORDS']
		ecl_coords = params['ECL_COORDS']
		comments = params['COMMENTS']		

	@staticmethod
	def route_alert(value, topic):

		Utils.log('Alert received from ' + topic + ' channel.', 'info')

		# Kafka notices
		if topic == 'igwn.gwalert':
			Processor.process_igwn_gwalert(value)

		elif topic == 'gcn.notices.icecube.lvk_nu_track_search':
			Processor.process_lvk_nu_track_search(value)

		elif topic == 'gcn.notices.swift.bat.guano':
			Processor.process_swift_bat_guano(value)

		# Fermi notices
		elif topic == 'gcn.classic.text.FERMI_GBM_ALERT':
			Processor.process_fermi_gbm_alert(value)

		elif topic == 'gcn.classic.text.FERMI_GBM_FIN_POS':
			Processor.process_fermi_gbm_fin_pos(value)

		elif topic == 'gcn.classic.text.FERMI_GBM_FLT_POS':
			Processor.process_fermi_gbm_flt_pos(value)

		elif topic == 'gcn.classic.text.FERMI_GBM_GND_POS':
			Processor.process_fermi_gbm_gnd_pos(value)

		elif topic == 'gcn.classic.text.FERMI_GBM_POS_TEST':
			Processor.process_fermi_gbm_pos_test(value)

		elif topic == 'gcn.classic.text.FERMI_GBM_SUBTHRESH':
			Processor.process_fermi_gbm_subthresh(value)

		elif topic == 'gcn.classic.text.FERMI_LAT_MONITOR':
			Processor.process_fermi_lat_monitor(value)

		#elif topic == 'gcn.classic.text.FERMI_LAT_OFFLINE':
		#	Processor.process_fermi_lat_offline(value)

		elif topic == 'gcn.classic.text.FERMI_LAT_POS_TEST':
			Processor.process_fermi_lat_pos_test(value)

		elif topic == 'gcn.classic.text.FERMI_POINTDIR':
			Processor.process_fermi_pointdir(value)

		# IceCube notices
		elif topic == 'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE':
			Processor.process_icecube_astrotrack_bronze(value)

		#elif topic == 'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD':
		#	Processor.process_icecube_astrotrack_gold(value)

		#elif topic == 'gcn.classic.text.ICECUBE_CASCADE':
		#	Processor.process_icecube_cascade(value)

		# Swift notices
		elif topic == 'gcn.classic.text.SWIFT_ACTUAL_POINTDIR':
			Processor.process_swift_actual_pointdir(value)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_LC':
			Processor.process_swift_bat_grb_lc(value)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK':
			Processor.process_swift_bat_grb_pos_ack(value)

		elif topic == 'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST':
			Processor.process_swift_bat_grb_pos_test(value)

		elif topic == 'gcn.classic.text.SWIFT_BAT_QL_POS':
			Processor.process_swift_bat_ql_pos(value)

		elif topic == 'gcn.classic.text.SWIFT_BAT_SCALEDMAP':
			Processor.process_swift_bat_scaledmap(value)

		#elif topic == 'gcn.classic.text.SWIFT_BAT_TRANS':
		#	Processor.process_swift_bat_trans(value)

		elif topic == 'gcn.classic.text.SWIFT_FOM_OBS':
			Processor.process_swift_fom_obs(value)

		elif topic == 'gcn.classic.text.SWIFT_POINTDIR':
			Processor.process_swift_pointdir(value)

		elif topic == 'gcn.classic.text.SWIFT_SC_SLEW':
			Processor.process_swift_sc_slew(value)

		elif topic == 'gcn.classic.text.SWIFT_TOO_FOM':
			Processor.process_swift_too_fom(value)

		elif topic == 'gcn.classic.text.SWIFT_TOO_SC_SLEW':
			Processor.process_swift_too_sc_slew(value)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST':
			Processor.process_swift_uvot_dburst(value)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_DBURST_PROC':
			Processor.process_swift_uvot_dburst_proc(value)

		#elif topic == 'gcn.classic.text.SWIFT_UVOT_EMERGENCY':
		#	Processor.process_swift_uvot_emergency(value)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART':
			Processor.process_swift_uvot_fchart(value)

		elif topic == 'gcn.classic.text.SWIFT_UVOT_FCHART_PROC':
			Processor.process_swift_uvot_fchart_proc(value)

		#elif topic == 'gcn.classic.text.SWIFT_UVOT_POS':
		#	Processor.process_swift_uvot_pos(value)

		#elif topic == 'gcn.classic.text.SWIFT_UVOT_POS_NACK':
		#	Processor.process_swift_uvot_pos_nack(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_CENTROID':
			Processor.process_swift_xrt_centroid(value)

		#elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE':
		#	Processor.process_swift_xrt_image(value)

		#elif topic == 'gcn.classic.text.SWIFT_XRT_IMAGE_PROC':
		#	Processor.process_swift_xrt_image_proc(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_LC':
			Processor.process_swift_xrt_lc(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_POSITION':
			Processor.process_swift_xrt_position(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM':
			Processor.process_swift_xrt_spectrum(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPECTRUM_PROC':
			Processor.process_swift_xrt_spectrum_proc(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER':
			Processor.process_swift_xrt_sper(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_SPER_PROC':
			Processor.process_swift_xrt_sper_proc(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX':
			Processor.process_swift_xrt_threshpix(value)

		elif topic == 'gcn.classic.text.SWIFT_XRT_THRESHPIX_PROC':
			Processor.process_swift_xrt_threshpix_proc(value)

		else:
			print(topic, type(topic))
			params = Processor.decode_classic_notice(value, verbose=True)
			print(params)