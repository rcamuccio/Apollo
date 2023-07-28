from astropy import units as u
from astropy.io import ascii
from astropy.table import Table, QTable, vstack
from astropy.time import Time
from base64 import b64decode
from gcn_kafka import Consumer
from io import BytesIO
from pprint import pprint
from scipy.stats import norm

import astropy
import astropy_healpix as ah
import configparser
import healpy as hp
import json
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import time

def generate_targets(skymap, detection_time=None, latitude=25.995789, num_remain=50, plot=True):

	print(" Generating targets")

	# --- Read GLADE catalog
	glade_df = pd.read_csv("GLADE+.txt", 
		delimiter=" ", 
		usecols=[2, 8, 9, 32], 
		header=None, 
		names=["GWGC", "ra", "dec", "dist"], 
		low_memory=False)

	# --- Drop NaNs from catalog
	glade_df = glade_df.dropna()

	# --- Slice catalog by declination
	lim_dec = glade_df.dec > (-90 + latitude)
	new_glade_df = glade_df[lim_dec]

	# --- Slice catalog by right ascension
	right_now = detection_time or Time.now()
	right_now = Time(right_now)

	sun = astropy.coordinates.get_sun(right_now)

	horizon = -15 * u.degree
	min_height = 30 * u.degree

	alpha_obs_min = (sun.ra - horizon + min_height).degree
	alpha_obs_max = (sun.ra + horizon - min_height).degree

	circum = 90.0 - abs(new_glade_df.dec) < abs(latitude)
	alfa_min = new_glade_df.ra > float(alpha_obs_min)
	alfa_max = new_glade_df.ra <= float(alpha_obs_max)

	case_1 = (alfa_min & alfa_max) | circum
	case_2 = (alfa_min | alfa_max) | circum

	if alpha_obs_max > alpha_obs_min:
		final_glade_df = new_glade_df[case_1]
	else:
		final_glade_df = new_glade_df[case_2]

	hpx, nside, cred_region_90, ra_max, dec_max, pd_ipix_max, distmu, distsigma, distnorm = read_skymap(skymap, resolution="flat", save=False)

	# --- Calculate likelihoods
	deg2rad = np.pi / 180.0
	thetas = np.pi / 2.0 - final_glade_df["dec"] * deg2rad
	phis = final_glade_df["ra"] * deg2rad
	dists = final_glade_df["dist"]
	ipix = hp.ang2pix(nside, thetas, phis)

	probs = hpx[ipix] * distnorm[ipix] * norm(distmu[ipix], distsigma[ipix]).pdf(dists) * dists ** 2

	# --- Append and sort likelihoods
	final_df = final_glade_df.copy()
	final_df["likelihood"] = list(probs)
	sample_df = final_df.copy()
	sorted_df = sample_df.sort_values(by="likelihood", ascending=False)

	# --- Write top 50 targets to file
	n = 50
	df = sorted_df.head(n)
	df.to_csv("df.csv", index=False)

	if plot:
		graph_targets(df, hpx)
 
	return df

def graph_targets(df, hpx, save=True):

	print(" Graphing targets")

	numfig = 1
	fig = plt.figure(numfig, figsize=(10, 5))
	graph_title = "Test"

	hp.mollview(hpx, title=graph_title, flip="astro", unit="likelihood", fig=numfig, cmap=plt.cm.gist_heat_r)
	fig.axes[1].texts[0].set_fontsize(8)

	ra_targets = df["ra"]
	dec_targets = df["dec"]

	ra_targets_array = ra_targets.to_numpy()
	dec_targets_array = dec_targets.to_numpy()

	hp.projscatter(ra_targets_array, dec_targets_array, lonlat=True, color="green", marker=".")
	hp.graticule()

	if save:
		print(" Saving figure")
		plt.savefig("Targets.png", dpi=1200)
	else:
		plt.show()

	plt.close()

	return

def parse_notice(record):

	# --- Read JSON record
	record = json.loads(record)

	alert_type = record["alert_type"]								# str	
	time_created = record["time_created"]							# str
	superevent_id = record["superevent_id"]							# str
	url = record["urls"]["gracedb"]									# str

	# --- Filter real and mock events
	if superevent_id[0] == "S":
		print(" NOTE: This is a REAL notice!")

	elif superevent_id[0] == "M":
		print(" NOTE: This is a MOCK notice!")

	else:
		print(" WARNING: This is an unrecognized notice!")
		return

	# --- Filter alert type
	if alert_type == "RETRACTION":
		print("", alert_type, "for", superevent_id)
		print(" Event time :", time_created)
		print(" Event URL  :", url)
		print()
		return

	elif alert_type == "EARLYWARNING":
		print("", alert_type, "for", superevent_id)
		print(" Event time :", time_created)
		print(" Event URL  :", url)
		print()

	elif alert_type == "PRELIMINARY":
		print("", alert_type, "for", superevent_id)
		print(" Event time :", time_created)
		print(" Event URL  :", url)
		print()

	elif alert_type == "INITIAL":
		print("", alert_type, "for", superevent_id)
		print(" Event time :", time_created)
		print(" Event URL  :", url)
		print()

	elif alert_type == "UPDATE":
		print("", alert_type, "for", superevent_id)
		print(" Event time :", time_created)
		print(" Event URL  :", url)
		print()

	else:
		print(" WARNING: Unrecognized alert type!")
		return
	
	significant = record["event"]["significant"]					# bool
	time = record["event"]["time"]									# str
	far = record["event"]["far"]									# float
	
	instruments = record["event"]["instruments"]					# list
	group = record["event"]["group"]								# str
	pipeline = record["event"]["pipeline"]							# str
	search = record["event"]["search"]								# str

	has_ns = record["event"]["properties"]["HasNS"]					# float
	has_remnant = record["event"]["properties"]["HasRemnant"]		# float
	has_mass_gap = record["event"]["properties"]["HasMassGap"]		# float

	bns_class = record["event"]["classification"]["BNS"]			# float
	nsbh_class = record["event"]["classification"]["NSBH"]			# float
	bbh_class = record["event"]["classification"]["BBH"]			# float
	terra_class = record["event"]["classification"]["Terrestrial"]	# float

	print("\033[1m" + " [Event]" + "\033[0m" + " --------------------------------------------")
	print()
	print("   Significant    :", significant)
	print("   Time           :", time)
	print("   FAR            :", far)
	print("   Instruments    :", instruments)
	print("   Group          :", group)
	print("   Pipeline       :", pipeline)
	print("   Search         :", search)
	print()
	print("\033[1m" + " [Properties]" + "\033[0m" + " --------------------------------------------")
	print()
	print("   HasNS          :", has_ns)
	print("   HasRemnant     :", has_remnant)
	print("   HasMassGap     :", has_mass_gap)
	print()
	print("\033[1m" + " [Classification]" + "\033[0m" + " --------------------------------------------")
	print()
	print("   BNS            :", bns_class)
	print("   NSBH           :", nsbh_class)
	print("   BBH            :", bbh_class)
	print("   Terrestrial    :", terra_class)
	print()

	# --- Read skymap	
	skymap_str = record.get("event", {}).pop("skymap")

	if skymap_str:

		skymap_bytes = b64decode(skymap_str)
		skymap = Table.read(BytesIO(skymap_bytes))

		hdr_object, hdr_dateobs, hdr_mjdobs, hdr_distmean, hdr_diststd, ra_max, dec_max, area_90 = read_skymap(skymap, resolution="multi", save=False)

		print("\033[1m" + " [Skymap]" + "\033[0m" + " --------------------------------------------")
		print()
		print("   Object                  :", hdr_object)
		print("   Obs Date                :", hdr_dateobs)
		print("   Obs MJD                 :", hdr_mjdobs)
		print("   Mean distance           :", "%.1f" % hdr_distmean, "+/-", "%.1f" % hdr_diststd, "Mpc")
		print("   Most probable location  :", "\u03B1", "=", "%.3f" % ra_max, "deg")
		print("                            ", "\u03B4", "=", "%.3f" % dec_max, "deg")
		print("   Area of 90%             :", area_90)
		print()

	return

def read_skymap(skymap, resolution, save=False):

	if resolution == "flat":

		hpx, distmu, distsigma, distnorm = hp.read_map(skymap, field=range(4))

		npix = len(hpx)
		nside = hp.npix2nside(npix)		# lateral resolution of HEALPix map
		sky_area = 4 * 180**2 / np.pi
		pix_per_deg = sky_area / npix

		# --- Construct map of credible level per pixel
		i = np.flipud(np.argsort(hpx))
		sorted_credible_levels = np.cumsum(hpx[i])
		credible_levels = np.empty_like(sorted_credible_levels)
		credible_levels[i] = sorted_credible_levels

		# --- Calculate area of 90 percent credible region
		cred_region_90 = np.sum(credible_levels <= 0.9) * hp.nside2pixarea(nside, degrees=True)

		# --- Find highest probability pixel
		ipix_max = np.argmax(hpx)

		# --- Calculate highest probability pixel on the sky
		theta_max, phi_max = hp.pix2ang(nside, ipix_max)
		ra_max = np.rad2deg(phi_max)
		dec_max = np.rad2deg(0.5 * np.pi - theta_max)

		# --- Calculate probability density per square degree at that location
		pd_ipix_max = hpx[ipix_max] / hp.nside2pixarea(nside, degrees=True)

		return hpx, nside, cred_region_90, ra_max, dec_max, pd_ipix_max, distmu, distsigma, distnorm 

	elif resolution == "multi":

		# --- Read saved skymap as QTable
		try:
			skymap = QTable.read(skymap)
			i = np.argmax(skymap["PROBDENSITY"])
			skymap[i]["PROBDENSITY"].to_value(u.deg**-2)

		except AttributeError:
			i = np.argmax(skymap["PROBDENSITY"])

		# --- Read header information
		hdr_object = skymap.meta["OBJECT"]		# unique identifier for this event
		hdr_date = skymap.meta["DATE-OBS"]		# UTC of observation
		hdr_mjd = skymap.meta["MJD-OBS"]		# MJD of observation
		hdr_distmean = skymap.meta["DISTMEAN"]	# posterior mean distance [Mpc]
		hdr_diststd = skymap.meta["DISTSTD"]	# posterior std distance [Mpc]

		# --- Most probable sky location
		uniq = skymap[i]["UNIQ"]

		level, ipix = ah.uniq_to_level_ipix(uniq)
		nside = ah.level_to_nside(level)

		ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
		ra_max = ra.deg
		dec_max = dec.deg

		# --- Find 90% probability region
		skymap.sort("PROBDENSITY", reverse=True)

		level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])
		pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

		prob = pixel_area * skymap["PROBDENSITY"]
		cumprob = np.cumsum(prob)
		i = cumprob.searchsorted(0.9)

		area_90 = pixel_area[:i].sum()
		area_90.to_value(u.deg**2)

		# --- Save 90% probability region footprint as MOC
		if save:
			skymap = skymap[:i]
			skymap.sort("UNIQ")
			skymap = skymap["UNIQ",]
			skymap.write("90percent.moc.fits", overwrite=True)

		return hdr_object, hdr_date, hdr_mjd, hdr_distmean, hdr_diststd, ra_max, dec_max, area_90

	else:
		print(" WARNING: Unrecognized skymap!")
		return

def read_notice(event_id, alert_type):

	path = "superevents/" + event_id + "/" + event_id + "-" + alert_type + ".json"

	with open(path, "r") as f:
		record = f.read()

	parse_notice(record)

	return

def retrieve_notice(client_id, client_secret):

	consumer = Consumer(client_id=client_id, client_secret=client_secret)
	consumer.subscribe(["igwn.gwalert"])

	print(" Listening for alerts...")
	while True:
		for message in consumer.consume():
			print("\033[1m"+ " Alert received!")
			parse_notice(message.value())

def test_notice(test_type):

	if test_type == "earlywarning":
		path = "mock_tests/MS181101ab-earlywarning.json"

	elif test_type == "preliminary":
		path = "mock_tests/MS181101ab-preliminary.json"

	elif test_type == "initial":
		path = "mock_tests/MS181101ab-initial.json"

	elif test_type == "update":
		path = "mock_tests/MS181101ab-update.json"

	elif test_type == "retraction":
		path = "mock_tests/MS181101ab-retraction.json"

	with open(path, "r") as f:
		record = f.read()

	parse_notice(record)

if __name__ == "__main__":

	start_time = time.time()

	os.system("clear")
	print("[APOLLO]")
	print()

	# --- Read configuration file
	config = configparser.ConfigParser()
	config.read("config.ini")

	# --- READ NOTICE

	run_read_notice = config["Run"]["read_notice"]

	if run_read_notice == "yes":

		alert_type = config["Directory"]["alert_type"]
		event_id = config["Directory"]["event_id"]

		read_notice(event_id, alert_type)

	# --- RETRIEVE NOTICES ---
	
	run_retrieve_notice = config["Run"]["retrieve_notice"]

	if run_retrieve_notice == "yes":
		
		client_id = config["Client"]["client_id"]
		client_secret = config["Client"]["client_secret"]

		retrieve_notice(client_id, client_secret)

	# --- GENERATE TARGETS ---

	run_generate_targets = config["Run"]["generate_targets"]

	if run_generate_targets == "yes":
		
		latitude = config["Targets"]["latitude"]
		num_remain = config["Targets"]["num_remain"]

		latitude = float(latitude)
		num_remain = int(num_remain)

		event_id = config["Directory"]["event_id"]

		skymap = "superevents/" + event_id + "/bayestar.fits.gz,1"

		generate_targets(skymap, detection_time=None, latitude=latitude, num_remain=num_remain, plot=True)
	
	# --- TEST NOTICES ---

	run_test_notice = config["Run"]["test_notice"]
	if run_test_notice == "yes":
		test_type = config["Test"]["test_type"]
		test_notice(test_type)

	# --- TEST SKYMAPS ---

	run_read_flatres_skymap = config["Run"]["read_flatres_skymap"]
	if run_read_flatres_skymap == "yes":
		flatres_skymap = config["Test"]["flatres_skymap"]
		read_skymap(flatres_skymap, resolution="flat", save=False)

	run_read_multires_skymap = config["Run"]["read_multires_skymap"]
	if run_read_multires_skymap == "yes":
		multires_skymap = config["Test"]["multires_skymap"]
		read_skymap(multires_skymap, resolution="multi", save=False)

	end_time = time.time()
	total_time = end_time - start_time

	print()
	print("[APOLLO ended in", "%.1f" % total_time, "seconds]")