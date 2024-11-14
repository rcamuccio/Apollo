from config import Configuration
from libraries.processor import Processor
from libraries.utils import Utils
from gcn_kafka import Consumer
import logging
logging.getLogger('gcn').setLevel(logging.WARNING)

consumer = Consumer(client_id=Configuration.CLIENT_ID,
					client_secret=Configuration.CLIENT_SECRET)

available_topics = ['igwn.gwalert',										# finish
					'gcn.notices.icecube.lvk_nu_track_search',			# finish
					'gcn.notices.swift.bat.guano',						# finish
					'gcn.classic.text.FERMI_GBM_ALERT',					#
					'gcn.classic.text.FERMI_GBM_FIN_POS',				#
					'gcn.classic.text.FERMI_GBM_FLT_POS',				#
					'gcn.classic.text.FERMI_GBM_GND_POS',				#
					'gcn.classic.text.FERMI_GBM_POS_TEST',				#
					'gcn.classic.text.FERMI_GBM_SUBTHRESH',				#
					'gcn.classic.text.FERMI_LAT_MONITOR',				#
					'gcn.classic.text.FERMI_LAT_OFFLINE',				
					'gcn.classic.text.FERMI_LAT_POS_TEST',				#
					'gcn.classic.text.FERMI_POINTDIR',					#
					'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',		#
					'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
					'gcn.classic.text.ICECUBE_CASCADE',
					'gcn.classic.text.SWIFT_ACTUAL_POINTDIR',			#
					'gcn.classic.text.SWIFT_BAT_GRB_LC',				#
					'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK',			#
					'gcn.classic.text.SWIFT_BAT_GRB_POS_TEST',			#
					'gcn.classic.text.SWIFT_BAT_QL_POS',				#
					'gcn.classic.text.SWIFT_BAT_SCALEDMAP',				#
					'gcn.classic.text.SWIFT_BAT_TRANS',
					'gcn.classic.text.SWIFT_FOM_OBS',					#
					'gcn.classic.text.SWIFT_POINTDIR',					#
					'gcn.classic.text.SWIFT_SC_SLEW',					#
					'gcn.classic.text.SWIFT_TOO_FOM',					#
					'gcn.classic.text.SWIFT_TOO_SC_SLEW',				#
					'gcn.classic.text.SWIFT_UVOT_DBURST',				#
					'gcn.classic.text.SWIFT_UVOT_DBURST_PROC',			#
					'gcn.classic.text.SWIFT_UVOT_EMERGENCY',
					'gcn.classic.text.SWIFT_UVOT_FCHART',				#
					'gcn.classic.text.SWIFT_UVOT_FCHART_PROC',			#
					'gcn.classic.text.SWIFT_UVOT_POS',
					'gcn.classic.text.SWIFT_UVOT_POS_NACK',
					'gcn.classic.text.SWIFT_XRT_CENTROID',				#
					'gcn.classic.text.SWIFT_XRT_IMAGE',
					'gcn.classic.text.SWIFT_XRT_IMAGE_PROC',
					'gcn.classic.text.SWIFT_XRT_LC',					#
					'gcn.classic.text.SWIFT_XRT_POSITION',				#
					'gcn.classic.text.SWIFT_XRT_SPECTRUM',				#
					'gcn.classic.text.SWIFT_XRT_SPECTRUM_PROC',			#
					'gcn.classic.text.SWIFT_XRT_SPER',					#
					'gcn.classic.text.SWIFT_XRT_SPER_PROC',				#
					'gcn.classic.text.SWIFT_XRT_THRESHPIX',				#
					'gcn.classic.text.SWIFT_XRT_THRESHPIX_PROC']		#

consumer.subscribe(available_topics)

while True:

	for message in consumer.consume(timeout=1):

		message_error = message.error()
		message_topic = message.topic()
		message_value = message.value()

		if message_error:
			Utils.log(message_error, 'info')
			Utils.log(message_topic, 'info')
			Utils.log(message_value, 'info')

		else:
			Processor.route_alert(message_value, message_topic)


unavailable_topics = ['gcn.classic.text.FERMI_GBM_LC',
					'gcn.classic.text.FERMI_GBM_TRANS',
					'gcn.classic.text.FERMI_LAT_GND',
					'gcn.classic.text.FERMI_LAT_POS_DIAG',
					'gcn.classic.text.FERMI_LAT_POS_INI',
					'gcn.classic.text.FERMI_LAT_POS_UPD',
					'gcn.classic.text.FERMI_LAT_TRANS',	
					'gcn.classic.text.FERMI_SC_SLEW',
					'gcn.classic.text.SWIFT_BAT_ALARM_LONG',
					'gcn.classic.text.SWIFT_BAT_ALARM_SHORT',
					'gcn.classic.text.SWIFT_BAT_GRB_ALERT',
					'gcn.classic.text.SWIFT_BAT_GRB_LC_PROC',
					'gcn.classic.text.SWIFT_BAT_GRB_POS_NACK',
					'gcn.classic.text.SWIFT_BAT_KNOWN_SRC',
					'gcn.classic.text.SWIFT_BAT_MONITOR',
					'gcn.classic.text.SWIFT_BAT_SLEW_POS',
					'gcn.classic.text.SWIFT_BAT_SUB_THRESHOLD',
					'gcn.classic.text.SWIFT_BAT_SUBSUB',
					'gcn.classic.text.SWIFT_FOM_PPT_ARG_ERR',
					'gcn.classic.text.SWIFT_FOM_SAFE_POINT',
					'gcn.classic.text.SWIFT_FOM_SLEW_ABORT',
					'gcn.classic.text.SWIFT_XRT_EMERGENCY']