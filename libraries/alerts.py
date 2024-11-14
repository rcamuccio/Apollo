from config import Configuration
from libraries.utils import Utils
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import requests
import smtplib
import logging
logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

class Alerts:

	@staticmethod
	def alert_team(alert_channel='test', alert_type='test', event_name='test'):
		''' This function sends text messages and emails to the team when an alert is detected.

		:parameter alert_channel - 
		:parameter alert_type - The example alert type for the text to send; the default is to simply send a test alert.
		:parameter event_name - 

		:return - Nothing is returned, however a message is sent to devices.

		'''

		# log in to the email server using the appropriate credentials
		server = smtplib.SMTP_SSL(Configuration.SMTP, Configuration.PORT)
		server.ehlo()
		server.login(Configuration.EMAIL, Configuration.PAS)

		# set up the message based on the type of alert
		msg = MIMEMultipart()
		msg['From'] = Configuration.EMAIL
		msg['To'] = ', '.join(Configuration.MAILING_LIST)

		if alert_channel == 'igwn_gwalert':

			if alert_type == 'M-PRELIMINARY':
				msg['Subject'] = event_name + ': MOCK Preliminary Alert [' + alert_channel + ']\n'
				body = 'A MOCK GW-event has been detected by the ' + alert_channel + ' alert channel. This is a preliminary alert. Check queue and results.\n'

			elif alert_type == 'M-INITIAL':
				msg['Subject'] = event_name + ': MOCK Initial Alert [' + alert_channel + ']\n'
				body = 'A MOCK GW-event has been detected by the ' + alert_channel + ' alert channel. This is an initial alert. Check queue and results.\n'

			elif alert_type == 'M-UPDATE':
				msg['Subject'] = event_name + ': MOCK Update Alert [' + alert_channel + ']\n'
				body = 'A MOCK GW-event has been detected by the ' + alert_channel + ' alert channel. This is an updated alert. Check queue and results.\n'

			elif alert_type == 'M-RETRACTION':
				msg['Subject'] = event_name + ': MOCK Retraction Alert [' + alert_channel + ']\n'
				body = 'This MOCK alert has been retracted by LVK.\n'

			elif alert_type == 'S-PRELIMINARY':
				msg['Subject'] = event_name + ': REAL Preliminary Alert [' + alert_channel + ']\n'
				body = 'A REAL GW-event has been detected by LVK. Check queue and results.\n'

			elif alert_type == 'S-INITIAL':
				msg['Subject'] = event_name + ': REAL Initial Alert [' + alert_channel + ']\n'
				body = 'A REAL GW-event has been detected by LVK. Check queue and results.\n'

			elif alert_type == 'S-UPDATE':
				msg['Subject'] = event_name + ': REAL Update Alert [' + alert_channel + ']\n'
				body = 'A REAL GW-event has been detected by LVK. Check queue and results.\n'

			elif alert_type == 'S-RETRACTION':
				msg['Subject'] = event_name + ': REAL Retraction Alert [' + alert_channel + ']\n'
				body = 'This REAL alert has been retracted by LVK.\n'

			else:
				msg['Subject'] = event_name + ': UNKNOWN (' + alert_type + ') Alert [' + alert_channel + ']\n'
				body = 'An UNKNOWN (' + alert_type + ') alert has been detected by LVK. Check queue and results.\n'

		else:
			Utils.log('Unknown alert channel (' + alert_channel + ').', 'info')

		# convert message to text and send
		msg.attach(MIMEText(body, 'plain'))
		sms = msg.as_string()
		server.sendmail(Configuration.EMAIL, Configuration.MAILING_LIST, sms)

		# clean up
		server.quit()

		# log that the message was sent
		Utils.log('The following alert message was sent to the team: ' + body, 'info')

		return