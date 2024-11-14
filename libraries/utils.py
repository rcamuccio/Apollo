from config import Configuration
import logging
import os

class Utils:

	@staticmethod
	def create_directories(directory_list):

		for path in directory_list:
			if os.path.exists(path) is False:
				os.mkdir(path)

	@staticmethod
	def log(statement, level):

		# create logger
		logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', filename=Configuration.LOGS_DIRECTORY + 'main.log', filemode='a')
		logger = logging.getLogger()

		if not getattr(logger, 'handler_set', None):
			logger.setLevel(logging.DEBUG)

			# create console handler and set level to debug
			ch = logging.StreamHandler()
			ch.setLevel(logging.DEBUG)

			# create formatter
			formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')

			# add formatter to ch
			ch.setFormatter(formatter)

			# add ch to logger
			logger.addHandler(ch)

			# 'set' handler
			logger.handler_set = True

		if level == 'info':
			logger.info(statement)

		if level == 'debug':
			logger.debug(statement)

		if level == 'warning':
			logger.warning(statement)

		if level == 'error':
			logger.error(statement)

		if level == 'critical':
			logger.critical(statement)