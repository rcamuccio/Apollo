from config import Configuration
from libraries.utils import Utils

delphi_list = Configuration.DELPHI_DIRECTORIES
Utils.create_directories(delphi_list)

Utils.log('APOLLO setup complete.', 'info')