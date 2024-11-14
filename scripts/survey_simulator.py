from config import Configuration
from libraries.priority import Priority
from libraries.utils import Utils
import numpy as np
import pandas as pd

# get the updated field list if necessary
sep = Observatory.get_field_separation(Configuration.OBSERVATORY)
field_list = Priority.field_generator(sep)

# commissioning field locations
# Alessi 13
al_13_ra = 52.0500000
al_13_dec = -35.9000000
al_13_ext = 0.7

# Blanco 1
bl_1_ra = 1.0291667
bl_1_dec = -29.8333333
bl_1_ext = 0.7

# 47 Tucanae
tuc_47_ra = 6.0223292
tuc_47_dec = -72.0814444
tuc_47_ext = 44. / 60.

# generate the field lists for each commissioning field
tuc_47_list, tuc_47_ang = Priority.find_field(field_list, tuc_47_ra, tuc_47_dec, moon_phase=2, ang_extend=tuc_47_ext, program='commissioning')
bl_1_list, bl_1_ang = Priority.find_field(field_list, bl_1_ra, bl_1_dec, moon_phase=2, ang_extend=bl_1_ext, program='commissioning')
al_13_list, al_13_ang = Priority.find_field(field_list, al_13_ra, al_13_dec, moon_phase=2, ang_extend=al_13_ext, program='commissioning')

if (type(tuc_47_list) is str) | (type(bl_1_list) is str) | (type(al_13_list) is str):
	Utils.log('Uh oh...check your field list. It looks like there are no fields to observe.', 'info')
else:
	comms_list = pd.concat([tuc_47_list, bl_1_list, al_13_list]).reset_index(drop=True)
	field_plan = Priority.night_field_selector(field_list, comms_list)
	Utils.log('Field plan has been generated for tonight.', 'info')

Utils.log('See you later, alligator.', 'info')