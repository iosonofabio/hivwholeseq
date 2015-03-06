# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/03/15
content:    Support module for filenames.
'''
# Modules
from ..filenames import analysis_data_folder



# Globals
purifying_data_folder = analysis_data_folder+'purifying_selection/'



# Functions
def get_fitness_cost_entropy_filename(region, format='pickle'):
    '''Get the filename of the fitness costs for subtype entropy classes'''
    fn = purifying_data_folder+'fitness_cost_entropy_'+region+'.'+format
    return fn

