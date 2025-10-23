"""
This python script generates automatically the INI file that can be used for running the simulation of the moduel TumorNet.
The user could only modify the parameters from the INI file already created or create a new one from scratch following the same structure. 
"""
# Third party libraries
import configparser

# Create an empty dictionary to store the config_orbit.ini file information
config = configparser.ConfigParser()

config['SIMULATION'] = {
    'nx': '80',
    'ny': '80',
    'neighborhood': 'moore', # Options: 'moore' or 'von_neumann'
    'boundary': 'reflective', # Options: 'reflective' or 'periodic'
    'dt': '1.0',
    'steps': '240',
    'ps': '0.5', # Stem cell probability
    'alpha': '0.0',
    'prolif_capacity': '5',
    'mean_cycle': '24.0',
    'sd_cycle': '2.0',
    'diff_constant': '1.0',
    'decay': '0.0',
    'initial_nutrient': '1.0',
    'uptake_nonstem': '0.02',
    'uptake_stem': '0.01',
    'diffusion_substeps': '5',
    'chemotaxis_beta': '3.0',
}

config['THERAPY'] = {
    'therapy_on': 'True',
    'therapy_start': '10.0',
    'therapy_duration': '50.0',
    'therapy_period': '200.0',
    'therapy_kill_prob_nonstem': '0.5',
    'therapy_kill_prob_stem': '0.1',
    'therapy_reduce_prolif': '1.0',
    'therapy_apply_each_step': 'True'
}

config['INIT_SEED'] = {
    'init_seed': 'cluster', # Options: 'single', 'cluster' or 'random'
    'seed_count': '5', # Used only for 'random' init_seed
    'seed_coords': '' # Format: '(x1,y1);(x2,y2);...' Used only for specific coordinates
    }
config['OUTPUT'] = {
    'save_gif': 'True',
    'save_time_series_png': 'True',
    'save_time_series_csv': 'True',
    'output_dir': 'output',
    'output_gif': 'tumor.gif',
    'output_time_series_png': 'time_series.png',
    'output_time_series_csv': 'counts.csv',
    'frames_capture_every': '1',
    'fps': '10'
}

# Save the config file
with open('config.ini', 'w') as configfile:
  config.write(configfile)
print("Configuration file 'config.ini' generated successfully.")