# For installation
from setuptools import setup, find_packages

# Call setup
setup(name = "tumorNet", 
      description = '''\
        A Cellular Automaton Tumor Growth Simulator
        ''', 
      author = "Alan I. Palma", license = "MIT",
      version = "1.0.0",
      author_email = "alan.palma@yachaytech.edu.ec", 
      packages = find_packages(), 
      install_requires = ["numpy", "matplotlib", "scipy", "pandas",
                          "scienceplots", "argparse",
                          "pytest", "configparser", "imageio"])