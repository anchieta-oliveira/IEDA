""" VOL module for pyepm

This module contains the VOL class, which is used to store and manipulate data from the VOL model. The class has methods
to analyze and visualize data from the VOL model.

Classes
-------
cube
    Class to store and manipulate data from cubes.
dx
    Class to store and manipulate data from dx files.
"""
from ieda.VOL.cube import Cube
from ieda.VOL.dx import DX

__all__ = ['Cube', 'DX']
