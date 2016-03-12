#!/usr/bin/env ccp4-python

import matplotlib
matplotlib.use('Cairo')
import matplotlib.pyplot


class Plotter(object):
    """ Class to plot graphs """
    
    def __init__(self):
        matplotlib.rcParams.update({'family': 'Arial', 'font.size': 8})
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
    
    def initialise(self, **kwargs):
        fig = matplotlib.pyplot.figure(**kwargs)
        self.ax = fig.add_subplot(111)
        
    def axisLimits(self, length, offset=0):
        self.ax.set_xlim([offset-1, offset+length])
        self.ax.set_ylim([offset-1, offset+length])
    
    def axisTitles(self, x_axis, y_axis):
        self.ax.set_xlabel(x_axis)
        self.ax.set_ylabel(y_axis)    
    
    def addScatter(self, Xs, Ys, **kwargs):
        self.ax.scatter(Xs, Ys, **kwargs)
        
    def saveFig(self, figurefile):
        matplotlib.pyplot.savefig(figurefile, bbox_inches=0)