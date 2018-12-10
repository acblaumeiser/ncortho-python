#Construct a covariance model for a given ncRNA core set

#import RNA
import subprocess

class CmConstructor(object):
    def __init__(self, alignment, outpath):
        self.alignment = alignment
        self.outpath = outpath
    
    def construct(self):
        construct_command = ''
        ### cmbuild -n name of the cm
        subprocess.call(construct_command, shell=True)
        return None
    
    def calibrate(self):
        calibrate_command = ''
        subprocess.call(calibrate_command, shell=True)
        return None

class CovarianceModel(object):
    def __init__(self):
        None
        self.taxa = 'rat, hamster'
    def core_set(self):
        return None
