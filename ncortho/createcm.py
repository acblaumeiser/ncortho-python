#Construct a covariance model for a given ncRNA

#import RNA
import subprocess

class CmConstructor(object):
    def __init__(self):
        None

class CovarianceModel(object):
    def __init__(self):
        None
        self.taxa = 'rat, hamster'
    def core_set(self):
        return self.taxa

def construct_cm():
    cm = CovarianceModel()
    print cm.core_set()
    
construct_cm()