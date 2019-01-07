#Construct a covariance model for a given ncRNA core set alignment in Stockholm format

#Python imports
#import RNA
import argparse
import subprocess

#ncOrtho imports
from coreset import CoreSet

class CmConstructor(object):
    
    def __init__(self, alignment, outpath, name, cpu):
        self.alignment = alignment
        self.outpath = outpath
        self.name = name
        self.cpu = cpu
        self.model = '{0}/{1}.cm'.format(outpath, name)
    
    def construct(self):
        cmbuild = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmbuild'
        construct_command = '{3} -n {0} -o {1}/{0}_cmbuild.log {1}/{0}.cm {2}'.format(self.name, self.outpath, self.alignment, cmbuild)
        subprocess.call(construct_command, shell=True)
        #return None
    
    def calibrate(self):
        cmcalibrate = '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/binaries/cmcalibrate'
        calibrate_command = '{0} --cpu {1} {2}'.format(cmcalibrate, self.cpu, self.model)
        subprocess.call(calibrate_command, shell=True)
        #return None

def main():
    parser = argparse.ArgumentParser()
    #args = parser.parse_args()
    #alignment = args.alignment
    #output = args.output
    #cput = args.cpu
    
    #if not alignment:
        #create core set
        #None
        #test_coreset = CoreSet()
    
    test_alignment = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output/rna_aln.sto'
    test_outpath = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/example/output'
    test_name = 'mir-1'
    test_cpu = 4
    test_cc = CmConstructor(test_alignment, test_outpath, test_name, test_cpu)
    #test_cc.construct()
    #test_cc.calibrate()

if __name__ == '__main__':
    main()


'''
	system("$cmbuild -F $covariance_model $stockholm_aln");
#	system("$cmbuild $covariance_model $stockholm_aln");
	print "STEP 04\n";
	system("$cmcalibrate --cpu $cpu $covariance_model");


class CovarianceModel(object):
    def __init__(self):
        None
        self.taxa = 'rat, hamster'
    def core_set(self):
        return None
'''
