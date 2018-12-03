import argparse
parser = argparse.ArgumentParser()
#parser.add_argument('hw', help='echo the string you use here')
#parser.add_argument('square', help='display a square of a given number', type=int)
parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
args = parser.parse_args()
#print args.square**2
if args.verbose:
    print 'you want more?'