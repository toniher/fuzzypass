from fuzzypass import *
import argparse
import sys

def main(args):

	parser = inputparser.InputParser( args[1], args[2] )
	parser.read()
	
	fuzzout = fuzzifier.FuzzyFier( {} )
	exit()
	fuzzout.calculate( parser.list_seqs )
	
	print( fuzzout.list_fuzz )

if __name__== "__main__":
	main(sys.argv)



