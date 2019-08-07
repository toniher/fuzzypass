from fuzzypass import *
import argparse
import sys

def main(args):

	parser = inputparser.InputParser( args[1] )
	parser.read()
	
	fuzzout = fuzzifier.FuzzyFier( {} )
	fuzzout.calculate( parser.list_seqs )

	print( fuzzout.list_fuzz )

if __name__== "__main__":
	main(sys.argv)



