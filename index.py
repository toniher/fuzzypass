from fuzzypass import *
import argparse

def main():

	parser = inputparser.InputParser( "Tal", "Tal", "Tal" )
	parser.read()
	pearson = metrics.DistPearson()
	pearson.calculate()
	print( pearson.distance )
	print( pearson.rcoeff )

if __name__== "__main__":
	main()



