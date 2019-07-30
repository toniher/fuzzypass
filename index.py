from fuzzypass import *
import argparse

def main():

	parser = inputparser.InputParser( "Tal", "Tal", "Tal" )
	parser.read()
	dist = metrics.DistPearson().calculate()
	print( "R" )
	print( dist )

if __name__== "__main__":
	main()



