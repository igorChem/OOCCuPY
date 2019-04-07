import pandas as pd
import argparse
import os, sys

def mergeTables():
    
    df = pd.read_csv(os.getcwd()+"/"+args.data, delimiter= '\t', header=0)

    df2 = pd.read_csv(os.getcwd()+"/"+args.data2, delimiter= '\t', header=0)

    df[args.option] = pd.Series(df2[args.option])

    df.to_csv(args.output + '.tbl', sep='\t', index=False)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Script developed to merge table files.', formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--originalTable\t\t\t', action='store', required=False, dest='data', help="Specify the name of your file.")
    parser.add_argument('-p2', '--newTable\t\t\t', action='store', required=False, dest='data2', help="Specify the new file.")
    parser.add_argument('-t', '--analysisType\t\t', action='store', required=False, dest='option', type = str, help="Specify the column name.")
    parser.add_argument('-o', '--outputFile\t\t', action='store', required=False, dest='output', type = str, help="Choose the new file name.")
    args = parser.parse_args()
    mergeTables()