from pathlib import Path
from argparse import ArgumentParser
import subprocess

def main()->None:
    parser = ArgumentParser()
    parser.add_argument("-f","--filestem",type=str,required=True,help="Filestem to look for")
    parser.add_argument("-k","--array_id",type=str,required=False,default="0",help="This jobs array ID")
    parser.add_argument("-p","--padding",type=int,required=True,help="Number of digits to pad numbers to")
    parser.add_argument("-s","--single_barcode",action="store_true",help="Run in single barcode mode")
    parser.add_argument("-n","--num_files",type=int,required=True,help="Number of files to expect, numbered starting at 1")
    args = parser.parse_args()
    
    files = []

    if args.single_barcode:
        padded_num = "0"*(args.padding-len(args.array_id)) + args.array_id
        for i in range(1,args.num_files+1):
            numstr = str(i)
            numstr = "0"*(args.padding-len(numstr))+numstr
            files.append(f"{args.filestem}{numstr}.{padded_num}.tsv")
        outfile = f"{args.filestem}{padded_num}.tsv"
        subprocess.run(f"cat {' '.join(files)} > {outfile}",shell=True)
    else:
        for i in range(1,args.num_files+1):
            numstr = str(i)
            numstr = "0"*(args.padding-len(numstr))+numstr
            files.append(f"{args.filestem}{numstr}.tsv")
        outfile = f"{args.filestem}all.tsv.gz"
        print(files)
        subprocess.run(f"cat {' '.join(files)} | gzip -c -9 - > {outfile}",shell=True)


if __name__=="__main__":
    main()

