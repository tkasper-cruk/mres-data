from argparse import ArgumentParser
from pathlib import Path
import re

def setupParser()->ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument("-f","--fastq_dir",type=str,required=True,help="The path of the directory containing the FastQ files")
    parser.add_argument("-o","--outfile",type=str,required=True,help="The location to write the samplesheet to")
    parser.add_argument("-n","--base_name",type=str,required=True,help="The base name to give to the cells, will be used in the form base_name-BARCODE")
    parser.add_argument("--slx",type=str,required=False,default=r".*?",help="REGEX for the Pool ID, if none is given all fastqs present will be used")
    parser.add_argument("--sample",type=str,required=False,default=r".*?",help="REGEX for the Sample ID, if none is given all fastqs present will be used")
    parser.add_argument("--flowcell",type=str,required=False,default=r".*?",help="REGEX for the Flowcell ID, if none is given all fastqs present will be used")
    parser.add_argument("--lane",type=str,required=False,default=r"\d",help="REGEX for the Lane number, if none is given all fastqs present will be used")
    return parser

def pad_num(num:int,target_len:int)->str:
    numstr = str(num)
    return "0"*(target_len-len(numstr))+numstr

def main()->None:
    parser = setupParser()
    args = parser.parse_args()
    outfile = Path(args.outfile).absolute()
    fqdir = Path(args.fastq_dir)
    if not outfile.parent.is_dir():
        outfile.parent.mkdir(parents=True)
    barcode = r"(?P<barcode>[ACGT]{16})"
    match_statement = rf"(?P<slx>{args.slx})\.(?P<sample>{args.sample})-{barcode}\.{args.flowcell}\.s_{args.lane}\.r_1\.f(ast)?q(\.gz)?"
    counter = 0
    with outfile.open("wt") as outhandle:
        outhandle.write("Pool,Barcode,Sample name\n")
        for file in fqdir.iterdir():
            match = re.match(match_statement,file.name)
            if not match:
                continue
            sample_name = f"{match.group('sample')}_{match.group('slx')}_{pad_num(counter,6)}_{args.base_name}-{match.group('barcode')}"
            outhandle.write(f"{match.group('slx')},{match.group('sample')}-{match.group('barcode')},{sample_name}\n")
            counter += 1

if __name__=="__main__":
    main()