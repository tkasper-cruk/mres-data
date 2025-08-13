from argparse import ArgumentParser
from pathlib import Path
import re

def main():
    p = ArgumentParser()
    p.add_argument("-i","--input_dir_stem",type=str,help="Stem for the directories containing the input fastqs")
    p.add_argument("-o","--output_dir",type=str,help="Path to the output directory")
    p.add_argument("-f","--flowcells",type=str,help="Space-separated list of flowcell IDs. The first will be used for the output")
    p.add_argument("-l","--lanes",type=str,help="List of flowcell lanes, order must match flowcells list")
    p.add_argument("--slx",type=str,help="SLX ID")
    p.add_argument("--sample",type=str,help="Sample name")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    if not outdir.is_dir():
        outdir.mkdir(parents=True)
    files = {}
    
    seqruns = list(zip(args.flowcells.split(" "),args.lanes.split(" ")))
    for seqrun in seqruns:
        print(args.slx + r"[.]"+args.sample.strip()+r"-(?P<cell>[ACGT]{16})[.]"+seqrun[0]+r"[.]s_"+seqrun[1]+r"[.]r_(?P<direction>[12])[.]fq")
        file_pattern : re.Pattern = re.compile(
            args.slx + r"[.]"+args.sample.strip()+r"-(?P<cell>[ACGT]{16})[.]"+seqrun[0]+r"[.]s_"+seqrun[1]+r"[.]r_(?P<direction>[12])[.]fq"
        )
        seqdir = Path(f"{args.input_dir_stem}{seqrun[0]}_{seqrun[1]}")
        if not seqdir.is_dir():
            print(f"{seqdir} is not a directory")
        for file in seqdir.iterdir():
            name_match = file_pattern.match(file.name)
            if not name_match:
                print(file.name)
                continue
            cell = name_match.group("cell")
            readfile = name_match.group("direction")
            if cell not in files:
                files[cell] = {"1":{},"2":{}}
            files[cell][readfile][seqrun] = file
    for cell in files:
        if len(files[cell]["1"]) != len(seqruns) or len(files[cell]["2"]) != len(seqruns):
            print(cell, files[cell])
            continue
        for direction in files[cell]:
            outfile = outdir / f"{args.slx}.{args.sample.strip()}-{cell}.{seqruns[0][0]}.s_{seqruns[0][1]}.r_{direction}.fq"
            with outfile.open("wb") as outhandle:
                for file in files[cell][direction].values():
                    with open(file,"rb") as inhadle:
                        outhandle.write(inhadle.read())
                        #outhandle.write(b"\n")
            print(f"Done with R{direction} for {cell}")


if __name__ == "__main__":
    main()
