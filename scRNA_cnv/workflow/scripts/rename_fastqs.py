from pathlib import Path
import re

if __name__ == "__main__":
    #Target format: {BARCODE}_S1_L00{LANE}_R{read_num}_001.fastq.gz
    #original format: {EXPERIMENT}.{BARCODE}.{FLOWCELL}.s_{LANE}.r_{read_num}.fastq.gz
    #find all files matching targets (i.e. correct barcode) and format, and move them to the renamed folder
    orig_dir = Path(snakemake.input["raw_fq_dir"])
    target_dir = Path(snakemake.output[0])
    barcode = snakemake.params["s_id"]
    orig_pattern = r"SLX-[0-9]{1,5}\."+barcode+r"\.[A-Z0-9]+?\.s_(?P<lane>\d)\.(?P<type>[ri])_(?P<read>[12])\.f(ast)?q\.gz"
    print(orig_pattern)
    if not target_dir.is_dir():
        target_dir.mkdir()

    for file in orig_dir.iterdir():
        match = re.match(orig_pattern,file.name)
        if match:
            new_filepath = target_dir / f"{barcode}_S1_L00{match.group("lane")}_{match.group("type").upper()}{match.group("read")}_001.fastq.gz"
            file.rename(new_filepath)
        else:
            print(file)
