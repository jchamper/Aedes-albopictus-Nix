from argparse import ArgumentParser
import subprocess
import pandas as pd
import csv
import re
from pathlib import Path
import time

def parse_slim(slim_string):
    gen = []
    total = []
    sex_rate = []
    wtM = []
    drM = []
    wt_in_M = []
    dr_in_M = []

    lines = slim_string.split('\n')

    for line in lines:
        spaced_line = line.split()
        if line.startswith("OUT:: "):  # Mosquito output
            this_gen = spaced_line[1]
            this_total = spaced_line[2]
            this_sex_rate = spaced_line[3]
            this_wtM = spaced_line[4]
            this_drM = spaced_line[5]
            this_wt_in_M = spaced_line[6]
            this_dr_in_M = spaced_line[7]

            gen.append(this_gen)
            total.append(this_total)
            sex_rate.append(this_sex_rate)
            wtM.append(this_wtM)
            drM.append(this_drM)
            wt_in_M.append(this_wt_in_M)
            dr_in_M.append(this_dr_in_M)

    return gen, total, sex_rate, wtM, drM, wt_in_M, dr_in_M

def run_slim(command_line_args):
    """
    Runs SLiM using subprocess.
    Args:
        command_line_args: list; a list of command line arguments.
    return: The entire SLiM output as a string.
    """
    slim = subprocess.Popen(command_line_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
    out, err = slim.communicate()
    print(out)
# For debugging purposes:
# std.out from the subprocess is in slim.communicate()[0]
# std.error from the subprocess is in slim.communicate()[1]
# Errors from the process can be printed with:
    print(err)
    return out

def configure_slim_command_line(args_dict):
    """
    Sets up a list of command line arguments for running SLiM.
    Args:
        args_dict: a dictionary of arg parser arguments.
    Return
        clargs: A formated list of the arguments.
    """
# We're running SLiM, so the first arg is simple:
    clargs = "slim "
# The filename of the source file must be the last argument:
    source = args_dict.pop("source")
# Add each argument from arg parser to the command line arguments for SLiM:
    for arg in args_dict:
        print(arg)
        if isinstance(args_dict[arg], bool):
            clargs += f"-d {arg}={'T' if args_dict[arg] else 'F'} "
        else:
            clargs += f"-d {arg}={args_dict[arg]} "
        print(clargs)
# Add the source file, and return the string split into a list.
    clargs += source
    return clargs.split()
def main():
    """
    1. Configure using argparse.
    2. Generate the command line list to pass to subprocess through the run_slim() function.
    3. Run SLiM.
    4. Process the output of SLiM to extract the information we want.
    5. Print the results.
    """
# Get args from arg parser:
    parser = ArgumentParser()
    parser.add_argument('-src', '--source', default="Discrete_Homing_Nix_0.75.slim", type=str,
    help=r"SLiM file to be run.")
    parser.add_argument('-outfile', '--output_file', default="Starting_0.75.xlsx", type=str)

    args_dict = vars(parser.parse_args())

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    outfile = args_dict.pop("output_file")
    outfile_stamp = f"{outfile}_{timestamp}.xlsx"

# Next, assemble the command line arguments in the way we want to for SLiM:
    clargs = configure_slim_command_line(args_dict)
# Run the file with the desired arguments.
    slim_result = run_slim(clargs)
# Parse and analyze the result.
    parsed_result = parse_slim(slim_result)
    writeOutput(outfile_stamp, parsed_result)

def write_Output(file, data):
    f = open(file, 'a')
    for value in data:
        f.write(str(value) + ',')
    f.write('\n')
    f.close()

def writeOutput(file, data):

    # 定义表头
    headers = ["Generation", "Total", "Sex Rate", "WT Males", "Red Males", "WT Freq in Males", "DR Freq in Males"]

    # 创建数据框并确保数据为数字类型
    df = pd.DataFrame(data).T
    df = df.apply(pd.to_numeric, errors='coerce')  # 报错时会将无法转换的值变为 NaN

    df.to_excel(file, index=False, header=headers, sheet_name='Sheet1')

if __name__ == "__main__":
    main()
