import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--old", type=str)
parser.add_argument("--new", type=str)
args = parser.parse_args()

with open('convert_HiggsDNA_to_flashgg_fromMerged_Tprime.py','r') as f:
    program=f.readlines()
program[37] = program[37].replace(args.old,args.new)
program[38] = program[38].replace(args.old,args.new)
program[39] = program[39].replace(args.old,args.new)
#print(program[37])
#print(program[38])
#print(program[39])
with open('convert_HiggsDNA_to_flashgg_fromMerged_Tprime.py', 'w') as f:
    f.writelines(program)



