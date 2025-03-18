import pandas as pd
import sys

# Check if the script receives the file path as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_file>")
    sys.exit(1)

# Get the file path from the command-line argument
input_file = sys.argv[1]

# Read the input file
df = pd.read_csv(input_file, header=None, sep=r'\s+', index_col=0).T.dropna()

Helical = (df.sum() + df.mean()) / 360
Lk = round(Helical[1])
print('Helical:', Helical[1], 'Lk:', Lk)
Writhe = (Lk - Helical)

ind = ['"V{}"'.format(i) for i in range(1, len(Writhe) + 1)]
integers = ['{}'.format(i) for i in range(1, len(Writhe) + 1)]

fileW = '"",' + ','.join(ind) + '\n"1",' + ','.join(integers) + '\n"2",' + ','.join(Writhe.apply("{:.03f}".format))
with open('Circular/writhe.csv', 'w') as f:
    f.write(fileW)

fileT = '"",' + ','.join(ind) + '\n"1",' + ','.join(integers) + '\n"2",' + ','.join(Helical.apply("{:.03f}".format))
with open('Circular/twist.csv', 'w') as f:
    f.write(fileT)

