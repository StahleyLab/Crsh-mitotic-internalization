import numpy as np
import os
import sys
import argparse
import tifffile

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--npyDir", help="Path to directory containing _seg.npy files")
parser.add_argument("--tifDir", help="Path to directory where you want converted tiffs to be saved")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def main():
    
    if not os.path.exists(args.tifDir):
        os.makedirs(args.tifDir)

    seg_npy_file = "_seg.npy"
    npy_files = os.listdir(args.npyDir)

    for f in npy_files:
        f_path = os.path.join(args.npyDir, f)
        data_to_load = np.load(f_path, allow_pickle=True).item()

        labels = data_to_load['masks']
        filename = data_to_load['filename']
        tifffile.imwrite(os.path.join(args.tifDir, filename), labels)


if __name__ == "__main__":
    main()