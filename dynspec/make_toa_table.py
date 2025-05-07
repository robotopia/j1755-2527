import argparse
import numpy as np
from astropy.table import QTable

def main():
    parser = argparse.ArgumentParser(description="Construct a LaTeX table from ToA data (default=toas.ecsv) and print to stdout")
    parser.add_argument('--toas_file', default="toas.ecsv", help="Choose a different input file from the default 'toas.ecsv")

    args = parser.parse_args()

    table = QTable.read(args.toas_file)

    print("Telescope & Freq & ToA & Fluence & $\\sigma$ \\\\")
    print(" & (MHz) & (MJD) & ($\\times 10^3\\,$Jy\\,s) & (s) \\\\")
    for row in table:
        toa_err = row['ToA_err']
        toa_precision = int(np.floor(-np.log10(toa_err.to('d').value))) + 1
        toa_str = f"{row['ToA']:.{toa_precision+1}f}({int(np.round(row['ToA_err'].to('d').value*10**(toa_precision+1))):02.0f})"

        fluence_err = row['fluence_err']
        fluence_precision = int(np.floor(-np.log10(fluence_err.to('Jy s').value/1e3))) + 1
        fluence_str = f"{row['fluence'].to('Jy s').value/1e3:.{fluence_precision}f}({int(np.round(fluence_err.to('Jy s').value*10**(fluence_precision-3))):.0f})"
        
        width = row['width']
        width_err = row['width_err']
        width_precision = int(np.floor(-np.log10(width_err.to('s').value))) + 1
        width_str = f"{width.to('s').value:.{width_precision}f}({int(np.round(width_err.to('s').value*10**(width_precision))):.0f})"

        print(f"{row['telescope']} & {row['freq'].to('MHz').value:.0f} & {toa_str} & {fluence_str} & {width_str}")

if __name__ == '__main__':
    main()
