from aod import AOD1b
import sys
import pyshtools as sht

if __name__ == "__main__":
    print("Runstring:  monthly_avg_AOD1B ../AOD1B_2022-01_06.tgz")
    
    a = AOD1b(None)
    a.fname = sys.argv[1]
    a.read_tar()
    date = a.epochs[0]
    year = date.year
    month = date.month
    for t in ['atm','ocn','glo','oba']:
        print("Averaging ",t," for file",a.fname)
        ave = a.average(t)#print(a.type)
        sht.shio.shwrite(f'AOD1B_{year}-{month:02d}_{t}.sph.gz', ave, header=f'AOD1B monthly {year:04d}-{month:02d} {t}', header2='Stokes coefficients')
