import numpy as np
import pandas as pd
import sys
import argparse
from scipy import interpolate

'''
ASSIGN-EDIT.py is meant to parse outputs of Population_code.py and then output as an appropriate tsv file.

Filename for input into Args_and_parent_scraper MUST be in the format of the following:
full/path/to/dir/TYPE_SHAPE_MODEL_PARAM_CORR_[ARBITRARY_SURVEY_LIST].tsv

AND be in a folder with other files that will be used for finding correlations GGN and other arbitrary shape in 
scraper file args. 
'''

# ---SETUP BEGINS HERE---

parser = argparse.ArgumentParser()
parser.add_argument('--BS20', help='True if you want to drop C', default=False)
parser.add_argument('--Hostlib', help='The filepath of the host library')
parser.add_argument('--Hostlib_new', help='The filepath to the directory to write the new survey to')
parser.add_argument('--Args_and_parent_scraper', help="The full file path for parent_scraper method + args grabbing")

args = parser.parse_args()

filename = args.Parent_scraper  # Used to find args + file_path for parent_scraper, which looks for other tsvs
if '/' in filename:
    file_split = filename.split('/')
    file_path = '/'.join(file_split[:-2]) + '/'
    file_name = file_split[-1]
else:
    file_split = filename.split('\\')
    file_path = '\\'.join(file_split[:-2]) + '\\'
    file_name = file_split[-1]

process_args = file_name.split('_')
TYPE = process_args[0]
SHAPE = process_args[1]
MODEL = process_args[2]
PARAM = process_args[3]

# TODO: Properly use CORR to replace where we hardcode look for MASS
tsv_to_hostlib_lookup = {"HOST_LOGMASS": "LOGMASS", "zHD": "ZTRUE", "zCMB": "ZTRUE"}
CORR = process_args[4]
if CORR in tsv_to_hostlib_lookup.keys():
    CORR = tsv_to_hostlib_lookup[CORR]

SURVEY = '_'.join(process_args[5:]).replace('.tsv', '')

if len(SURVEY) == 0:
    print("Malformed path in the file name to parse arguments from. Please check that file conforms to expected "
          "standard.")
    sys.exit(-1)

BS20 = args.BS20
HOSTLIB = args.Hostlib
HOSTLIB_NEW = args.Hostlib_new + SURVEY + '-' + TYPE + '-' + MODEL + '-' + SHAPE + '.HOSTLIB-BS20'

if HOSTLIB is None or HOSTLIB_NEW is None:
    print("No filepath in which to write or no filepath in which to read. Exiting program.")
    sys.exit(-1)

# --- SETUP ENDS HERE ---

# HOSTLIB = '/project2/rkessler/SURVEYS/SDSS/USERS/BAP37/MASS_paper/AS/HOSTLIBS/F1/ASXC.F1.HOSTLIB' # The HOSTLIB that you want to add information to. It will not be overwritten.
# HOSTLIB_NEW = '/project2/rkessler/SURVEYS/SDSS/USERS/BAP37/MASS_paper/AS/HOSTLIBS/F1/'+SURVEY+'-'+TYPE+'-'+MODEL+'.HOSTLIB' #The name and filepath of the new HOSTLIB to be generated

# HOSTLIB = '/project2/rkessler/SURVEYS/SDSS/USERS/BAP37/MASS_paper/AS/HOSTLIBS/F1/BS20-FOUND.F1.HOSTLIB'
# HOSTLIB_NEW = '/project2/rkessler/SURVEYS/SDSS/USERS/BAP37/MASS_paper/AS/HOSTLIBS/F2/' + SURVEY + '-' + TYPE + '-' + MODEL + '-' + SHAPE + '.HOSTLIB-BS20'

# No need to edit anything below here!

def parent_scraper(SURVEY, TYPE, MODEL, PARAM, SHAPE='GGN'):
    scrape_file_location = file_path + TYPE + '_' + SHAPE + '_' + MODEL + "_" + PARAM + "_" + SURVEY + '.tsv'
    #    filename = '/project2/rkessler/SURVEYS/DES/USERS/rebeccachen/photoz_bias/ParentPops/output/'+TYPE+'_'+SHAPE+'_'+MODEL+'_'+PARAM+'_'+SURVEY+'.tsv'
    # dfpops = pd.read_csv(filename, delim_whitespace=True)
    return scrape_file_location


def ggn(x, a, x0, sigma_l, sigma_r, n):
    if x < x0:
        return a * np.exp(-np.abs(x - x0) ** n / (n * sigma_l ** n))
    else:
        return a * np.exp(-np.abs(x - x0) ** n / (n * sigma_r ** n))


def dgauss(x, a1, x0_1, sigma_1, a2, x0_2, sigma_2):
    return a1 * np.exp(-np.abs(x - x0_1) ** 2 / (2 * sigma_1 ** 2)) + a2 * np.exp(
        -np.abs(x - x0_2) ** 2 / (2 * sigma_2 ** 2))


def NAndR(filename):
    with open(filename) as fp:
        for i, line in enumerate(fp):
            if line.startswith('VARNAMES:'):
                line = line.replace(',', ' ')
                line = line.replace('\n', '')
                Names = line.split()
            elif line.startswith('GAL'):
                Startrow = i
                break
    return Names, Startrow

    # Need to load in means/values/etc


# load in as mass, x1mean, x1l, x1r
# and for c, mass, cmean, cl, cr

# x1l += 0.2

# Note for future Brodie add an error trap if one of these doesn't exist!
filename_c = parent_scraper(SURVEY, TYPE, MODEL, 'c', 'GGN')
filename_x = parent_scraper(SURVEY, TYPE, MODEL, 'x1', SHAPE)

dfpopsc = pd.read_csv(filename_c, delim_whitespace=True)
dfpopsx = pd.read_csv(filename_x, delim_whitespace=True)

print(list(dfpopsx))

if SHAPE == 'GGN':
    cmean, cl, cr, mass = dfpopsc.means.values, dfpopsc.stdl.values, dfpopsc.stdr.values, dfpopsc.MASS.values
    # RC: changed n to nfake here temporarily to get it to work
    xmean, xl, xr, mass, xn = dfpopsx.means.values, dfpopsx.stdl.values, dfpopsx.stdr.values, dfpopsx.MASS.values, dfpopsx.nfake.values
elif SHAPE == 'DGaussian':
    cmean, cl, cr, mass = dfpopsc.means.values, dfpopsc.stdl.values, dfpopsc.stdr.values, dfpopsc.MASS.values
    xa1, xmean1, xstd1, xa2, xmean2, xstd2 = dfpopsx.a1.values, dfpopsx.mean1.values, dfpopsx.std1.values, dfpopsx.a2.values, dfpopsx.mean2.values, dfpopsx.std2.values

# The process of writing the events to the HOSTLIB begins below.

bins = np.arange(-.5, .51, .01)
xbins = np.arange(-5, 5.1, .1)

x1dic = {}
cdic = {}

cpdf = {}
x1pdf = {}

namesdic = {}
for i in mass:
    namesdic[str(i)] = i

for i in range(len(cmean)):
    temp = []
    for qc in ((bins[1:] + bins[:-1]) / 2):
        temp.append(ggn(qc, 1, cmean[i], cl[i], cr[i],
                        2))  # Right now, the colour distribution is always an Asymmetric Gaussian.
    temp = np.array(temp)
    ysum = np.cumsum(temp)
    y = temp
    x = (bins[1:] + bins[:-1]) / 2.
    cdf = interpolate.interp1d(ysum / np.amax(ysum), x)
    cdic[str(mass[i]) + '_c'] = cdf
    cpdf[str(mass[i]) + '_c'] = [x, y]

    tempx = []
    for qx in ((xbins[1:] + xbins[:-1]) / 2):
        if SHAPE == 'GGN':  # will need to expand this to account for any shape. This for loop is what handles the shape - once the information is stored here we don't need to worry about shape anymore.
            tempx.append(ggn(qx, 1, xmean[i], xl[i], xr[i], xn[i]))
        elif SHAPE == 'DGaussian':
            tempx.append(dgauss(qx, xa1[i], xmean1[i], xstd1[i], xa2[i], xmean2[i], xstd2[i]))
        else:
            print("Undefined Shape! Quitting.")
            quit()
    tempx = np.array(tempx)
    ysum = np.cumsum(tempx)
    y = tempx
    x = (xbins[1:] + xbins[:-1]) / 2.
    cdf = interpolate.interp1d(ysum / np.amax(ysum), x)
    x1dic[str(mass[i]) + '_x'] = cdf
    x1pdf[str(mass[i]) + '_x'] = [x, y]


def find_nearest(array, value):  # For a given random value, will find the closest thing in the TSV.
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


Names1, StartRow1 = NAndR(HOSTLIB)
df2 = pd.read_csv(HOSTLIB, header=None, skiprows=StartRow1, names=Names1, delim_whitespace=True)

df2['X1'] = -9  # Make sure that any previous values are overwritten, just in case.
df2['C'] = -9

clist = []
for i in range(len(df2)):
    try:
        loc = find_nearest(mass, df2.LOGMASS.values[
            i])  # Will need to make sure that we can read any arbitrary value, not just LOGMASS.
    except TypeError:
        print(
            'Looks like you might not be using the correct TSV file! Make sure that you are pointing to the ITERATIVE == True option')  # bad error trap
        quit()
    try:
        clist.append(cdic[str(loc) + '_c'](np.random.uniform(0, 1, 1))[0])
    except ValueError:
        clist.append(cdic[str(loc) + '_c'](np.random.uniform(0.01, 0.99, 1))[0])

clist = np.array(clist)

x1list = []
for i in range(len(df2)):
    loc = find_nearest(mass, df2.LOGMASS.values[i])
    try:
        x1list.append(x1dic[str(loc) + '_x'](np.random.uniform(0, 1, 1))[0])
    except ValueError:
        print(loc, 'This bin presented an issue in interpolating! Trying again.')
        x1list.append(x1dic[str(loc) + '_x'](np.random.uniform(0.01, 0.99, 1))[0])

x1list = np.array(x1list)

df2['C'] = clist
df2['X1'] = x1list

if BS20:
    df2 = df2.drop(['C'], axis=1)
# df2 = df2.drop(['sSFR', 'SNMAGSHIFT'], axis=1)

print("Done!")
df2.to_csv(HOSTLIB_NEW, sep=' ', header=True, index=False, float_format='%g')
