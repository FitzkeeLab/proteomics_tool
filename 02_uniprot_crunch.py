#
# title: 02_uniprot_crunch.py
# author: Nick Fitzkee (nfitzkee@chemistry.msstate.edu)
# date:   July 27, 2025
# summary: Generate basic proteomics output
#

import os.path, math

# The first bin contains all entries less than the stated value. The
# last bin contains all entries greater than the state value, and it
# should be repeated from the i-1 value. The general rule for internal
# bins is that each bin contains values greater than equal to the i-1
# bin, but less than the i bin.

# So the MW value below contains a bin for < 20000, a bin for >= 20000
# but less than 40000, etc. all the way up to a final bin for >= 100000.

MW_BINS = [20000, 40000, 60000, 80000, 100000, 100000]
PI_BINS = [4, 5, 6, 7, 8, 9, 10, 10]

# Transform to apply to abundances; this converts log2 abundances
# into raw abundances
TRANSFORM = lambda x: math.pow(2.0, x)

# This is a null transform; it leaves abundances as-is
#TRANSFORM = lambda x: x

# Value to scale final values by. Normally they add up to 1.0, but this
# can be used to scale them to be 100.
SCALE = 100.0

# String to use for "other" functions, to be placed at end
OTHER_STR = 'other'

def open_uplist(fname):
    """
    Open a list of UniProt IDs and associated data; supports comments
    
    Args:
        fname (str): A File name of UniProt IDs, # indicates the comments
        
    Returns:
        Dict of UniProt IDs: a list of MW, pI, and Function
    
    """
    
    result = {}

    with open(fname) as uplist:
        l = uplist.readline()

        while l:
            l = l.strip()
            if not l or l[0] == '#':
                l = uplist.readline()
                continue

            l = l.split('#')[0]
            my_data_str, prot_func = l.split('(')
            prot_func = prot_func[:-1] # remove final ")"
            uniprot_id, mw, pI = my_data_str.split()
            mw = float(mw)
            pI = float(pI)
                        
            result[uniprot_id] = [mw, pI, prot_func]
            
            l = uplist.readline()

    return result

def open_ablist(fname):
    """
    Open a paired list of UniProt IDs and MS Abundances

    Args:
        fname (str): A file name containing UniProt IDs and some numerical amount

    Returns
    -------
        Dict of UniProt IDs mapped to amount
        
    """
    
    result = {}
    total  = 0.0

    with open(fname) as uplist:
        l = uplist.readline()

        while l:
            l = l.strip()
            if not l or l[0] == '#':
                l = uplist.readline()
                continue

            l = l.split('#')[0]
            uniprot_id, amt = l.split()
            
            try:
                amt = float(amt)
                amt = TRANSFORM(amt) # See note about TRANSFORM above
                total = total + amt
                
                result[uniprot_id] = amt
            
            except ValueError:
                pass
            
            
            l = uplist.readline()

    return result, total
    

def main(uplist, ablist):
    prot_data        = open_uplist(uplist)
    prot_amt, total  = open_ablist(ablist)
    
    #print(prot_amt)
    #print(prot_data)
    
    base_fn = os.path.splitext(ablist)[0]
    pI_hist_file = '%s-pi_hist.out.txt' % base_fn
    mw_hist_file = '%s-mw_hist.out.txt' % base_fn
    proc_hist_file = '%s-proc_hist.out.txt' % base_fn
    
    pI_hist = []
    mw_hist = []
    proc_dict = {}
    
    for i in PI_BINS:
        pI_hist.append(0.0)
    
    for i in MW_BINS:
        mw_hist.append(0.0)
        
    uniprot_ids = list(prot_amt.keys())
    uniprot_ids.sort()
    
    for id in uniprot_ids:
        mw, pI, prot_func = prot_data[id]
        amt = prot_amt[id]
        
        # Add to MW Histogram
        if mw < MW_BINS[0]:
            mw_hist[0] = mw_hist[0] + amt/total
        elif mw >= MW_BINS[-1]:
            mw_hist[-1] = mw_hist[-1] + amt/total
        else:
            for i in range(len(MW_BINS)-2):
                bin_min = MW_BINS[i]
                bin_max = MW_BINS[i+1]
                
                if mw >= bin_min and mw < bin_max:
                    mw_hist[i+1] = mw_hist[i+1]+amt/total
        
        # Add to pI Histogram            
        if pI < PI_BINS[0]:
            pI_hist[0] = pI_hist[0] + amt/total
        elif pI >= PI_BINS[-1]:
            pI_hist[-1] = pI_hist[-1] + amt/total
        else:
            for i in range(len(PI_BINS)-2):
                bin_min = PI_BINS[i]
                bin_max = PI_BINS[i+1]
                
                if pI >= bin_min and pI < bin_max:
                    pI_hist[i+1] = pI_hist[i+1]+amt/total

        # Construct process/function list
        
        if prot_func not in proc_dict:
            proc_dict[prot_func] = amt/total
        else:
            proc_dict[prot_func] = proc_dict[prot_func] + amt/total
        

    pif = open(pI_hist_file, 'w')
    mwf = open(mw_hist_file, 'w')
    prf = open(proc_hist_file, 'w')

    #print(mw_hist)
    #print(pi_hist)
    
    # Write the final pI Histogram
    
    pif.write('< %i, %8.3f\n' % (int(PI_BINS[0]), pI_hist[0]*SCALE))
    
    for i in range(len(PI_BINS)-2):
        pif.write('%i-%i, %8.3f\n' % (int(PI_BINS[i]), int(PI_BINS[i+1]), pI_hist[i+1]*SCALE))
    
    pif.write('> %i, %8.3f\n' % (int(PI_BINS[-1]), pI_hist[-1]*SCALE))
    
    # Write the final MW Histogram
    
    mwf.write('< %i, %8.3f\n' % (int(MW_BINS[0]/1000), mw_hist[0]*SCALE))
    
    for i in range(len(MW_BINS)-2):
        mwf.write('%i-%i, %8.3f\n' % (int(MW_BINS[i]/1000), int(MW_BINS[i+1]/1000), mw_hist[i+1]*SCALE))
    
    mwf.write('> %i, %8.3f\n' % (int(MW_BINS[-1]/1000), mw_hist[-1]*SCALE))
    
    # Write the final process file
    
    prs = list(proc_dict.keys())
    prs.sort()
    
    for pr in prs:
        if pr != OTHER_STR:
            prf.write('%s, %8.3f\n' % (pr, proc_dict[pr]*SCALE))
    
    if OTHER_STR in proc_dict:
        prf.write('%s, %8.3f\n' % (OTHER_STR, proc_dict[OTHER_STR]*SCALE))
    
    pif.close()
    mwf.close()
    prf.close()
        
if __name__ == '__main__':
    import os, sys

    try:
        my_uniprot   = sys.argv[1]
        my_abundance = sys.argv[2]
    except:
        print('usage: %s <culled uniprot data file> <uniprot with abundancies>' % os.path.split(sys.argv[0])[1])
        sys.exit(1)

    main(my_uniprot, my_abundance)
    
              
