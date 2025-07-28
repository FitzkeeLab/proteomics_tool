#
# title: uniprot_fetch.py
# author: Nick Fitzkee (nfitzkee@chemistry.msstate.edu)
# date:   July 20, 2025
# summary: Fetch proteomics data given a list of UniProtKB IDs
#

import urllib.request, json
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def open_list(fname):
    """
    Open a list of UniProt IDs; supports comments
    
    Args:
        fname (str): A File name of UniProt IDs, # indicates the comments
        
    Returns:
        List of strings: a list of UniProt IDs contained in the file
    
    """
    
    result = []

    with open(fname) as uplist:
        l = uplist.readline()

        while l:
            l = l.strip()
            if not l or l[0] == '#':
                l = uplist.readline()
                continue

            result.append(l)
            l = uplist.readline()

    return result


def protein_info(up):
    """
    Given a UniProt ID, Proteomics-relevant information
    
    Args:
        up (str): A UniProt ID
        
    Returns:
        Dict: a dictionary datastructure with protein information
    """
    
    # First, get the protein sequence
    url = 'https://www.uniprot.org/uniprot/%s.fasta' % up
    with urllib.request.urlopen(url) as page:
        fasta_seq = page.read().decode('utf-8')
        aaseq = ''.join(fasta_seq.split('\n')[1:])
        
    #print(aaseq)
    
    # Second, get the known functions (see below)
    functions = []
    processes = []
    
    # note the split funciton, which removes isoforms for IDs
    
    url = 'https://rest.uniprot.org/uniprotkb/%s?fields=go_f&format=json' % up.split('-')[0]
    with urllib.request.urlopen(url) as page:
        go_terms = page.read().decode('utf-8')
        go_json = json.loads(go_terms)
        go_list = go_json.get('uniProtKBCrossReferences', [])
        
        for go_elem_list in go_list:
            for go_key in go_elem_list['properties']:
                if go_key['key'] == 'GoTerm':
                    go_func_string = (go_key['value'])
                    
                    # C - cellular location
                    # F - molecular function
                    # P - biological process
                    
                    # (using go_f doesn't seem to get just 'F' terms)
                    
                    if go_func_string[0] == 'P':
                        processes.append(go_func_string.split(':')[1].strip())
    
                    if go_func_string[0] == 'F':
                        functions.append(go_func_string.split(':')[1].strip())

    #print(functions)
        
    my_protparam_analysis = ProteinAnalysis(aaseq)

    mw = my_protparam_analysis.molecular_weight()
    isoelec_pi = my_protparam_analysis.isoelectric_point()
    
    result = {'uniprot':       up,
              'seq':           aaseq,
              'go_functions':  functions,
              'go_processes':  processes,
              'pI':            isoelec_pi,
              'molec_mass':    mw
              }
    
    return result


def main(uplist):
    uniprot_list = open_list(uplist)

    prots = {}
    
    proc_file = 'go_processes.txt'
    func_file = 'go_functions.txt'
    
    pf = open(proc_file, 'w')
    ff = open(func_file, 'w')
    
    for upid in uniprot_list:
        info = protein_info(upid)
        prots[upid] = info
    
        
        pf.write('%-10s %10.2f %6.3f ' % (upid, info['molec_mass'], 
                                    info['pI']))
        ff.write('%-10s %10.2f %6.3f ' % (upid, info['molec_mass'], 
                                    info['pI']))
        print('%-10s %10.2f %6.3f ' % (upid, info['molec_mass'], 
                                    info['pI']))
        
        fns = info['go_functions']
        fns.sort()
        
        for fn in fns:
            ff.write('(%s) ' % fn)
        
        pcs = info['go_processes']
        pcs.sort()
        
        for pc in pcs:
            pf.write('(%s) ' % pc)
        
        pf.write('\n')
        ff.write('\n')
        
if __name__ == '__main__':
    import os, sys

    try:
        fn = sys.argv[1]
    except:
        print('usage: %s <uniprot list file>' % os.path.split(sys.argv[0])[1])
        sys.exit(1)

    main(fn)
    
              
