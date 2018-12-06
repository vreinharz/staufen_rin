import argparse
from collections import defaultdict
from itertools import product
import os
import pickle
from pprint import pprint
import re


from tqdm import tqdm

PATH_IL = os.path.join('..', 'Data', 'clusters_IL.pickle')
PATH_STEMS = os.path.join('..', 'Data', 'hStau1_SRS_loops/posplus/loops_in_SRS_flanking2nt.txt')


def read_stauffen(file_path):
    """Reads the file and creates a generator of triplets with 
        i) the "name" (line that starts with the > character
       ii) A tuple for the sequence
      iii) A tuple for the secondary structure
    """
    with open(file_path) as f:
        #This is a neat trick to read 3 lines at a time
        #If the number of lines is not a multiple of 3 it will
        #NOT yield the last ones
        for name, seq, ss in zip(f, f, f):
            name = name.strip()
            seq = tuple(seq.split())
            ss = tuple(ss.split())
            yield name, seq, ss


def read_restrain_noncanonical(cluster_path):
    """Dumb filter function to read cluster and keep only conserved clusters that have at least 
    one non canonical interaction"""
    with open(cluster_path, 'rb') as f:
        data = pickle.load(f)
    non_canonical = [c if any(e[2]['label'] != 'CWW' and e[2]['label'] != 'B53' for e in c['graph'].edges(data=True))
                     else None for c in data]
    return non_canonical


def graph2sequence(g, exact_len=False, bps=0):
    """Makes a regex from a graph, only needs to match the non-canonical positions.
    if exact_len is True then we enforce the length of each loop too.
    
    Since every interior loop we have is surrounded by 2 base pairs, we don't count them and just 
    put a "." in the regex. The argument "bps" tells how many of them need an exact match
    """

    #We just grab all the nodes
    l = {x:g.node[x]['nt'] for x in sorted(g.nodes())}
    #Look at the canonical base pairse positions
    cww = [(e[0], e[1]) for e in g.edges(data=True) if e[2]['label'] == 'CWW' and e[0] < e[1]]
    #By design first position must be limit of the 5' side of the IL
    left_min = min(x[0] for x in cww)
    left_max = max(x[0] for x in cww)
    #resp. second gives limit on the 3' side of the IL
    right_min = min(x[1] for x in cww)
    right_max = max(x[1] for x in cww)

    #make bad pattern, () are just to help retrieve positions later
    left = ''.join(f'({l[i]})' if i in l else '.' for i in range(left_min, left_max + 1))
    right = ''.join(f'({l[i]})' if i in l else '.' for i in range(right_min, right_max + 1))
    

    #We just replace all the "dots" for positions without non-canonical interactions with .*
    #inside the expression
    if not exact_len:
        left = '.*'.join([x for x in left.split('.') if x])
        right = '.*'.join([x for x in right.split('.') if x])


    #Since only two base pairs, if bps >=2 we just keep everything
    if bps < 2:
        #3 positions per base pair since there is two () also
        nb_pos_to_remove = 6 - 3 * bps
        nb_gaps = '.' * (2 - bps) 
        
        left = nb_gaps + left[nb_pos_to_remove:-nb_pos_to_remove] + nb_gaps
        right = nb_gaps + right[nb_pos_to_remove:-nb_pos_to_remove] + nb_gaps

    return left, right


def match_graph_stauffen(cluster, seq, bps, exact_len, verbose=False):
    """For a cluster, will check if there is a march from
    any of its graphs in the seq. Since every graph and seq are IL
    we separate them into left/right part. We must have a match in both
    to return the results.

    If verbose is False, it just return the first match, 
    if it is True  it returns a list with all of them
    """

    def _match_sequence(match):
        """Small helper function to create the match sequence"""
        seq = ['.'] * match.endpos
        for i in range(1, match.lastindex + 1):
            seq[match.span(i)[0]] = match.group(i)
        return ''.join(seq)


    matches =[]
    #If there is nothing in the cluster (mostly for no non-canonical)
    #we just return an empty list
    if cluster is None:
        return matches

    #we simply go through all graphs in that cluster
    #Get the sequences with graph2sequence
    #try to match if we do we get the sequence that matched and add it to "matches"
    for g in cluster['l_graphs']:
        left, right = graph2sequence(g, exact_len=exact_len, bps=bps)
        l_match = re.match(left, seq[0])
        r_match = re.match(right, seq[1])
        #This might seems weird but we want to ensure that (a) we have a match
        #And (b) that we still have nucleotides in the match (if we remove the 
        #canonical base pairs the only non-canonical interactions migth be between the 
        #opening / starting bps)
        if l_match and r_match and l_match.groups() and r_match.groups():
            l_match, r_match = _match_sequence(l_match), _match_sequence(r_match)
            matches.append((l_match, r_match))
            if verbose is False:
                return matches
    return matches 


def number_match_per_cluster(stems_path, cluster_path, exact_len, bps, output_path=None):
    """Simply checks for each cluster the number of sequence it can match.
    Prints to STDOUT if output_path is None, else write to it
    """
    clusters = read_restrain_noncanonical(cluster_path)

    output = [['Cluster', 'Matches']]
    for i, c in enumerate(clusters):
        output.append([i, 0])
        for name, seq, ss in read_stauffen(stems_path):
            if match_graph_stauffen(c, seq, exact_len, bps):
                output[-1][1] += 1
    output = '\n'.join('\t'.join(map(str, x)) for x in output)
    
    if output_path is None:
        print(output)
    else:
        with open(output_path, 'w') as f:
            f.write(output)


def all_match_per_cluster(stems_path, cluster_path, exact_len, bps, output_path=None):
    """
    Prints to STDOUT if output_path is None, else write to it
    """
    clusters = read_restrain_noncanonical(cluster_path)


    results = []
    for i, c in enumerate(clusters):
        for name, seq, ss in read_stauffen(stems_path):
            #we are verbose, we want all the things
            out = match_graph_stauffen(c, seq, exact_len, bps, True)
            if out:
                header = '\t'.join((str(i), name, '\t'.join(seq)))
                results.extend(['\t'.join((header, '\t'.join(match)))
                                for match in out])
    output = '\t'.join(("Cluster", "Name", "Seq5'", "Seq3'", "Match5'", "Match3'"))
    output += '\n' + '\n'.join(results)

    if output_path is None:
        print(output)
    else:
        with open(output_path, 'w') as f:
            f.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Identify double stranded RNA sequences that can\
        accomodate conserved non-canonical interaction patterns')

    parser.add_argument('-i', '--input', required=True,
                        help='Path to file with stauffen motifs')
    parser.add_argument('-c', '--cluster',required=True,
                        help='Path to file with graphs clusters')
    parser.add_argument('-o', '--output',default=None,
                         help='File to write output, if not provided\
                         it is written to stdout')
    parser.add_argument('-e', '--exact_len', default=False, type=bool, 
                        choices=(True, False),
                        help="When true the match must also have the exact\
                        same length")
    parser.add_argument('-b', '--bps', type=int, default=0, choices=(0, 1, 2),
                         help='How many base pairs surrounding the IL do we\
                         want a sequence match for. We have at most 2')

    parser.add_argument('-v', '--verbose', type=bool, default=False, 
                        choices=(True, False), 
                        help='If verbose we print all the matches for each sequence for each cluster')


    if not args.verbose:
        number_match_per_cluster(args.input, args.cluster, args.exact_len, args.bps, args.output)
    else:
        all_match_per_cluster(args.input, args.cluster, args.exact_len, args.bps, args.output)
