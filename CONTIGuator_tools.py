#!/usr/bin/env python

# Copyright 2013 Mitchell Stanton-Cook Licensed under the
# Educational Community License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may
# obtain a copy of the License at
#
# http://www.osedu.org/licenses/ECL-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.


"""
CONTIGuator tools
=================

Set of tools to manipulate CONTIGuator runs. At the moment:
    * Adds unmapped contigs (excluded contigs) to CONTIGuator PseudoContig.embl 
      file

Requires:
    * BioPython
"""

from Bio import SeqIO, SeqFeature
from Bio.Alphabet import generic_dna

import sys, os, traceback, argparse
import time
import ast


def fix_args(args):
    """
    Ensure all arguments with paths are absolute & have simplification removed
    Just apply os.path.abspath & os.path.expanduser

    :param args: the arguments given from argparse
    :returns: an updated args
    """
    args.EMBL  = os.path.expanduser(args.EMBL)
    args.Fasta = os.path.expanduser(args.Fasta) 
    return args


def get_record(embl_path):
    """
    Returns the 1st record (asusmes a pseudo embl = single record file)

    :param embl_path: full path as a string to the EMBL file

    :type embl_path: string

    :returns: a BioPython record
    """
    record = []
    with open(embl_path, "rU") as fin:
        record = list(SeqIO.parse(fin, "embl"))[0]
    return record


def parse_location(location_string):
    """
    Returns the position for which we can begin appending features

    At the moment it is not very elegant and proberly does not need to parse
    the record.features[-1].location string
    
    :param location_string: a string formatted like [xxx:yyy](z)
    """
    start, tmp = location_string[1:].split(':')
    end,   tmp = tmp.split(']')
    strand     = tmp[1:-1]
    # swap start end if we have a -'ve strand
    # TODO check that the following is valid -
    if strand == '-':
        end = start
    return int(end)


def append_excluded_features(fasta_path, original_record, position):
    """
    Append the excluded features given in a fasta file to original record at pos

    :param fasta_path: full path as a string to the fasta file
    :param original_record: the BioPython record
    :param position: the position to start appending features from

    :type fasta_path: string
    :type original_record: BioPython record
    :type position: int
    
    :returns: the updated record
    """
    with open(fasta_path, 'rU') as fasta:
        fa_records = SeqIO.parse(fasta, "fasta")
        for fa_record in fa_records:
            # Append the new sequence to the existing
            original_record.seq = original_record.seq + fa_record.seq
            location = SeqFeature.FeatureLocation(
                    SeqFeature.ExactPosition(position), 
                    SeqFeature.ExactPosition(position+len(fa_record)))
            qual = {'method': ['CONTIGuator/Excluded'], 
                    'systematic_id': [fa_record.id]}
            feature = SeqFeature.SeqFeature(location,type="Contig", 
                                        qualifiers=qual)
            original_record.features.append(feature)
            position = position+len(fa_record)
        return original_record


def write_new_record(updated_record, out=None):
    """
    Write the updated record with the appended features

    :param updated_record:
    :param out: [def = None] full path as a string to write the record to

    :type updated_record: BioPython record
    :type out = Noe or string
    """
    #Needed to set the alphabet to allow writing?
    updated_record.seq.alphabet = generic_dna
    if out == None:
        out = "PseudoContig_Excluded.embl"
    with open(out, "w") as fout:
        SeqIO.write(updated_record, fout, "embl")


def core(args):
    """
    Driver module
    """
    args = fix_args(args)
    print args.EMBL
    rec = get_record(args.EMBL)
    # Get the last features location
    pos = parse_location(str(rec.features[-1].location))
    updated_record = append_excluded_features(args.Fasta, rec, pos)
    write_new_record(updated_record, args.out)


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(description=
                'Add CONTIGuator excluded contigs to PseudoContig.embl file', 
                epilog='By: Mitchell Stanton-Cook, E: m.stantoncook@gmail.com')
        parser.add_argument("EMBL", help="Path to PseudoContig.embl")
        parser.add_argument("Fasta", help="Path to Excluded.fsa")
        parser.add_argument('-o','--out',action='store',
                        default=None, help=('Output PseudoContig_Excluded.embl '
                                            'to this file and location'))
        parser.add_argument('-v','--verbose',action='store_true',
                        default=False, help=('Be loud & noisy'))
 
        parser.set_defaults(func=core)
        args = parser.parse_args()
        if args.verbose:
            print "Executing @ " + time.asctime()
        args.func(args)
        if args.verbose:
            print "Ended @ " + time.asctime()
            print 'Exec time minutes %f:' % ((time.time() - start_time) / 60.0)
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1) 
