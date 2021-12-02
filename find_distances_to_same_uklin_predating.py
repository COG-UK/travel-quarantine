#!/usr/bin/env python3
"""

"""

import sys
import argparse
import logging
import csv

__version__ = '0.1'
__date__ = '2020-09-10'
__author__ = ''

# --------------------------------------------------------------------------------------------------

def get_args():
    """
    Parge arguments
    Parameters
    ----------
    no inputs
    Returns
    -------
    oArgs: obj
        arguments object
    """

    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    oParser = argparse.ArgumentParser(description=sDescription)

    oParser.add_argument('--cogmetadata',
                         '-m',
                         required=True,
                         help='COG cog_<DATE>_metadata.csv file')

    oParser.add_argument('--phemetadata',
                         '-p',
                         required=True,
                         help='civet input cvs with additional columns')

    oParser.add_argument('--cogalign',
                         '-a',
                         required=True,
                         help='cog_<DATE>_alignment.fasta file. GOES REALLY BAD IF THIS IS NOT ONE LINE FASTA!')

    oArgs = oParser.parse_args()
    return oArgs

# --------------------------------------------------------------------------------------------------

def main():
    """
    This is the main function, innit?
    """

    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                        level=logging.DEBUG)

    logging.debug("----------- STARTING ----------")

    args = vars(get_args())

    logging.debug("Args: %s", args)

    # read in the phe metadata fromt he civet input file, like region etc.
    phe_metadata = {}
    with open(args['phemetadata'], 'r') as phemd:
        reader = csv.DictReader(phemd)
        for row in reader:
            phe_metadata[row['name']] = {'country': row['region'], 'casecontact': row['casecontact']}

    # read in cogmetadata and store uklin and epiweek per samples
    qry_samples = {}
    with open(args['cogmetadata'], 'r') as cogmetadata:
        reader = csv.DictReader(cogmetadata)
        for row in reader:
            for query_id in phe_metadata:
                if query_id in row['sequence_name']:
                    uklin = row['uk_lineage']
                    epiweek = int(row['epi_week'])
                    qry_samples[query_id] = {'uklin': uklin, 'epiweek': epiweek}

    logging.debug("Read %i query samples", len(qry_samples))

    # output a table
    sys.stdout.write("#query_id\tregion\tcasecontact\tuklin\ttotal\twk_b4\t2wk_b4\t2prec_wks\t4prec_wks\tmin_d\t2follow\t4follow\tPreceding Dist\n")

    # open cog metadata again
    with open(args['cogmetadata'], 'r') as mdfile:
        reader = csv.DictReader(mdfile)

        # iterate through samples
        for qryid in qry_samples:
            # reset values that are outputted for each query sequence.
            cnt = 0
            wkb4_1 = 0
            wkb4_2 = 0
            prec_wks2 = 0
            prec_wks4 = 0
            foll_wks2 = 0
            foll_wks4 = 0

            # get the seq of the query samples from the alignment
            qry_seq = get_seq_from_aln([qryid], args['cogalign'])
            assert len(qry_seq) == 1 # better check ...

            this_sample_epi_week = qry_samples[qryid]['epiweek']

            # iterate through the COG metadata
            beforeseqnames = []
            foll2wks_seqnames = []
            foll4wks_seqnames = []
            for row in reader:
                if row['uk_lineage'] == qry_samples[qryid]['uklin']:
                    # count samples w/ same lineage
                    cnt += 1
                    row_epi_week = int(row['epi_week'])
                    # count all older samples by epiweek
                    if row_epi_week <= this_sample_epi_week - 1:
                        wkb4_1 += 1
                    # count all samples that are two weeks older
                    if row_epi_week <= this_sample_epi_week - 2:
                        wkb4_2 += 1
                    # count all samples in the preceeding 2 epi weeks
                    if row_epi_week < this_sample_epi_week and row_epi_week >= this_sample_epi_week - 2:
                        prec_wks2 += 1
                    # count all samples in the preceeding 4 epi weeks
                    if row_epi_week < this_sample_epi_week and row_epi_week >= this_sample_epi_week - 4:
                        prec_wks4 += 1
                        # store other seqname
                        beforeseqnames.append(row['sequence_name'])
                    # count all samples in the following 2 epi weeks
                    if row_epi_week > this_sample_epi_week and row_epi_week <= this_sample_epi_week + 2:
                        foll2wks_seqnames.append(row['sequence_name'])
                    # count all samples in the following 4 epi weeks
                    if row_epi_week > this_sample_epi_week and row_epi_week <= this_sample_epi_week + 4:
                        foll4wks_seqnames.append(row['sequence_name'])

            # calculate minimum distance to preceding samples
            mindis = -1
            minoth = 1000
            prec4 = ""
            if len(beforeseqnames) > 0:
                # get stored min dist from for this query sample
                otherseqs = get_seq_from_aln(beforeseqnames, args['cogalign'])
                y = list(qry_seq.values())[0]
                mindis = min([calc_pw_dist(y, x) for x in otherseqs.values()])
                #look at all seq in the preceeding 4 weeks
                if len(otherseqs) <= 1000:
                    for ky1, vl1 in otherseqs.items():
                        minoth = 10000
                        for ky2, vl2 in otherseqs.items():
                            if ky2 != ky1:
                                odist = calc_pw_dist(vl1, vl2)
                                if odist < minoth:
                                    minoth = odist
                        if minoth < 10000:
                            prec4 += "{},".format(minoth)
                else:
                    prec4 = "-1"



            # calculate how many samples are within 2 SNPs in the follwing 2 and 4 weeks
            foll_wks2 = 0
            foll_wks4 = 0
            if len(foll4wks_seqnames) > 0:
                # logging.debug("%i in the following 4 weeks", len(foll4wks_seqnames))
                # logging.debug("%i in the following 2 weeks", len(foll2wks_seqnames))
                # get all seqs in the 4 weeks after from alignment - 2 weeks are a strict subset
                otherseqs = get_seq_from_aln(foll4wks_seqnames, args['cogalign'])
                # y is now the query sequence
                y = list(qry_seq.values())[0]
                all_dists = {}
                wks2_dists = []
                for (nme, sq1) in otherseqs.items():
                    nme = nme.replace(">", "").strip()
                    dst = calc_pw_dist(y, sq1)
                    # all_dists now contains all seqs from the 4 weeks after
                    all_dists[nme] = dst
                    if nme in foll2wks_seqnames:
                        # wks2_dists contains a subset of those
                        wks2_dists.append(dst)
                # logging.debug("4 weeks all dists: %s", all_dists.values())
                # logging.debug("2 weeks all dists: %s", wks2_dists)
                foll_wks4 = sum([1 if d <= 2 else 0 for d in all_dists.values()])
                foll_wks2 = sum([1 if d <= 2 else 0 for d in wks2_dists])

            # output table row for this query sample
            # logging.debug("%i samples with UK lineage %s, %i earlier ones, same as %s", cnt, qry_samples[qryid]['uklin'], earlier, qryid)
            sys.stdout.write("%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\n" % (qryid,
                                                                                   phe_metadata[qryid]['country'],
                                                                                   phe_metadata[qryid]['casecontact'],
                                                                                   qry_samples[qryid]['uklin'],
                                                                                   cnt,
                                                                                   wkb4_1,
                                                                                   wkb4_2,
                                                                                   prec_wks2,
                                                                                   prec_wks4,
                                                                                   mindis,
                                                                                   foll_wks2,
                                                                                   foll_wks4,
                                                                                   prec4))

            # reset the iterator of the COG metadatafile for next query sample
            mdfile.seek(0)
            next(reader)

    logging.debug("----------- FINISHED ----------")
    return 0

# end of main --------------------------------------------------------------------------------------

def calc_pw_dist(s1, s2):
    """
    Calculate the pw distance between 2 seqs

    Parameters
    ------
    s1: str
        sequence 1
    s2: str
        seq2

    Returns
    -------
    d: int
        pw SNP dist bnetween them
    """

    assert len(s1) == len(s2)

    s1 = s1.upper()
    s2 = s2.upper()
    d = 0
    for a, b in zip(s1, s2):
        if not a in ['A', 'C', 'G', 'T']:
            continue
        if not b in ['A', 'C', 'G', 'T']:
            continue
        if a != b:
            d += 1

    return d

# ---------------------------------------------------------------------------------------------------

def get_seq_from_aln(qryids, cogalign):
    """
    Get the seqquence for this ID from the alignment.

    GOES TERRIBLY WRONG IF THIS ISN'T ONE LINE FASTA as in
    >head1
    AACGTAGCTAC
    >head2
    ACTGTACGTAC
    ...

    Parameters
    ------
    qryids: list
        list of COGID to be mached to seq header in alignment
    cogalign: str
        alignemnt fasta file name

    Returns
    -------
    answer: dict
        {header: seq, head2: seq, ...}
        empty if none found
    """

    answer = {}

    with open(cogalign, 'r') as f:
        while True:
            header = f.readline()
            seq = f.readline()
            if not seq:
                break  # EOF

            for qid in qryids:
                if qid in header:
                    answer[header.strip()] = seq.strip()

    return answer

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
