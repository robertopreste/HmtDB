#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
import os
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import logging
from datetime import datetime
# import argparse
from optparse import OptionParser
from bioinf.seqs import Alignment
from bioinf.sitevar import SiteVar, Variability, out_delimited, out_classic, gen_relative_pos
from bioinf.tab_tor import TabTorCln
from bioinf.distances import dists
from bioinf.comp import NucComposition
from app import db
from app.static.dbdata import last_update


PROGNAME = "Nt Var"
PROGDESC = "Nucleotide Variability Calculation"
PROGVER = "0.1"
PROGLOG = """
NtVar 
Based on the original nt_var/sitevar written by Francesco Rubino.
Edited and optimized by Roberto Preste.

Programma per l'analisi di variabilità nucleotidica.
"""

parser = OptionParser(prog=PROGNAME, description=PROGDESC, version="%%prog %s" % PROGVER)

parser.add_option("-F", "--from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")
# # reference sequence
# input/output files
# parser.add_option("-i", "--input-file", action="store", type="string",
#     dest="infile", help="alignment file to be processed (default: %default)", default="infile")
# parser.add_option("-o", "--output-file", action="store", type="string", dest="outfile",
#                   help="Alignment file to be processed (default: infile.out)")
# reference sequence related
parser.add_option("-H", "--healthy", action="store_true", dest="on_healthy",
                  help="Launch the program on healthy or patient genomes (default: False).",
                  default=False)
parser.add_option("-C", "--continent", action="store", type=str, dest="on_continent",
                  help="Perform continent-specific calculation (default: no).", default=None)
parser.add_option("-M", "--complete", action="store_true", dest="on_complete",
                  help="Perform calculations only on complete genomes (default: False).",
                  default=False)
parser.add_option("-r", "--reference-sequence", action="store", type="int", dest="refseq",
                  help="Reference sequence position in alignment file (default: %default).",
                  default=1)
parser.add_option("-x", "--reference-include", action="store_true", dest="refinc",
                  help="Include anderson/reference sequence in the process (default: %default).",
                  default=False)
parser.add_option("-R", "--reference-positions", action="store_false", dest="refpos",
                  help="Don't use anderson/reference sequence relative positions in output file (default: False).",
                  default=True)
#dloop - regions related
parser.add_option("-l", "--dloop", action="store_true", dest="dloop",
                  help="Dloop and coding: must set beginning and end of coding region (default: %default).",
                  default=False)
parser.add_option("-s", "--coding-start", action="store", type="int", dest="coding_beg",
                  help="Start (first) position for coding region.", default=0)
parser.add_option("-e", "--coding-end", action="store", type="int", dest="coding_end",
                  help="End (last) position for coding region.", default=0)
# complete - incomplete sequences
# if --dloop option is used you must specify the number of consecutive gaps
# which will be used to assume if a sequence is complete or not
parser.add_option("-n", "--gaps-number", action="store", type="int", dest="gaps_number",
                  help="Start position for coding region (default: %default).", default=50)
# customize Weights for gap, transition, transversion, 2/distMean add
parser.add_option("-g", "--weight-gap", action="store", type="float", dest="gap",
                  help="Value used for gap penalty (default: %default).", default=0.0)
parser.add_option("-p", "--weight-transition", action="store", type="float", dest="transition",
                  help="Value used for transition penalty (default: %default).", default=1.0)
parser.add_option("-q", "--weight-transversion", action="store", type="float", dest="transversion",
                  help="Value used for transversion penalty (default: %default).", default=2.0)
parser.add_option("-a", "--weight-ambiguity", action="store", type="float", dest="ambiguity",
                  help="Value used for ambiguity penalty (default: %default).", default=0.1)
parser.add_option("-m", "--weight-site-gap", action="store_false", dest="wgaps",
                  help="Disable addition of 2/distMean for a site containing a gap (default: enabled).",
                  default=True)
# distance type
parser.add_option("-d", "--distance-formula", action="store", type="string", dest="distance",
                  help="Distance formula to be used: %s (default: %%default)." % ', '.join(dists.keys()),
                  default='kimura')
# verbose
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                  help="Be verbose (default: %default).", default=False)
parser.add_option("-V", "--more-verbose", action="store_true", dest="more_verbose",
                  help="Be more verbose (default: %default).", default=False)
# output
parser.add_option("-c", "--classic-output", action="store_true", dest="classic",
                  help="Uses classic output format (default: %default).", default=False)
parser.add_option("-X", "--exponential", action="store_true", dest="exp",
                  help="In classic format uses exponential numbers (default: %default).", default=False)
parser.add_option("-D", "--delimiter", action="store", type=str, dest="delimiter",
                  help="Output file delimiter (default: ,).", default=",")
#tab and tor files
parser.add_option("-t", "--tab-file", action="store_true", dest="tab",
                  help="Export Tab file (default: %default).", default=False)
# tab options
parser.add_option("-f", "--tab-file-small", action="store_true", dest="tabsmall",
                  help="Export Tab file with max 255 columns (for import in Spredsheet Program) (default: %default).",
                  default=False)
parser.add_option("-w", "--tab-min-var", action="store_false", dest="tab_max",
                  help="Use first 255 less variant sites to produce Tab File (default: False).",
                  default=True)
parser.add_option("-k", "--tab-sort-var", action="store_false", dest="tab_sort_pos",
                  help="Sort sites by variability instead of positions in Tab File (default: by position).",
                  default=True)
parser.add_option("-T", "--tor-file", action="store_true", dest="tor",
                  help="Export Tor File (default: %default).", default=False)
parser.add_option("-L", "--no-cleanup", action="store_false", dest="cleanup",
                  help="Don't use cleanup (default: False).", default=True)

(opts, args) = parser.parse_args()

logging.basicConfig(filename="log_calc_sitevar_{}.log".format(opts.from_date),
                    format="%(levelname)s\t%(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="a", level=logging.DEBUG)


query_healthy = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE '%s%%')"

query_healthy_only_complete = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE '%s%%' AND completeGenome = 'Y')"

query_patient = "SELECT 'PA_EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE 'PA_%s%%')"

query_patient_only_complete = "SELECT 'PA_EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE 'PA_%s%%' AND completeGenome = 'Y')"

query_compl = "SELECT (SELECT haplotypeHmdb FROM Genome WHERE GenomeId = TBAlign.genomeId) AS haplotype, alignment from GenAlignment as TBAlign WHERE genomeId IN (SELECT genomeId FROM Genome WHERE genomeType = '%s')"

query_only_complete = "SELECT (SELECT haplotypeHmdb FROM Genome WHERE GenomeId = TBAlign.genomeId) AS haplotype, alignment from GenAlignment as TBAlign WHERE genomeId IN (SELECT genomeId FROM Genome WHERE genomeType = '%s' AND completeGenome = 'Y')"


def check_opts(opts):
    if opts.dloop:
        if opts.coding_beg < 1 or opts.coding_end <= opts.coding_beg:
            logging.critical("Invalid values for start or end of coding region.")
            sys.exit(1)
    # todo: qui ci vanno le opts per verbose e more_verbose


def check_progress(msg):
    if len(msg) == 3:
        par = (msg[0] + 1, msg[1], "all" if msg[2] else "complete")
        logging.info("Analysing region (%i-%i) on %s sequences..." % par)

    if len(msg) == 5:
        par = (msg[0] + 1, msg[1], msg[2] + 1, msg[3], "all" if msg[-1] else "complete")
        logging.info("Analysing region (%i-%i,%i-%i) on %s sequences..." % par)


def main():
    started = datetime.now()
    check_opts(opts)
    alg = Alignment()
    sitev = SiteVar()
    refseq = None

    os.system("mkdir -p sitevars")

    # Load alignment

    if opts.on_continent:
        if opts.on_complete:
            if opts.on_healthy:
                query = query_healthy_only_complete % opts.on_continent
                logging.info("Carico allineamento dei genomi healthy del continente {}...".format(opts.on_continent))
                outname = "healthy_complete_{}".format(opts.on_continent)
            else:
                query = query_patient_only_complete % opts.on_continent
                logging.info("Carico allineamento dei genomi patient del continente {}...".format(opts.on_continent))
                outname = "patient_complete_{}".format(opts.on_continent)
        else:
            if opts.on_healthy:
                query = query_healthy % opts.on_continent
                logging.info("Carico allineamento dei genomi healthy del continente {}...".format(opts.on_continent))
                outname = "healthy_{}".format(opts.on_continent)
            else:
                query = query_patient % opts.on_continent
                logging.info("Carico allineamento dei genomi patient del continente {}...".format(opts.on_continent))
                outname = "patient_{}".format(opts.on_continent)

        res = db.engine.execute(query)
        for haplotype, sequence in res:
            alg.add_seq(str(haplotype), str(sequence))

        if len(alg) == 1:
            logging.warning("Allineamento del continente {} saltato.".format(opts.on_continent))
            return None, None
        elif len(alg) == 0:
            logging.critical("Haplotype Hmdb mancanti.")
            sys.exit(1)
        else:
            logging.info("Allineamento del continente {} caricato.".format(opts.on_continent))

    else:  # allineamento tutti
        if opts.on_complete:
            if opts.on_healthy:
                query = query_only_complete % "N"
                logging.info("Carico allineamento dei genomi healthy...")
                outname = "healthy_complete"
            else:
                query = query_only_complete % "P"
                logging.info("Carico allineamento dei genomi patient...")
                outname = "patient_complete"
        else:
            if opts.on_healthy:
                query = query_compl % "N"
                logging.info("Carico allineamento dei genomi healthy...")
                outname = "healthy"
            else:
                query = query_compl % "P"
                logging.info("Carico allineamento dei genomi patient...")
                outname = "patient"

        res = db.engine.execute(query)
        for haplotype, sequence in res:
            alg.add_seq(str(haplotype), str(sequence))
        logging.info("Allineamento caricato.")

    logging.debug("Loaded {} sequences, with lengths {}.".format(len(alg), alg.seq_len))

    if opts.dloop:
        regions = [(0, opts.coding_beg - 1, opts.coding_end, alg.seq_len, False),
                   (opts.coding_beg - 1, opts.coding_end, True)]
    else:
        regions = [(0, alg.seq_len, True)]

    # Reference sequence analysis

    if opts.refpos:
        start_block = datetime.now()
        logging.info("Calculating relative positions...")
        aseq = str(alg[opts.refseq - 1])
        try:
            refseq = gen_relative_pos(aseq.replace("-", ""), aseq)
        except IndexError:
            logging.error("An error occurred, reverting to normal positions.")
        if not opts.refinc:
            sitev.exclude = opts.refseq - 1
            logging.info("Removed reference sequence from analysis.")
        end_block = datetime.now()
        logging.debug("Reference positions execution time: {}".format(str(end_block - start_block)))
        # todo: forse qui sopra si potrebbe togliere str()

    # Variability processing

    sitev.set_weights(opts.gap, opts.transition, opts.transversion, opts.ambiguity)
    sitev.alg = alg
    sitev.regions = regions
    sitev.gapnumber = opts.gaps_number

    if opts.wgaps:
        sitev.wgaps = True
    else:
        logging.info("Disabling addition of 2/distMean for a site containing a gap.")
        sitev.wgaps = False

    dist_mod = None
    logging.debug("Using {} distance.".format(opts.distance))
    try:
        from bioinf import c_distances
        dist_mod = c_distances
    except ImportError:
        logging.error("Reverting to Python distance functions.")
        from bioinf import distances
        dist_mod = distances
    finally:
        dist_func = getattr(dist_mod, dists[opts.distance] + "R", None)
        sitev.set_dist_type(dist_func)

    try:
        from bioinf import c_variability
        sitev.set_var_type(c_variability.Variability)
    except ImportError:
        logging.error("Reverting to Python variability class.")
        sitev.set_var_type(Variability)

    start_block = datetime.now()
    sitev.run_analysis()
    end_block = datetime.now()
    logging.debug("Variability analysis execution time: {}.".format(str(end_block - start_block)))
    # todo: forse qui sopra si potrebbe togliere str()

    # Composition and output file generation

    logging.info("Calculating nucleotide composition...")
    start_block = datetime.now()
    comp = NucComposition()
    comp.alg = alg
    if opts.dloop:
        comp.set_ranges((0, opts.coding_beg - 1, opts.coding_end),
                        (opts.coding_beg - 1, opts.coding_end, alg.seq_len),
                        (list(sitev.comp_seqs), list(sitev.comp_seqs.union(sitev.part_seqs)),
                         list(sitev.comp_seqs)))

    comp.calc_comp()
    end_block = datetime.now()
    logging.debug("Nucleotide composition calculation time: {}.".format(str(end_block - start_block)))
    # todo: forse qui sopra si potrebbe togliere str()

    outfile = "sitevars/sitevar_{}".format(outname)
    logging.info("Saving output to file: {}...".format(outfile))
    if opts.classic:
        out_classic(comp, sitev.rates, outfile, refseq, exp=opts.exp)
    else:
        out_delimited(comp, sitev.rates, outfile, refseq, delimiter=opts.delimiter, exp=opts.exp)

    # Tab and Tor file generation

    if opts.tab:
        start_block = datetime.now()
        logging.info("Generating tab file...")
        t = TabTorCln()
        t.alg = alg
        t.cln_create(opts.cleanup)
        t.rates = sitev.rates
        t.ref_seq = opts.refseq - 1
        t.tab_create()
        end_block = datetime.now()
        logging.info("Saving tab file to {}.tab...".format(outfile))
        t.tab_write(outfile + ".tab", 0)
        logging.debug("Tab file creation execution time: %s." % str(end_block - start_block))

        logging.info("Saving duplicates file to {}.red...".format(outfile))
        with open(outfile + ".red", "w") as f:
            for lst in alg.get_duplicates():
                f.write(alg[lst[0]].name + ": " + ", ".join([alg[idx].name for idx in lst[1:]]) + "\n")

        if opts.tabsmall:
            logging.info("Saving small tab file to {}.tab255...".format(outfile))
            t.in_rates = sitev.rates
            t.tab_write(outfile + ".tab255", 255, opts.tab_sort_pos, opts.tab_max)

        if opts.tor:
            logging.info("Saving tor file to {}.tor...".format(outfile))
            t.tor_create(outfile + ".tor")

    ended = datetime.now()
    logging.info("Total execution time: {}.".format(str(ended - started)))
    # todo: forse qui sopra si potrebbe togliere str()

# todo: magari splittare questo unico run() in diverse funzioni


if __name__ == '__main__':
    main()

