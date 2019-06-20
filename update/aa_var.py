#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from optparse import OptionParser
from bioinf.seqs import Alignment
from glob import glob
import os
import logging
import subprocess
from app import db
from app.static.dbdata import last_update

PROGNAME = "Aa Var"
PROGDESC = "Aminoacid Variability Calculation"
PROGVER = "0.1"
PROGLOG = """
AaVar 
Based on the original aa_var/sitevar written by Francesco Rubino.
Edited and optimized by Roberto Preste.

Programma per l'analisi di variabilità aminoacidica.
"""

parser = OptionParser(prog=PROGNAME, description=PROGDESC, version="%%prog %s" % PROGVER)

parser.add_option("-F", "--from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")
parser.add_option("-H", "--healthy", action="store_true", dest="on_healthy",
                  help="Launch the program on healthy or patient genomes (default: False).",
                  default=False)
parser.add_option("-f", "--frame", dest="frame", default=1, type=int, help="Frame (default: 1).")
parser.add_option("-t", "--table", dest="table", default=2, type=int, help="Table (default: 2).")
parser.add_option("-s", "--step", dest="step", default=1, type=int, help="""Step: 1 will 
                    produce the aa_var_temp.sh which needs to be run using qsub, 2 will use its 
                    results. """)

(opts, args) = parser.parse_args()

logging.basicConfig(filename="log_aa_var_{}.log".format(opts.from_date),
                    format="%(levelname)s\t%(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)

query_healthy = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE '%s%%')"

query_patient = "SELECT 'PA_EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE 'PA_%s%%')"

query_compl = "SELECT (SELECT haplotypeHmdb FROM Genome WHERE GenomeId = TBAlign.genomeId) AS haplotype, alignment from GenAlignment as TBAlign WHERE genomeId IN (SELECT genomeId FROM Genome WHERE genomeType = '%s')"

query_position = "SELECT geneName, startPosition, endPosition FROM Locus L WHERE locusType = 'CDS'"


def load_alignment(genome_type, continent=None):
    alg = Alignment()
    if continent:
        if genome_type == "N":
            query = query_healthy % continent
            logging.info("Carico allineamento dei genomi healthy del continente {}...".format(continent))
        else:
            query = query_patient % continent
            logging.info("Carico allineamento dei genomi patient del continente {}...".format(continent))

        res = db.engine.execute(query)
        for haplotype, sequence in res:
            alg.add_seq(str(haplotype), str(sequence))
        logging.info("Allineamento del continente {} caricato.".format(continent))

    else:
        if genome_type == "N":
            query = query_compl % genome_type
            logging.info("Carico allineamento dei genomi healthy...")
        else:
            query = query_compl % genome_type
            logging.info("Carico allineamento dei genomi patient...")

        res = db.engine.execute(query)
        for haplotype, sequence in res:
            alg.add_seq(str(haplotype), str(sequence))
        logging.info("Allineamento caricato.")

    return alg


def find_alg_pos(r_pos, alg_seq):
    n_gap = 0

    for a_pos, nuc in enumerate(alg_seq):
        if nuc == "-":
            n_gap += 1

        if a_pos - n_gap + 1 == r_pos:
            return a_pos

    return None


def cut_alignment(alg):
    query = query_position
    res = db.engine.execute(query)
    # tmp_file = "tmp_%s"
    path = os.path.dirname(os.path.abspath(__file__))

    if not os.path.isdir(os.path.join(path, "tmp_aa_var")):
        os.mkdir(os.path.join(path, "tmp_aa_var"))

    for gene, start, end in res:
        logging.info("Gene {}".format(gene))
        logging.info("Calcolo sul genoma {}...".format(alg[0].name))
        a_start = find_alg_pos(start, str(alg[0]))
        a_end = find_alg_pos(end, str(alg[0]))
        logging.info("Mappa da {} a {} (posizioni dell'allineamento)...".format(a_start + 1,
                                                                                a_end + 1))
        nuc_add = ""

        if gene in ["MT-ND2", "MT-ND3", "MT-ND4", "MT-CO3", "MT-CYB"]:
            logging.info("Aggiunta di AA alla fine delle sequenze...")
            nuc_add = "AA"
        elif gene == "MT-ND1":
            logging.info("Aggiunta di A alla fine delle sequenze...")
            nuc_add = "A"

        gene_alg = Alignment()
        for comp_seq in alg:
            gene_alg.add_seq(str(comp_seq.name), str(comp_seq)[a_start: a_end + 1] + nuc_add)
        # gene_file = tmp_file % gene
        # logging.info("Scrittura del file %s..." % gene_file)
        # gene_alg.write_file(gene_file + ".nuc")
        logging.info("Scrittura del file {}...".format(os.path.join(path, "tmp_aa_var", gene)))
        gene_alg.write_file(os.path.join(path, "tmp_aa_var", gene + ".nuc"))


def write_batch(batch_file, fase):
    path = os.path.dirname(os.path.abspath(__file__))

    with open(batch_file, "w") as f:
        if fase == 0:
            ext_in = "*.nuc"
            ext_out = ".aa"
        elif fase == 1:
            ext_in = "*.aa"
            ext_out = ""
        ext_rev = ".nuc_rev"
        f.write("#!/bin/sh\n")

        for fname in glob(os.path.join(path, "tmp_aa_var", ext_in)):
            base_name = os.path.join(path, "tmp_aa_var", os.path.splitext(os.path.basename(fname))[0])
            base_dir = os.path.join(path, "mitvar")
            emb_dir = os.path.join(path, "emboss/emboss")
            out_file = base_name + ext_out
            f.write("\necho 'File %s'\n" % os.path.basename(out_file))
            if fase == 0:
                # todo: devo capire come vengono presi i nomi di questi geni (se dal database o da qualche altra parte)
                if os.path.basename(fname).startswith("nd6") or os.path.basename(fname).startswith("MT-ND6"):
                    logging.info("Utilizzo di revseq per il gene ND6...")
                    out_file = base_name + ext_rev
                    logging.info("revseq -sequence {} -outseq {}".format(fname, out_file))
                    f.write("{}/revseq -sequence {} -outseq {}\n".format(emb_dir, fname, out_file))
                    fname_tmp = out_file
                    out_file = base_name + ext_out
                    # options = " ".join(["%s %s" % (par, val) for par, val in {"-frame": opts.frame, "-table": opts.table}])
                    options = "-frame {} -table {}".format(str(opts.frame), str(opts.table))
                    logging.info("transeq -sequence {} -outseq {} -trim {}".format(fname_tmp,
                                                                                   out_file,
                                                                                   options))
                    f.write("{}/transeq -sequence {} -outseq {} -trim {}\n".format(emb_dir,
                                                                                   fname_tmp,
                                                                                   out_file,
                                                                                   options))
                else:
                    # options = " ".join(["%s %s" % (par, val) for par, val in {"-frame": opts.frame, "-table": opts.table}])
                    options = "-frame {} -table {}".format(str(opts.frame), str(opts.table))
                    logging.info("transeq -sequence {} -outseq {} -trim {}".format(fname, out_file,
                                                                                   options))
                    f.write("{}/transeq -sequence {} -outseq {} -trim {}\n".format(emb_dir, fname,
                                                                                   out_file,
                                                                                   options))
            elif fase == 1:
                if os.path.basename(fname).startswith("nd2") or os.path.basename(fname).startswith("MT-ND2"):
                    logging.info("Cambio dell'aa per ND2...")
                    a = Alignment()
                    a.load_file(fname)
                    for seq in a:
                        seq[0] = "M"
                    a.write_file(fname)
                f.write("{}/mitvarprot -IN={} -OUT={}\n".format(base_dir, fname, out_file))
    logging.info("Esecuzione del batch: le righe seguenti sono prodotte dal batch.")
    logging.info("Verificare la presenza di errori.")
    # debug: questo subprocess.call fa bloccare tutto perché ogni comando mitvarprot presente in aa_var_temp.sh viene eseguito singolarmente,
    # debug: bisognerebbe lanciarlo con qsub, o meglio ancora lanciare un qsub per ogni comando mitvarprot
    # test: per ora provo con un solo qsub - NON VA DIO CANE
    if fase == 0:
        if opts.step == 1:
            subprocess.call("sh " + batch_file, cwd=".", shell=True)
    # subprocess.call("qsub -q bigmpi2@sauron.recas.ba.infn.it -d /lustrehome/preste/porcoddiodihmtdb/HmtDB_flask/db_update " + batch_file, cwd=".", shell=True)
    # test: rimuovo il subprocess, che invece viene chiamato manualmente usando il cazzo di qsub, e
    # test: inserisco un exit
    if fase == 1:
        if opts.step == 1:
            logging.info("Exiting. Please run aa_var_temp.sh using qsub and then rerun aa_var.py -s 2")
            sys.exit()


def load_from_file():
    path = os.path.dirname(os.path.abspath(__file__))
    var = {}
    ext = ".rates.tdt"
    for fname in glob(os.path.join(path, "tmp_aa_var", "*" + ext)):
        gene = os.path.basename(fname).replace(ext, "")
        logging.info("Caricamento della variabilità del gene {}...".format(gene))

        with open(fname, "U") as f:
            rates = []
            f.readline()
            for line in f:
                value = line.split("\t")[1]
                if value == "nan":
                    rates.append(0.0)
                else:
                    rates.append(float(value))
            var[gene] = rates

    return var


def perform_ops(genome_type, batch_file, continent=None):
    alg = load_alignment(genome_type, continent)

    # alcuni continenti non hanno sequenze (per i pazienti), inserendo 2 anderson si dovrebbe
    # avere il risultato desiderato -> solo 0
    if len(alg) < 2:
        alg.add_seq(str(alg[0].name), str(alg[0]))

    cut_alignment(alg)
    logging.info("Scrittura del batch per l'esecuzione di transeq/revseq...")
    write_batch(batch_file, 0)
    logging.info("Scrittura del batch per l'esecuzione di mitvar_prot...")
    write_batch(batch_file, 1)

    return load_from_file()


def delete_table(genome_type):
    delete_query = "SELECT * FROM AaVariability WHERE genomeType = '%s'" % genome_type
    logging.info("Cancellazione tabella AaVariability, genomi {}...".format(genome_type))
    db.engine.execute(delete_query)
    db.session.commit()

    return


def load_table(var_tot, var_cont, genome_type):
    logging.info("Inizio caricamento dei dati...")
    s_cont = sorted(var_cont)

    base_query = "INSERT INTO AaVariability (aaPos, geneName, rcrsAa, varAa_intrahs, varAa_intermam, %s, genomeType) VALUES (?, ?, ?, ?, ?, %s, ?)"

    sub_base = (", ".join(["varAa_%s" % x.lower() for x in s_cont]),
                ", ".join(["?" for x in s_cont]))
    base_query = base_query % sub_base
    logging.info("Query di inserimento: ")
    logging.info("%s" % base_query)

    path = os.path.dirname(os.path.abspath(__file__))
    righe = []
    for gene in var_tot:
        f = open(os.path.join("mitvar/var_aa_intermam", "aa_var_" + gene + ".txt"))
        seqs = Alignment()
        seqs.load_file(os.path.join(path, "tmp_aa_var", gene + ".aa"))
        rcrs = str(seqs[0])
        f.readline()
        var_mam = [float(line.split("\t")[1]) for line in f]

        for pos in xrange(len(var_tot[gene])):
            dati = []
            dati.append(pos + 1)
            dati.append(gene)
            dati.append(rcrs[pos])
            # dati.append(Decimal("%.4f" % var_tot[gene][pos]))
            dati.append(var_tot[gene][pos])
            # dati.append(Decimal("%.4f" % var_mam[pos]))
            dati.append(var_mam[pos])
            for continent in s_cont:
                # dati.append(Decimal("%.4f" % var_cont[continent][gene][pos]))
                dati.append(var_cont[continent][gene][pos])
            dati.append(genome_type)
            righe.append(dati)
    db.engine.execute(base_query, righe)
    db.session.commit()

    return


def main():
    batch_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "aa_var_temp.sh")

    if opts.on_healthy:
        logging.info("Allineamento dei genomi completi healthy...")
        genome_type = "N"
    else:
        logging.info("Allineamento dei genomi completi patient...")
        genome_type = "P"

    rates_tot = perform_ops(genome_type, batch_file)
    rates_cont = {}

    for cont in ["AF", "AM", "AS", "EU", "OC"]:
        logging.info("Continente {}.".format(cont))
        rates_cont[cont] = perform_ops(genome_type, batch_file, cont)

    delete_table(genome_type)
    load_table(rates_tot, rates_cont, genome_type)


if __name__ == '__main__':
    main()

