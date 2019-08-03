#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import argparse
from bioinf.seqs import Alignment
from zipfile import ZipFile
import os
import logging
import csv
from app import db
from app.static.dbdata import last_update

parser = argparse.ArgumentParser(description="""Modulo per generare i file zip contenenti la 
    variabilità e gli allineamenti totali e continente-specifici dei genomi sul database. Se 
    lanciato senza argomenti si presuppone che nello step precedente sia stata usata la data 
    dell'ultimo update del database, altrimenti è possibile specificare la data usata per la 
    ricerca.""")
parser.add_argument("-from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")

args = parser.parse_args()

logging.basicConfig(filename="log_zip_data_{}.log".format(args.from_date),
                    format="%(levelname)s\t%(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)

query_alg = "SELECT (SELECT haplotypeHmdb FROM Genome G WHERE G.genomeId = A.genomeId) as HAPLOTYPE_HMDB, Alignment FROM GenAlignment A WHERE A.genomeId IN (SELECT individualId FROM IndividualsData WHERE countryId IN (SELECT countryId FROM Country %s)) %s ORDER BY HAPLOTYPE_HMDB"

query_cont = " WHERE continentCode = '%s' "
query_gentype = " AND (SELECT genomeId FROM Genome WHERE genomeId = A.genomeId AND genomeType = '%s') IS NOT NULL "

query_var = "SELECT nucleotidePosition, insertionPosition, var_%s, compVar_%s FROM NtVariability WHERE genomeType = '%s' ORDER BY nucleotidePosition, insertionPosition"


def extract_alg(genome_type, continent=None):
    if continent:
        query = query_alg % (query_cont % continent, query_gentype % genome_type)
    else:
        query = query_alg % ("", query_gentype % genome_type)

    res = db.engine.execute(query)
    alg = Alignment()

    for name, seq in res:
        alg.add_seq(str(name), str(seq))
    logging.info("Numero di sequenze nell'allineamento: {}.".format(len(alg)))

    if genome_type == "P":
        base_name = "alg_pa_"
    else:
        base_name = "alg_"
    if continent:
        base_name += continent
    else:
        base_name += "tot"

    base_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), base_name)

    alg.write_file(base_name + ".fasta")
    logging.info("Salvataggio nel file {}.zip...".format(base_name))
    alg_zip = ZipFile(base_name + ".zip", "w")
    alg_zip.write(base_name + ".fasta", os.path.basename(base_name) + ".fasta")
    os.remove(base_name + ".fasta")
    alg_zip.close()
    os.system("mv {}.zip zips/".format(base_name))


def extract_var(genome_type, continent=None):
    if continent:
        query = query_var % (continent.lower(), continent.lower(), genome_type)
    else:
        query = query_var % ("tot", "tot", genome_type)

    res = db.engine.execute(query)

    if genome_type == "P":
        base_name = "var_pa_"
    else:
        base_name = "var_"
    if continent:
        base_name += continent
    else:
        base_name += "tot"

    base_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), base_name)

    with open(base_name + ".csv", "w") as f:
        f.write("Position,Variability,\%A,\%C,\%T,\%G,\%GAP\n")

        for pos, ins, var, comp in res:
            pos = ".".join((str(pos), str(ins)))
            try:
                # expl: perché nel db il campo compVar è composto da float separati da ;
                comp = comp.split(";")
                f.write(",".join([pos, str(var)] + comp) + "\n")
            except AttributeError, e:
                logging.warning("Problema nel parsing dei dati: {}; {}.".format(comp, e))
                logging.warning("Nessun dato per questa sezione.")
                break

    logging.info("Salvataggio nel file {}.zip...".format(base_name))
    alg_zip = ZipFile(base_name + ".zip", "w")
    alg_zip.write(base_name + ".csv", os.path.basename(base_name) + ".csv")
    os.remove(base_name + ".csv")
    alg_zip.close()
    os.system("mv {}.zip zips/".format(base_name))


def main():
    os.system("mkdir -p zips/")
    for genome_type in ["N", "P"]:
        for cont in ["AF", "AM", "AS", "EU", "OC"]:
            logging.info("Estrazione dell'allineamento dei genomi {} per il continente {}...".format(genome_type, cont))
            extract_alg(genome_type, cont)
            logging.info("Estrazione della variabilità dei genomi {} per il continente {}...".format(genome_type, cont))
            extract_var(genome_type, cont)
        logging.info("Estrazione dell'allineamento completo dei genomi {}...".format(genome_type))
        extract_alg(genome_type)
        logging.info("Estrazione della variabilità totale dei genomi {}...".format(genome_type))
        extract_var(genome_type)

    # os.system("mv zips/ ../app/static/")

# todo: creazione dei file direttamente nella cartella zips dell'applicazione, per ora lo fa in una cartella zips locale


if __name__ == '__main__':
    main()


