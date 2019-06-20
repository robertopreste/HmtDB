#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from app import db
from app.site.models import GenAlignment
from app.static.dbdata import last_update
from bioinf.seqs import Alignment
import csv
import os
import logging
import argparse


parser = argparse.ArgumentParser(description="""Modulo per eseguire l'allineamento automatizzato 
    dei nuovi genomi trovati con quelli preesistenti. Se lanciato senza argomenti si presuppone che 
    nello step precedente sia stata usata la data dell'ultimo update del database, altrimenti è 
    possibile specificare la data usata per la ricerca.""")
parser.add_argument("-from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")
parser.add_argument("-mafft_exec", dest="mafft_exec", default="mafft-linux/bin/mafft", type=str,
    help="""Eseguibile Mafft da utilizzare; deve trovarsi nella stessa cartella di questo script. 
    Formato: folder/with/mafft. Default: mafft-linux/bin/mafft""")

args = parser.parse_args()

logging.basicConfig(filename="log_align_seqs_{}.log".format(args.from_date),
                    format="%(levelname)s %(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)

update_dir = "HmtDB_update_{}".format(args.from_date)
fasta_dir = "{}/Fasta_Results_{}".format(update_dir, args.from_date)
gb_dir = "{}/GenBank_Results_{}".format(update_dir, args.from_date)

alg_dir = "multialgs/latest"
new_alg_dir = "multialgs/algs_{}".format(args.from_date)
# expl: total_alg è l'allineamento totale dell'ultimo aggiornamento effettuato
total_alg = "{}/total_alg.fasta".format(alg_dir)
# expl: to_align sono le sequenze nel multifasta appena generato
to_align = "{}/multi_seqs.fasta".format(fasta_dir)

query_alg = "SELECT genomeId, alignment FROM GenAlignment"


def get_alg_from_db():
    """Carica l'allineamento di tutti i genomi contenuti nel db e lo salva in multialgs/latest."""

    os.system("mkdir -p {}".format(alg_dir))

    alg = Alignment()
    res = db.engine.execute(query_alg)

    logging.info("Loading existing alignment from db...")
    for genome_id, sequence in res:
        alg.add_seq(str(genome_id), str(sequence))
    logging.info("Done. Loaded {} sequences in the existing alignment.".format(len(alg)))

    logging.info("Writing alignment to {}...".format(total_alg))
    alg.write_file(total_alg)
    logging.info("Done.")


def upload_alg_on_db():
    """Carica il nuovo allineamento sul db (caricando solo i nuovi genomi allineati)."""

    alg = Alignment()
    logging.info("Loading new alignment...")
    alg.load_file("{}/alg_{}.fasta".format(new_alg_dir, args.from_date))
    logging.info("Done.")

    last_id = int(GenAlignment.query.order_by(GenAlignment.alignmentId.desc()).first().alignmentId)

    logging.info("Loading list of new genomes found...")
    with open("{}/new_hmtdb_ids.txt".format(update_dir), "rb") as list_file:
        results = [el.strip() for el in list_file]
    logging.info("Done. {} aligned genomes will be added to the db.".format(len(results)))

    for n, gen_id in enumerate(results, start=1):
        logging.info("Loading alignment for genome {}...".format(gen_id))
        new_id = last_id + n
        g = GenAlignment(alignmentId=new_id, genomeId=gen_id, alignment=str(alg[gen_id]))
        db.session.add(g)
    logging.info("Uploading new alignments to the db...")
    db.session.commit()
    logging.info("Done.")


def align_seqs():
    """Lancia l'allineamento automatico usando mafft --add. Il template usato è quello con tutte
    le sequenze allineate generate durante l'ultimo aggiornamento, e il subject è il file
    multifasta con le sequenze nuove da allineare, generato nello step precedente."""

    logging.info("Performing automatic multialignment using {} as template and {} as subject sequences...".format(total_alg, to_align))

    os.system("mkdir -p {}".format(new_alg_dir))
    # expl: qui lancio mafft per fare l'allineamento
    os.system("{} --thread 4 --add {} --keeplength {} > {}/alg_{}.fasta".format(args.mafft_exec, to_align, total_alg,
                                                        new_alg_dir, args.from_date))

    logging.info("Alignment complete.")

    return


def main():
    get_alg_from_db()
    align_seqs()
    upload_alg_on_db()


if __name__ == '__main__':
    main()
