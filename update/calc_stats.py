#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import logging
import argparse
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from app import db
from app.site.models import Genome, Stats
from app.static.dbdata import last_update


parser = argparse.ArgumentParser(description="""Modulo per calcolare le statistiche dei genomi 
    presenti nel database dopo l'aggiornamento. Se lanciato senza argomenti si presuppone che 
    nello step precedente sia stata usata la data dell'ultimo update del database, altrimenti è 
    possibile specificare la data usata per la ricerca.""")
parser.add_argument("-from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")

args = parser.parse_args()

logging.basicConfig(filename="log_calc_stats_{}.log".format(args.from_date),
                    format="%(levelname)s %(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)


def calc_genomes(continent, gen_compl, gen_type):
    """Calculates the total number of genomes belonging to the defined continent, the specified
    type (normal/patient) and completeness.
    :param continent ["AF", "AM", "AS", "EU", "OC"]
    :param gen_compl ["Y", "N"]
    :param gen_type ["N", "P"]
    """

    like_q = continent + "%" if gen_type == "N" else "PA_" + continent + "%"
    q = Genome.query.filter(Genome.completeGenome == gen_compl,
                            Genome.genomeType == gen_type,
                            Genome.haplotypeHmdb.like(like_q)).all()
    logging.info("\tcomplete genomes: {} genome type: {} - found {} genomes".format(gen_compl,
                                                                                    gen_type,
                                                                                    len(q)))

    return len(q)


def update_stats(cont_tup):
    """Updates the values for the provided continent.
    :param cont_tup ("Africa", "AF")
    """

    logging.info("Calculating number of genomes belonging to {}...".format(cont_tup[0]))
    el_y_n = Stats(continentName=cont_tup[0], completeGenome="Y", genomeType="N",
                   total=calc_genomes(cont_tup[1], "Y", "N"))
    el_n_n = Stats(continentName=cont_tup[0], completeGenome="N", genomeType="N",
                   total=calc_genomes(cont_tup[1], "N", "N"))
    el_y_p = Stats(continentName=cont_tup[0], completeGenome="Y", genomeType="P",
                   total=calc_genomes(cont_tup[1], "Y", "P"))
    el_n_p = Stats(continentName=cont_tup[0], completeGenome="N", genomeType="P",
                   total=calc_genomes(cont_tup[1], "N", "P"))
    logging.info("Done.\n")

    db.session.add_all([el_y_n, el_n_n, el_y_p, el_n_p])
    db.session.commit()


def main():
    continents = [("Africa", "AF"),
                  ("America", "AM"),
                  ("Asia", "AS"),
                  ("Europe", "EU"),
                  ("Oceania", "OC"),
                  ("Undefined Continent", "XX")]

    # todo: reset the Stats table (unnecessary, just for this first time, next times use UPDATE)
    logging.info("Resetting the Stats table...")
    db.session.query(Stats).delete()
    db.session.commit()
    logging.info("Done.\n")

    for cont_tup in continents:
        update_stats(cont_tup)


if __name__ == '__main__':
    main()
