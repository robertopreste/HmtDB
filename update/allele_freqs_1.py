#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
import os
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import argparse
import logging
import pandas as pd
from app import db
from app.site.models import NtVariability
from app.static.dbdata import last_update

parser = argparse.ArgumentParser(description="""Modulo per il calcolo delle frequenze alleliche 
    sani e pazienti e continente-specifiche.""")
parser.add_argument("-f", "--from_date", dest="from_date", default=last_update, type=str,
                    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi 
                    pubblicati. Formato mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo 
                    update del database.""")
parser.add_argument("-H", "--healthy", action="store_true", dest="on_healthy",
                    help="Launch the program on healthy genomes (default: False).")
parser.add_argument("-C", "--complete", action="store_true", dest="on_complete",
                    help="Launch the program on complete genomes (default: False).")
parser.add_argument("-l", "--local", action="store_true", dest="local",
                    help="Do not work on db, but only create local files (default: False).",
                    default=False)

args = parser.parse_args()

logging.basicConfig(filename="log_allele_freqs_{}.log".format(args.from_date), format="%(levelname)s\t%(asctime)s: %(message)s", datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)


QUERY_TOTAL = "SELECT (SELECT haplotypeHmdb FROM Genome WHERE genomeId = TBAlign.genomeId) AS haplotype, alignment FROM GenAlignment as TBAlign WHERE genomeId IN (SELECT genomeId FROM Genome WHERE genomeType = '{}' {})"
QUERY_CONT = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE '{}_%%' {})"
QUERY_CONT_PA = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE GENOMEID = 0 UNION SELECT (SELECT haplotypeHmdb FROM GENOME WHERE GENOMEID = TBL_ALIGN.GENOMEID) AS HAPLOTYPE, ALIGNMENT FROM GenAlignment AS TBL_ALIGN WHERE GENOMEID IN (SELECT GENOMEID FROM GENOME WHERE haplotypeHmdb LIKE 'PA_{}_%%' {})"
QUERY_REF = "SELECT 'EU_XX_0000' AS 'HAPLOTYPE', ALIGNMENT FROM GenAlignment WHERE genomeId = 0"

base_path = os.path.dirname(os.path.abspath("."))  # HmtDB_flask/
algs_path = os.path.join(base_path, "update", "tmp_algs")
freq_path = os.path.join(base_path, "update", "all_freqs")
if not os.path.isdir(algs_path):
    os.mkdir(algs_path)
if not os.path.isdir(freq_path):
    os.mkdir(freq_path)


def dump_algs(genome_type, complete_genome):
    """Download alignments from the database and save them as csv files."""

    logging.info("Saving {} alignment to {}...".format(genome_type, algs_path))
    os.system("sqlite3 -header -csv {} \"{}\" > {}".format(os.path.join(base_path, "hmtdb.db"),
                                                       QUERY_TOTAL.format(genome_type, complete_genome),
                                                       os.path.join(algs_path, "alg_tot_{}.csv".format(genome_type))))
    logging.info("Done.")

    logging.info("Saving {} continent-specific alignments to {}...".format(genome_type, algs_path))
    for cont in ["AF", "AM", "AS", "EU", "OC"]:
        logging.info("Continent {}...".format(cont))
        if genome_type == "N":
            os.system("sqlite3 -header -csv {} \"{}\" > {}".format(os.path.join(base_path, "hmtdb.db"),
                                                               QUERY_CONT.format(cont, complete_genome),
                                                               os.path.join(algs_path, "alg_{}_{}.csv".format(cont, genome_type))))
        else:
            os.system("sqlite3 -header -csv {} \"{}\" > {}".format(os.path.join(base_path, "hmtdb.db"),
                                                               QUERY_CONT_PA.format(cont, complete_genome),
                                                               os.path.join(algs_path, "alg_{}_{}.csv".format(cont, genome_type))))
    logging.info("Done.")

    logging.info("Saving reference genome alignment...")
    os.system("sqlite3 -header -csv {} \"{}\" > {}".format(os.path.join(base_path, "hmtdb.db"),
                                                       QUERY_REF,
                                                       os.path.join(algs_path, "refj.csv")))
    logging.info("Done.")


def convert_seq_to_df(genome_type, continent):
    """Manipulate the dumped csv alignments and convert them to a dataframe in which each nt
    position is reported as a single column."""

    logging.info("Loading {} {} sequences...".format(continent, genome_type))
    df = pd.read_csv(os.path.join(algs_path, "alg_{}_{}.csv".format(continent, genome_type)),
                     names=["id", "alg_seq"], skiprows=1)
    logging.info("Done.")
    logging.info("Loading reference genome alignment...")
    refj = pd.read_csv(os.path.join(algs_path, "refj.csv"), names=["id", "alg_seq"], skiprows=1)
    logging.info("Done.")

    # Move the reference genome on top
    logging.info("Merging reference and aligned sequences...")
    df.drop(df[df["id"] == "EU_XX_0000"].index, inplace=True)
    df = pd.concat([refj, df], ignore_index=True)
    logging.info("Done.")

    # Create the reference positions
    logging.info("Creating reference positions...")
    indexes = []
    n = 0  # nucleotide position
    i = 0  # insertion position
    for el in list(refj.iloc[0]["alg_seq"]):
        if el != "-":  # non-gap position
            i = 0
            n += 1
        else:  # gap position
            i += 1
        indexes.append("{}.{}_{}".format(n, i, el))
    logging.info("Done. Saved {} indexes.".format(len(indexes)))

    # Conversion of sequences
    logging.info("Starting conversion of alignments to dataframe...")
    # fr_df = pd.DataFrame(columns=indexes)
    df.set_index("id", inplace=True)
    fr_df = df.alg_seq.str.split("", expand=True)  # .fillna("-")  # .drop(0, axis=1)
    fr_df.drop(0, axis=1, inplace=True)
    fr_df.fillna("-", inplace=True)
    # logging.info("Positions: {}, indexes: {}.".format(len(fr_df.columns), len(indexes)))
    if len(fr_df.columns) < len(indexes):
        indexes = indexes[:len(fr_df.columns)]
    else:
        while len(fr_df.columns) != len(indexes):
            fr_df.drop(fr_df.columns[len(fr_df.columns) - 1], axis=1, inplace=True)
    fr_df.columns = indexes

    # for row in df.set_index("id").itertuples():
    #     if len(tuple(list(row.alg_seq))) == len(indexes):
    #         nt_data = tuple(list(row.alg_seq))
    #     else:
    #         nt_data = list(row.alg_seq)
    #         while len(nt_data) < len(indexes):
    #             nt_data.append("-")
    #         nt_data = tuple(nt_data)
    #
    #     fr_df = pd.concat([fr_df, pd.DataFrame.from_records([nt_data], columns=indexes, index=[row.Index])])
    logging.info("Done.")

    fr_df.to_csv(os.path.join(algs_path, "fr_df_{}_{}.csv".format(continent, genome_type)), index=False)
    logging.info("Dataframe for {} {} sequences saved to {}.".format(continent, genome_type,
                                                                     os.path.join(algs_path,
                                                                                  "fr_df_{}_{}.csv".format(
                                                                                      continent,
                                                                                      genome_type))))


def update_freq(nt_pos, ins_pos, genome_type, cont, freq_val):
    """Update the db row(s) containing the given nucleotide and insertion position with the
    provided allele frequency. """

    q = NtVariability.query.filter(NtVariability.nucleotidePosition == int(nt_pos),
                                   NtVariability.insertionPosition == int(ins_pos),
                                   NtVariability.genomeType == genome_type).all()
    if len(q) > 0:
        for el in q:
            if cont == "AF":
                el.compVar_af = freq_val
            elif cont == "AM":
                el.compVar_am = freq_val
            elif cont == "AS":
                el.compVar_as = freq_val
            elif cont == "EU":
                el.compVar_eu = freq_val
            elif cont == "OC":
                el.compVar_oc = freq_val
            else:
                el.compVar_tot = freq_val
            db.session.commit()


def calculate_freqs(genome_type, continent):
    """Calculate allele frequencies from the dataframes created by the convert_seq_to_df() function."""

    df = pd.read_csv(os.path.join(algs_path, "fr_df_{}_{}.csv".format(continent, genome_type)))
    var_df = pd.DataFrame(columns=["A", "C", "G", "T", "gap", "oth"])

    logging.info("Starting allele frequencies calculation for {} {} sequences...".format(continent, genome_type))
    for col in df.columns:
        freqs = df[col].value_counts(normalize=True)
        freq_A = freqs.get("A", 0.0)
        freq_C = freqs.get("C", 0.0)
        freq_G = freqs.get("G", 0.0)
        freq_T = freqs.get("T", 0.0)
        freq_gap = freqs.get("-", 0.0)
        freq_oth = 1.0 - (freq_A + freq_C + freq_G + freq_T + freq_gap)
        var_df = pd.concat([var_df, pd.DataFrame({"A": freq_A, "C": freq_C, "G": freq_G, "T": freq_T,
                                                  "gap": freq_gap, "oth": freq_oth}, index=[col])])
    logging.info("Done.")

    var_df.reset_index(inplace=True)
    var_df.rename({"index": "position"}, axis=1, inplace=True)

    # Converting columns for compliance with the old sitevar
    var_df["site"] = var_df["position"].str.split("_", expand=True)[0]
    var_df["nucl_type"] = var_df["position"].str.split("_", expand=True)[1]
    var_df = var_df[["site", "nucl_type", "A", "C", "G", "T", "gap", "oth"]]

    logging.info("Saving final csv...")
    if args.on_complete:
        final_path = os.path.join(freq_path, "sitevar_{}_{}_compl.csv".format(continent, genome_type))
    else:
        final_path = os.path.join(freq_path, "sitevar_{}_{}.csv".format(continent, genome_type))

    var_df.to_csv(final_path, index=False)
    logging.info("Allele frequencies for {} {} sequences saved to {}.".format(continent, genome_type,
                                                                              final_path))

    if not args.local:
        logging.info("Uploading allele frequencies for {} {} sequences on db...".format(continent, genome_type))
        for idx in var_df.index:
            freq_val = ";".join([str(var_df.loc[idx]["A"]), str(var_df.loc[idx]["C"]),
                                 str(var_df.loc[idx]["G"]), str(var_df.loc[idx]["T"]),
                                 str(var_df.loc[idx]["gap"] + var_df.loc[idx]["oth"])])
            nt_pos = var_df.loc[idx]["site"].split(".")[0]
            ins_pos = var_df.loc[idx]["site"].split(".")[1]
            update_freq(nt_pos, ins_pos, genome_type, continent, freq_val)

        logging.info("Done.")


def main():
    if args.on_healthy:
        genome_type = "N"
    else:
        genome_type = "P"

    complete_genome = ""
    if args.on_complete:
        complete_genome = "AND completeGenome = 'Y'"

    dump_algs(genome_type, complete_genome)

    for cont in ["tot", "AF", "AM", "AS", "EU", "OC"]:
        convert_seq_to_df(genome_type, cont)
        calculate_freqs(genome_type, cont)
    # convert_seq_to_df(genome_type, "tot")
    # calculate_freqs(genome_type, "tot")


if __name__ == '__main__':
    main()
