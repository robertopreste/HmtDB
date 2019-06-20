#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from Bio import Entrez, SeqIO, GenBank, Medline
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from app.static.dbdata import last_update
from app import db
from app.site import models
import csv
import os
import logging
import argparse
import pandas as pd
from app.site.models import Genome, Country, IndividualsData, Sources, GenAlignment, Reference, GenomeSnp, Insertion, Deletion


parser = argparse.ArgumentParser(description="""Modulo per eseguire il caricamento dei nuovi genomi 
    sul database. Se lanciato senza argomenti si presuppone che nello step precedente sia stata 
    usata la data dell'ultimo update del database, altrimenti è possibile specificare la data usata 
    per la ricerca.""")
parser.add_argument("-f", "--from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale sono stati ricercati i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")

args = parser.parse_args()

logging.basicConfig(filename="log_load_entries_{}.log".format(args.from_date),
                    format="%(levelname)s %(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)
last_update_month, last_update_year = args.from_date.split("_")

update_dir = "HmtDB_update_{}".format(args.from_date)
fasta_dir = "{}/Fasta_Results_{}".format(update_dir, args.from_date)
gb_dir = "{}/GenBank_Results_{}".format(update_dir, args.from_date)

query_country = "SELECT countryId, continentCode, countryCode FROM Country ORDER BY continentCode, countryCode"

query_genomes = "SELECT genomeID FROM Genome WHERE genomeType = '%s' AND individualId IN (SELECT individualId FROM IndividualsData WHERE countryId = '%d')"

query_alignments = "SELECT genomeId, alignment FROM GenAlignment WHERE genomeId IN (SELECT genomeId FROM Genome WHERE genomeType = '%s' AND genomeId <> 0 AND haplotypeHmdb <> 'PA_EU_XX_0001') ORDER BY genomeId"
# todo: forse questa non serve, viene usato solo per recuperare il numero di genomi allineati

insert_snp = "INSERT INTO GenomeSnp (genomeId, snpPosition, rcrsType, snpType) VALUES ('%d', '%d', '%s', '%s')"
insert_ins = "INSERT INTO Insertion (genomeId, position5P, sequence) VALUES ('%d', '%d', '%s')"
insert_del = "INSERT INTO Deletion (genomeId, fromPosition, toPosition) VALUES ('%d', '%d', '%d')"

last_id = int(Genome.query.order_by(Genome.genomeId.desc()).first().genomeId)
padding = "%%0%dd" % 4
# todo: questo padding si potrebbe sostituire con zfill()

country_dict_h = {}
country_dict_p = {}


def create_country_dict():

    query = query_country
    res = db.engine.execute(query)

    logging.info("Creating country dictionary...")

    for country_id, cont_code, country_code in res:
        # Normal genomes
        query_h = query_genomes % ("N", country_id)
        genomes_h = db.engine.execute(query_h)
        g_h = genomes_h.fetchall()

        country_dict_h["{}_{}".format(str(cont_code), str(country_code))] = len(g_h) if g_h else 0

        # Patient genomes
        query_p = query_genomes % ("P", country_id)
        genomes_p = db.engine.execute(query_p)
        g_p = genomes_p.fetchall()

        country_dict_p["PA_{}_{}".format(str(cont_code), str(country_code))] = len(g_p) if g_p else 0

    logging.info("Done. ")

    return


def get_snps(rif, inc):

    pos_a = 0
    n_gaps = 0
    alg_len = len(rif)
    gap = "-"
    l_snps = []
    l_ins = []
    l_dels = []

    while pos_a < alg_len:
        x = rif[pos_a]
        y = inc[pos_a]
        if x != y:
            if x != gap and y != gap:
                # SNP
                l_snps.append((pos_a - n_gaps + 1, x, y))
                pos_a += 1
            elif x == gap and y != gap:
                # Ins
                pos_i = pos_a - n_gaps
                ins_seq = [y]
                pos_a += 1
                n_gaps += 1
                x = rif[pos_a]
                y = inc[pos_a]
                while pos_a < alg_len - 1 and ((x == gap and y != gap) or (x == y == gap)):
                    if y != gap:
                        ins_seq.append(y)
                    pos_a += 1
                    n_gaps += 1
                    x = rif[pos_a]
                    y = inc[pos_a]
                if pos_a == alg_len - 1:
                    pos_a += 1
                l_ins.append((pos_i, "".join(ins_seq)))
            elif x != gap and y == gap:
                # Del
                pos_d = pos_a - n_gaps + 1
                pos_a += 1
                if pos_a < alg_len:
                    x = rif[pos_a]
                    y = inc[pos_a]
                    while pos_a < alg_len - 1 and ((x != gap and y == gap) or (x == y == gap)):
                        pos_a += 1
                        if x == y == gap:
                            n_gaps += 1
                        x = rif[pos_a]
                        y = inc[pos_a]
                    if pos_a == alg_len - 1:
                        pos_a += 1
                l_dels.append((pos_d, pos_a - n_gaps))
        else:
            # Permanenza, SNP presente nel riferimento
            # debug: questa la tolgo dato che lista_perm sembra sia sempre vuota!
            # if lista_perm:
            #     if pos_a - n_gaps + 1 in lista_perm:
            #         l_snps.append((pos_a - n_gaps + 1, x, y))
            pos_a += 1

            if x == gap:
                n_gaps += 1

    return l_snps, l_ins, l_dels


def load_genomes(gen_type):

    logging.info("Updating db with genomes of type '{}'...".format(gen_type))

    # Tabella GenBank
    if gen_type == "N":
        gb_table = pd.read_csv("{}/GenBank_Results_{}_Norm.csv".format(update_dir, args.from_date),
                               index_col="Counter", engine="python")
    else:
        gb_table = pd.read_csv("{}/GenBank_Results_{}_Pat.csv".format(update_dir, args.from_date),
                               index_col="Counter", engine="python")

    # Tabella PubMed
    ref_table = pd.read_csv("{}/PubMed_Results_{}.csv".format(update_dir, args.from_date),
                            index_col="Counter", engine="python")

    gb_table.fillna("", inplace=True)  # rimuovo gli NaN e li sostituisco con ""
    ref_table.fillna("", inplace=True)

    for n, row in enumerate(gb_table.itertuples(), start=1):

        logging.info("Collecting information for genome {} (HmtDB id {})".format(row.Locus, row.New_Id))

        # Genome Sequence dal relativo file fasta
        with open("{}/{}.fasta".format(fasta_dir, row.Locus)) as f:
            l = f.readlines()[1:]
            s = "".join([line.rstrip("\n") for line in l])
        genome_sequence = s

        # Complete Genome
        # debug: controllo se le prime o ultime 200 basi siano tutte N, in quel caso si dice parziale
        if genome_sequence[:200] == "N" * 200 or genome_sequence[-200:] == "N" * 200:
            complete_genome = "N"
        else:
            complete_genome = row.Gen_Complete

        # Positions
        start_pos = row.Start
        end_pos = row.End

        # Db di provenienza (solo Genbank?)
        reference_db = "Genbank"
        reference_db_id = row.Locus

        # Normal or Patient
        genome_type = "N" if row.Gen_Healthy == "Y" else "P"

        # Genome Id
        # genome_id = int(last_id + row.Index)
        genome_id = int(row.New_Id)

        # Alignment
        # alignment_id = GenAlignment.query.filter(GenAlignment.genomeId == genome_id).first().alignmentId

        # References
        ref_row = ref_table.loc[row.Index]
        ref_last_id = int(Reference.query.order_by(Reference.referenceId.desc()).first().referenceId)
        ref_pubmed_id = ref_row.PubmedId
        ref_authors = ref_row.Authors
        ref_title = ref_row.Title.decode("utf-8")
        ref_volume = ref_row.Volume if ref_row.Volume else None
        ref_year = int(ref_row.Year) if ref_row.Year else None
        ref_paper = ref_row.Paper
        ref_issue = ref_row.Issue if ref_row.Issue else None

        # Country
        # todo: rivedere questa routine sfruttando le nuove Iso Countries
        country = row.Country if row.Country else None
        if country:
            res = Country.query.filter(Country.countryName == country).first()
            country_id = res.countryId
        # print country_id, country
        # if gen_type == "N":
        #     # gen_num = country_dict_h["%s_%s" % (res.continentCode, res.countryCode)] + n
        #     gen_num = country_dict_h["%s_%s" % (res.continentCode, res.countryCode)] + 1
        #     haplotype_hmdb = "%s_%s_%s" % (res.continentCode, res.countryCode, padding % gen_num)
        #     country_dict_h["%s_%s" % (res.continentCode, res.countryCode)] += 1
        # else:
        #     # gen_num = country_dict_p["PA_%s_%s" % (res.continentCode, res.countryCode)] + n
        #     gen_num = country_dict_p["PA_%s_%s" % (res.continentCode, res.countryCode)] + 1
        #     haplotype_hmdb = "PA_%s_%s_%s" % (res.continentCode, res.countryCode, padding % gen_num)
        #     country_dict_p["PA_%s_%s" % (res.continentCode, res.countryCode)] += 1
        # print gen_num, haplotype_hmdb
        # if res:  # se la nazione esiste già nella lista
        #     country_id = res.countryId
            if gen_type == "N":
                # gen_num = country_dict_h["%s_%s" % (res.continentCode, res.countryCode)] + n
                gen_num = country_dict_h["{}_{}".format(res.continentCode, res.countryCode)] + 1
                haplotype_hmdb = "{}_{}_{}".format(res.continentCode, res.countryCode,
                                                   padding % gen_num)
                country_dict_h["{}_{}".format(res.continentCode, res.countryCode)] += 1
            else:
                # gen_num = country_dict_p["PA_%s_%s" % (res.continentCode, res.countryCode)] + n
                gen_num = country_dict_p["PA_{}_{}".format(res.continentCode, res.countryCode)] + 1
                haplotype_hmdb = "PA_{}_{}_{}".format(res.continentCode, res.countryCode,
                                                      padding % gen_num)
                country_dict_p["PA_{}_{}".format(res.continentCode, res.countryCode)] += 1
        else:  # se la nazione è NULL
            country_id = 210
            if gen_type == "N":
                gen_num = country_dict_h["XX_XX"] + 1
                haplotype_hmdb = "XX_XX_{}".format(padding % gen_num)
                country_dict_h["XX_XX"] += 1
            else:
                # gen_num = country_dict_p["PA_%s_%s" % (res.continentCode, res.countryCode)] + n
                gen_num = country_dict_p["PA_XX_XX"] + 1
                haplotype_hmdb = "PA_XX_XX_{}".format(padding % gen_num)
                country_dict_p["PA_XX_XX"] += 1
            # country_id = int(Country.query.order_by(Country.countryId.desc()).first().countryId) + 1
            # alla nuova country viene assegnato un nuovo id
            # todo: qui bisogna anche controllare come creare l'haplotype hmdb se la country è nuova, poi si pensa
            # todo: questo in realtà si può fare direttamente sfruttando sqlalchemy, che dovrebbe assegnare già un nuovo id con autoincrement
            # c = Country(countryId=country_id, countryName=country)
            # todo: bisogna trovare un modo per inserire automaticamente le country se non sono già nel db (servono il countryCode, il continente di appartenenza e il suo code

        # Individual data
        sex = row.Sex.strip('"') if row.Sex else None
        tissue_type = row.Tissue_Type.strip('"') if row.Tissue_Type else None
        devstage = row.Devstage.strip('"') if row.Devstage else None

        # Isolation source
        isolate = row.Isolate
        isolation_source = row.Isolation_Source

        if isolation_source != "":
            # se la source esiste già nella tabella Sources del db
            if Sources.query.filter(Sources.sourceName == isolation_source).first():
                source_id = int(Sources.query.filter(Sources.sourceName == isolation_source).first().sourceId)
            else:  # se si tratta invece di una nuova source
                source_id = int(Sources.query.filter(Sources.sourceName != "").order_by(Sources.sourceId.desc()).first().sourceId)
        else:  # se la source è vuota, viene dato il sourceId 299
            source_id = 299

        # Haplos
        haplotype_user = row.Haplotype.strip('"') if row.Haplotype else None
        haplogroup_user = row.Haplogroup.strip('"') if row.Haplogroup else None

        # New Reference entry
        # referenceId=ref_last_id + int(row.Index),
        logging.info("Updating Reference table for genome id {}...".format(genome_id))
        if ref_pubmed_id == "Unpublished":
            reference_Entry = Reference(pubmedId=None, author=None, title=None, volume=None,
                                        year=None, paper=None, issue=None, genomeId=genome_id)
            db.session.add(reference_Entry)
        else:
            reference_Entry = Reference(pubmedId=ref_pubmed_id, author=ref_authors, title=ref_title,
                                        volume=ref_volume, year=ref_year, paper=ref_paper,
                                        issue=ref_issue, genomeId=genome_id)
            db.session.add(reference_Entry)
        logging.info("Done.")

        # New Genome entry
        logging.info("Updating Genome table for genome id {}...".format(genome_id))
        genome_Entry = Genome(genomeId=genome_id, genomeSequence=str(genome_sequence),
                              completeGenome=complete_genome, startPosition=int(start_pos),
                              endPosition=int(end_pos), haplotypeUser=haplotype_user,
                              haplogroupUser=haplogroup_user, haplotypeHmdb=haplotype_hmdb,
                              referenceDb=reference_db, referenceDbId=reference_db_id,
                              genomeType=genome_type)
                              # alignmentId=alignment_id, sourceId=source_id)
                              # references=(ref_last_id + int(row.Index)))

        db.session.add(genome_Entry)
        logging.info("Done.")

        # New IndividualsData entry
        # todo: nell'IndividualsData dovrei inserire age, ma da dove la prendo?
        logging.info("Updating IndividualsData table for genome id {}...".format(genome_id))
        individuals_Entry = IndividualsData(individualId=genome_id, sex=sex,
                                            individualType=genome_type, countryId=country_id,
                                            genomeId=[genome_Entry])
                                            # genomeId=genome_id)
                                            # individualType=genome_type, countryId=country_id)

        db.session.add(individuals_Entry)
        logging.info("Done. ")

        # New GenomeSnp entry
        # expl: ref_alg sarebbe l'allineamento della sequenza di riferimento
        logging.info("Collecting SNPs and updating GenomeSnp table for genome id {}...".format(genome_id))
        ref_alg = str(GenAlignment.query.filter(GenAlignment.genomeId == 0).first().alignment)
        inc_alg = str(GenAlignment.query.filter(GenAlignment.genomeId == genome_id).first().alignment)
        l_snps, l_ins, l_dels = get_snps(ref_alg, inc_alg)

        # SNPs
        to_add = []
        for pos, ref_nt, snp_nt in l_snps:
            to_add.append(GenomeSnp(snpPosition=pos, snpType=snp_nt, rcrsType=ref_nt,
                                    genomeId=genome_id))
        db.session.add_all(to_add)

        # Ins
        to_add = []
        for pos, ins_seq in l_ins:
            to_add.append(Insertion(position5P=pos, sequence=ins_seq, genomeId=genome_id))
        db.session.add_all(to_add)

        # Dels
        to_add = []
        for f_pos, t_pos in l_dels:
            to_add.append(Deletion(fromPosition=f_pos, toPosition=t_pos, genomeId=genome_id))
        db.session.add_all(to_add)

        logging.info("Done. ")

    db.session.commit()
    logging.info("Done. ")


def main():
    create_country_dict()
    load_genomes("N")
    load_genomes("P")


if __name__ == '__main__':
    main()





