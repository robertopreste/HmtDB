#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from Bio import Entrez, SeqIO, GenBank, Medline
import sys
import csv
import os
import logging
import argparse
import time
import pandas as pd
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from app.static.dbdata import last_update
from app import db
from app.site.models import Genome, Country


parser = argparse.ArgumentParser(description="""Modulo per eseguire la ricerca di nuovi genomi, 
    il loro download e infine il parsing GenBank e MedLine. Se lanciato senza argomenti verranno 
    cercati i nuovi genomi pubblicati a partire dalla data dell'ultimo aggiornamento del database, 
    altrimenti è possibile specificare una data per la ricerca.""")
parser.add_argument("phase", type=int, help="""Fase del processo di ricerca nuovi genomi. Nella 
    fase 1, vengono cercati i nuovi genomi pubblicati dalla data specificata in -f (o da quella di 
    default); terminata questa fase vanno controllati i nuovi genomi trovati per inserire nel csv 
    GenBank temporaneo se ogni genoma è normale o paziente. Nella fase 2 queste informazioni vengono 
    usate per generare il multifasta con l'id corretto, che viene anche aggiunto nel csv GenBank 
    definitivo.""")
parser.add_argument("-f", "--from_date", dest="from_date", default=last_update, type=str,
    help="""Data a partire dalla quale ricercare i nuovi genomi pubblicati. 
    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del database.""")

args = parser.parse_args()

Entrez.email = "hmtdb.update@uniba.it"
Entrez.tool = "HmtDB updating protocol"
last_update_month, last_update_year = args.from_date.split("_")

update_dir = "HmtDB_update_{}".format(args.from_date)
fasta_dir = "{}/Fasta_Results_{}".format(update_dir, args.from_date)
gb_dir = "{}/GenBank_Results_{}".format(update_dir, args.from_date)
#pm_dir = "MedLine_Results_%s/" % last_update
id_list = []
id_set = set()
# todo: questi 2 file servono veramente? per ora li tolgo
# id_list_file = open("%s/id_list.txt" % update_dir, "a")
# pid_list_file = open("%s/pid_list.txt" % update_dir, "a")


def country_common_names(country):
    """Controlla se esistono altri nomi per country diversi da quello ufficiale. In caso positivo,
    viene restituito il nome ufficiale, altrimenti ritorna False.
    Da aggiornare periodicamente sulla base del contenuto del file notes.txt."""

    coun_dict = {
        "Russia": "Russian Federation",
        "Palestine": "State of Palestine",
        "United Kingdom": "United Kingdom of Great Britain and Northern Ireland",
        "Czech Republic": "Czechia",
        "USA": "United States of America",
        "Iran": "Islamic Republic of Iran",
        "Syria": "Syrian Arab Republic",
        "Tanzania": "United Republic of Tanzania",
        "Borneo": "Brunei Darussalam",
        "Venezuela": "Bolivian Republic of Venezuela"
    }

    if country in coun_dict:
        return coun_dict[country]
    else:
        return False


def search_new_genomes():
    """Ricerca i nuovi genomi mitocondriali umani pubblicati dalla data dell'ultimo aggiornamento.
    Ritorna una lista di Accession Number."""

    logging.info("Using the standard Entrez query to search for new genomes...")
    query_string = """((14000[SLEN] : 17000[SLEN]) AND "Homo sapiens"[Organism] AND mitochondrion[Filter]) NOT ("Homo sapiens ssp. Denisova"[Organism] OR ("Homo sapiens ssp. Denisova"[Organism] OR Denisova hominin[All Fields]) OR neanderthalensis[All Fields]) NOT "1980/01/01"[PDAT] : "{}/{}/01"[PDAT] NOT NC_012920.1[Accession]"""
    id_search = Entrez.esearch(db="nucleotide",
                               term=query_string.format(last_update_year, last_update_month),
                               retmax=100000,
                               idtype="acc")
    id_results = Entrez.read(id_search)["IdList"]
    # purtroppo bisogna splittare i risultati che sono nel nuovo formato AC0011223.1
    # bisognerebbe risolvere questo problema, magari aggiornando i dati nel database
    results = set()
    for el in id_results:
        # todo: controllare se è .2 bisogna aggiornare quello già presente nel db
        a, b = el.split(".")
        results.add(a)  # aggiungo solo il primo elemento, perché il secondo è semplicemente "1"
    results = list(results)
    logging.info("Found {} genomes. Their Accession Numbers are {}".format(len(results),
                                                                           [i for i in results]))

    # cerca se gli Accession Number trovati sono già nel database, e in questo caso li elimina dalla lista
    # debug: questo lo elimino per risparmiare tempo sul cazzo di recas...
    # for el in results:
    #     logging.info("Checking GenBankId {}...".format(el))
    #     if Genome.query.filter(Genome.referenceDbId == el).first():
    #         logging.info("Removing GenBankId {}...".format(el))
    #         results.remove(el)
    # debug: lo cambio così, dovrebbe essere più veloce
    logging.info("Retrieving GenBank Ids from the database...")
    total_q = pd.read_sql("SELECT referenceDbId FROM Genome", "sqlite:///../hmtdb.db")
    total_ids = total_q.referenceDbId.tolist()
    logging.info("Done.")
    logging.info("Checking for genomes already in the db...")
    for el in results:
        logging.info("Checking GenBank Id {}...".format(el))
        if el in total_ids:
            logging.info("Removing GenBank Id {}...".format(el))
            results.remove(el)
            logging.info("Done.")
    logging.info("Done.")

    # Cerca e salva in singoli file le entries trovate sia in formato GenBank che FASTA
    # qui non serve salvare tutti i file
    # basta salvarne uno alla volta, parsarlo (magari usando Biopython) e poi cancellarlo, e passare al seguente
    return results


def parse_genbank(gb_file, counter, note_file):
    """Effettua il parsing di un file GenBank. Ritorna una lista di due elementi, in cui il primo
    è così composto: [locus, starts, ends, isolate, isol_source, haplotype, haplogroup, sex,
    tissue_type, devstage, country] e il secondo contiene il pubmed id associato.
    Se i campi non vengono trovati nel file GenBank, il default inserito è None (campo vuoto)."""

    record = GenBank.read(gb_file)

    data_list = []

    locus = starts = ends = isolate = isol_source = haplotype = haplogroup = sex = tissue_type \
        = devstage = country = gen_complete = pubmed = None

    if record.locus:
        locus = record.locus

    if record.features[0].location:
        starts, ends = record.features[0].location.split("..")
        if (int(ends) - int(starts)) < 16000:
            gen_complete = "N"
        else:
            gen_complete = "Y"

    for feat in record.features[0].qualifiers:

        if feat.key == "/isolate=":
            isolate = feat.value.strip('"')
        if feat.key == "/isolation_source=":
            isol_source = feat.value.strip('"')
        if feat.key == "/haplotype=":
            haplotype = feat.value.strip('"')
        if feat.key == "/haplogroup=":
            haplogroup = feat.value.strip('"')
        if feat.key == "/sex=":
            if feat.value == "male":
                sex = "M"
            elif feat.value == "female":
                sex = "F"
            else:
                sex = feat.value.strip('"')
        if feat.key == "/tissue_type=":
            tissue_type = feat.value.strip('"')
        if feat.key == "/devstage=":
            devstage = feat.value.strip('"')
        if feat.key == "/country=":
            country = feat.value.strip('"')
            if ":" in country:
                country = country.split(":")[0]
            db_country = Country.query.filter(Country.countryName == country).first()
            if not db_country:
                if country_common_names(country):
                    country = country_common_names(country)
                else:
                    note_file.write("Country {} not in db! (GenBank file {}.gb)\n".format(country,
                                                                                          locus))

    for ref in record.references:

        if ref.pubmed_id:  # pubblicazione vera e non direct submission
            pubmed = ref.pubmed_id
            break
        else:
            continue

    if not pubmed:
        pubmed = "Unpublished"

    data_list.append([counter, locus, starts, ends, isolate, isol_source, haplotype, haplogroup,
                      sex, tissue_type, devstage, country, gen_complete, ""])
    # debug: quest'ultimo valore va poi completato a mano in base allo stato del genoma (healthy or patient, rispettivamente Y o N)
    data_list.append(str(pubmed))  # non si sa mai esce fuori un int

    return data_list


def parse_medline(record, counter):
    """Effettua il parsing di un file Medline. Ritorna una lista di elementi, così composta:
    [pubmedid, first author (et al), title, volume, year, paper, issue, counter (per associare
    l'entry a quella genbank].
    I campi first_page e abstract sono stati eliminati perché inutili, ed il campo authors in
    realtà riporta solo il primo autore.
    Se i campi non vengono trovati nel file Medline, il default inserito è None (campo vuoto)."""

    # test: rimuovo questo sotto e gli do direttamente l'elemento Medline.read()
    # record = Medline.read(pm_file)

    pubmed = authors = title = volume = year = paper = issue = None

    try:
        pubmed = record["PMID"]
    except:
        pubmed = ""

    try:
        authors = record["AU"][0] + " et al"  # prendo solo il primo autore
    except:
        authors = ""

    try:
        title = record["TI"]
    except:
        title = ""

    try:
        volume = record["VI"]
    except:
        volume = ""

    try:
        year = record["DP"].split(" ")[0]
    except:
        year = ""
    #elif record["EDAT"]:
    #    year = record["EDAT"]

    try:
        paper = record["TA"]
    except:
        paper = ""
    #elif record["TA"]:
    #    paper = record["TA"]

    try:
        issue = record["IP"]
    except:
        issue = ""

    data_list = [pubmed, authors, title, volume, year, paper, issue, counter]

    return data_list


def run_phase1():
    """Wrapper per le funzioni search_new_genomes, parse_genbank e parse_medline.
    Vengono prima cercati i nuovi genomi in Entrez, quindi scaricati i relativi file GenBank e
    FASTA, infine ogni file GenBank viene parsato per trovarne le informazioni e le referenze,
    rispettivamente salvate in un csv GenBank e in un csv PubMed."""

    os.system("if [ -d {} ]; then rm -Rf {}; fi".format(update_dir, update_dir))

    os.system("mkdir -p {}".format(update_dir))
    os.system("mkdir -p {}".format(fasta_dir))
    os.system("mkdir -p {}".format(gb_dir))
    #os.mkdir(pm_dir)

    # expl: resetto tutti questi file prima di tutto, usando il with context

    # Tabella con i dati GenBank dei nuovi genomi trovati
    with open("{}/GenBank_Results_{}.csv".format(update_dir, args.from_date), "wb") as gb_table:
        writer_gb = csv.writer(gb_table)
        field_names_gb = ["Counter", "Locus", "Start", "End", "Isolate", "Isolation_Source",
                          "Haplotype", "Haplogroup", "Sex", "Tissue_Type", "Devstage", "Country",
                          "Gen_Complete", "Gen_Healthy"]
        writer_gb.writerow(field_names_gb)
    gb_table = open("{}/GenBank_Results_{}.csv".format(update_dir, args.from_date), "ab")
    writer_gb = csv.writer(gb_table)

    # Tabella con i dati PubMed dei nuovi genomi trovati
    with open("{}/PubMed_Results_{}.csv".format(update_dir, args.from_date), "wb") as pm_table:
        writer_pm = csv.writer(pm_table)
        field_names_pm = ["PubmedId", "Authors", "Title", "Volume", "Year", "Paper", "Issue",
                          "Counter"]
        writer_pm.writerow(field_names_pm)
    pm_table = open("{}/PubMed_Results_{}.csv".format(update_dir, args.from_date), "ab")
    writer_pm = csv.writer(pm_table)

    logging.info("Searching and downloading new genomes...")
    # results = search_new_genomes()
    results_set = set(search_new_genomes())

    with open("{}/notes.txt".format(update_dir), "w") as notes_file:
        pass
    notes_file = open("{}/notes.txt".format(update_dir), "a")

    with open("{}/new_hmtdb_ids.txt".format(update_dir), "w") as new_hmtdb_file:
        pass

    with open("{}/pid_list.txt".format(update_dir), "w") as pid_list_file:
        pass
    pid_list_file = open("{}/pid_list.txt".format(update_dir), "a")

    last_id = int(Genome.query.order_by(Genome.genomeId.desc()).first().genomeId)

    # debug: qui abbiamo la lista di accession id da aggiungere al db
    for n, i in enumerate(results_set, start=1):
        logging.info("Saving {}...".format(i))

        # GenBank
        # Cerca e scarica il file .gb di ogni nuovo genoma trovato
        for num_retry in range(5):
            try:
                logging.info("Retrieving GenBank {}...".format(i))
                entry_gb = Entrez.efetch(db="nucleotide", id=i, retmax=100000, rettype="gb",
                                         retmode="text").read()
                assert (entry_gb != ""), "GenBank file for id {} is empty!".format(i)
                logging.info("Done. ")
                break
            except AssertionError:
                logging.error("Entry was empty. Retry {} in 30 seconds...".format(num_retry))
                time.sleep(30)
            except Exception:
                logging.error("There was an error. Retry {} in 30 seconds...".format(num_retry))
                time.sleep(30)
        else:
            logging.error("Could not retrieve GenBank {}. Stopping.".format(i))
            raise
        with open("{}/{}/{}.gb".format(os.getcwd(), gb_dir, i), "wb") as end_gb:
            end_gb.write(entry_gb)
        # entry_gb.close()

        # Fasta
        # debug: qui viene recuperato il fasta di ogni nuovo genoma trovato, e unito nel multifasta
        for num_retry in range(5):
            try:
                logging.info("Retrieving fasta from GenBank {}...".format(i))
                entry_fa = Entrez.efetch(db="nucleotide", id=i, retmax=100000, rettype="fasta",
                                         retmode="text").read()
                assert (entry_fa != ""), "Fasta file for id {} is empty!".format(i)
                logging.info("Done. ")
                break
            except AssertionError:
                logging.error("Entry was empty. Retry {} in 30 seconds...".format(num_retry))
                time.sleep(30)
            except Exception:
                logging.error("There was an error. Retry {} in 30 seconds...".format(num_retry))
                time.sleep(30)
        else:
            logging.error("Could not retrieve fasta from GenBank {}. Stopping.".format(i))
            raise
        with open("{}/{}/multi_seqs.fasta".format(os.getcwd(), fasta_dir), "a") as end_multi, \
                open("{}/{}/new_hmtdb_ids.txt".format(os.getcwd(), update_dir), "a") as new_hmtdb_file:
            with open("{}/{}/{}.fasta".format(os.getcwd(), fasta_dir, i), "wb") as end_fa:
                end_fa.write(entry_fa)
            logging.info("Adding Fasta sequence into multifasta file...")
            print "{}/{}.fasta".format(fasta_dir, i)
            record = SeqIO.read("{}/{}/{}.fasta".format(os.getcwd(), fasta_dir, i), "fasta")
            new_id = last_id + n
            end_multi.write(">" + str(new_id) + "\n")
            end_multi.write(str(record.seq) + "\n")
            logging.info("Entry {} has now id {}.".format(i, new_id))
            # salvo il nuovo id HmtDB per usarlo nella fase dell'allineamento
            new_hmtdb_file.write(str(new_id) + "\n")

        # entry_fa.close()
        id_list.append("{}".format(i))

        # Parsing GenBank
        # Fa il parsing del file GenBank scaricato poco sopra, e salva i risultati nei csv GenBank e PubMed
        with open("{}/{}/{}.gb".format(os.getcwd(), gb_dir, i), "rb") as gb_file:
            logging.info("Parsing GenBank file {}...".format(i))
            gb_parsed = parse_genbank(gb_file, n, notes_file)
            writer_gb.writerow(gb_parsed[0])  # scrive la riga nel csv genbank
            # Parsing Medline
            if gb_parsed[1] == "Unpublished":
                writer_pm.writerow(["Unpublished", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL",
                                    str(n)])
            else:
                for num_retry in range(5):
                    try:
                        logging.info("Retrieving PubmedId {}...".format(gb_parsed[1]))
                        # test: provo ad eseguire direttamente il Medline.read() dell'entry e fare l'assert su questo
                        entry_pm = Entrez.efetch(db="pubmed", id=gb_parsed[1], rettype="medline",
                                                 retmode="text")
                        pm_record = Medline.read(entry_pm)
                        assert (pm_record != ""), "Pubmed file for id {} is empty!".format(gb_parsed[1])
                        logging.info("Done. ")
                        break
                    except AssertionError:
                        logging.error("Entry was empty. Retry {} in 30 seconds...".format(num_retry))
                        time.sleep(30)
                    except Exception:
                        logging.error("There was an error. Retry {} in 30 seconds...".format(num_retry))
                        time.sleep(30)
                else:
                    logging.error("Could not retrieve PubmedId {}. Stopping.".format(gb_parsed[1]))
                    raise
                # debug: fa direttamente il parsing dell'entry Medline, senza prima salvarla in un file
                #with open("%s/%s%s.medline" % (os.getcwd(), pm_dir, gb_parsed[1]), "wb") as end_pm:
                #    end_pm.write(entry_pm.read())
                #with open("%s/%s%s.medline" % (os.getcwd(), pm_dir, gb_parsed[1]), "rb") as pm_file:
                #    pm_parsed = parse_medline(pm_file, n)
                #    writer_pm.writerow(pm_parsed)
                logging.info("Parsing Medline entry {}...".format(gb_parsed[1]))
                # todo: questa è una ripetizione dell'assert fatto sopra, bisogna sistemarlo
                # for num_retry in range(5):
                #     try:
                #         logging.info("Retrieving PubmedId {}...".format(gb_parsed[1]))
                #         entry_pm = Entrez.efetch(db="pubmed", id=gb_parsed[1], rettype="medline",
                #                                  retmode="text")
                #         break
                #     except Exception:
                #         logging.error("There was an error. Retry {} in 30 seconds...".format(num_retry))
                #         time.sleep(30)
                # else:
                #     logging.error("There was an error. Retry {} in 30 seconds...".format(num_retry))
                #     raise
                # pm_parsed = parse_medline(entry_pm, n)
                pm_parsed = parse_medline(pm_record, n)
                writer_pm.writerow(pm_parsed)
                # entry_pm.close()
                pid_list_file.write("{}\n".format(gb_parsed[1]))
        logging.info("Done.")

    logging.info("Retrieval and parsing of GenBank and Medline file complete.\n")
    gb_table.close()
    pm_table.close()
    notes_file.close()
    pid_list_file.close()

    logging.info("Saving new genomes' GenBank and HmtDB IDs...")
    with open("{}/new_gb_ids.txt".format(update_dir), "w") as new_gb_file:
        for el in results_set:
            new_gb_file.write(el + "\n")
    logging.info("Complete. Please check the new genomes and then relaunch this script in phase 2.")


def run_phase2():
    """Il file csv GenBank in cui sono state inserite le informazioni sui genomi sani/pazienti
    viene usato per generare gli id HmtDB per creare il multifasta da allineare."""

    last_id = int(Genome.query.order_by(Genome.genomeId.desc()).first().genomeId)

    logging.info("Loading list of new genomes found...")
    with open("{}/new_gb_ids.txt".format(update_dir), "rb") as new_gb_file:
        results = [el.strip() for el in new_gb_file]

    # prendo l'informazione se ogni genoma sia sano o meno
    health_dict = {}
    gb_file = pd.read_csv("{}/GenBank_Results_{}.csv".format(update_dir, args.from_date),
                          index_col="Locus", engine="python")
    for el in gb_file.index:
        health_dict[el] = gb_file["Gen_Healthy"][el]

    logging.info("Merging Fasta files in one single multifasta...")
    # with open("%s/multi_seqs.fasta" % fasta_dir, "a") as end_multi, \
    #         open("%s/new_hmtdb_ids.txt" % update_dir, "w") as new_hmtdb_file, \
    with open("{}/GenBank_Results_{}_Norm.csv".format(update_dir, args.from_date), "wb") as gb_healthy, \
            open("{}/GenBank_Results_{}_Pat.csv".format(update_dir, args.from_date), "wb") as gb_patient:
        w_h = csv.writer(gb_healthy)
        w_p = csv.writer(gb_patient)
        logging.info("Saving healthy genomes in GenBank_Results_{}_Norm.csv...".format(args.from_date))
        w_h.writerow(["Counter", "Locus", "Start", "End", "Isolate", "Isolation_Source",
                      "Haplotype", "Haplogroup", "Sex", "Tissue_Type", "Devstage", "Country",
                      "Gen_Complete", "Gen_Healthy", "New_Id"])
        logging.info("Saving patient genomes in GenBank_Results_{}_Pat.csv...".format(args.from_date))
        w_p.writerow(["Counter", "Locus", "Start", "End", "Isolate", "Isolation_Source",
                      "Haplotype", "Haplogroup", "Sex", "Tissue_Type", "Devstage", "Country",
                      "Gen_Complete", "Gen_Healthy", "New_Id"])
        for n, el in enumerate(results, start=1):
            # record = SeqIO.read("%s/%s.fasta" % (fasta_dir, el), "fasta")
            line = gb_file.loc[el]
            new_id = last_id + n
            if health_dict[el] == "Y":
                w_h.writerow([line.Counter, el, line.Start, line.End, line.Isolate,
                              line.Isolation_Source, line.Haplotype, line.Haplogroup, line.Sex,
                              line.Tissue_Type, line.Devstage, line.Country, line.Gen_Complete,
                              line.Gen_Healthy, new_id])
            elif health_dict[el] == "N":
                w_p.writerow([line.Counter, el, line.Start, line.End, line.Isolate,
                              line.Isolation_Source, line.Haplotype, line.Haplogroup, line.Sex,
                              line.Tissue_Type, line.Devstage, line.Country, line.Gen_Complete,
                              line.Gen_Healthy, new_id])
            # end_multi.write(">" + str(new_id) + "\n")
            # end_multi.write(str(record.seq) + "\n")
            # logging.info("Entry %s has now id %d." % (el, new_id))
            # salvo il nuovo id HmtDB per usarlo nella fase dell'allineamento
            # new_hmtdb_file.write(str(new_id) + "\n")
    logging.info("Done.\n")


if __name__ == '__main__':
    if args.phase == 1:
        logging.basicConfig(filename="log_genome_search_ph1_{}.log".format(args.from_date),
                            format="%(levelname)s %(asctime)s: %(message)s",
                            datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)
        run_phase1()

    elif args.phase == 2:
        logging.basicConfig(filename="log_genome_search_ph2_{}.log".format(args.from_date),
                            format="%(levelname)s %(asctime)s: %(message)s",
                            datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)
        run_phase2()

    else:
        print "Please insert either 1 or 2 in the phase option."
        sys.exit(1)

