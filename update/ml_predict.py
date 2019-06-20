#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import argparse
import csv
import logging
import os
import sys
import time
import pandas as pd
from Bio import Entrez, Medline
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from app.static.dbdata import last_update
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression
# from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import VotingClassifier
from sklearn.pipeline import Pipeline

parser = argparse.ArgumentParser(description="""Modulo per eseguire il processing dei file 
    necessari per usare modello di Machine Learning per la predizione dei genomi healthy e 
    patient a partire dai dati delle pubblicazioni. """)
parser.add_argument("-u", "--update_training", dest="update_training", default=True,
                    action="store_false", help="""Update the main_input.csv file using data coming 
                    from the database to reproduce it. Default: True.""")
parser.add_argument("-i", "--main_train_file", dest="main_train_file", type=str,
                    default="trainer_ml/main_input.csv", help="""File csv da usare per il training 
                    del modello di ML. Composto dalle colonne [PubmedID, GenomeType].""")
parser.add_argument("-o", "--main_test_file", dest="main_test_file", type=str,
                    default="HmtDB_update_{}/PubMed_Results_{}.csv".format(last_update, last_update),
                    help="""File csv contenente i pubmed Id delle pubblicazioni di cui eseguire la 
                    predizione con il modello di ML. Se specificato, deve contenere almeno una 
                    colonna PubmedId, altrimenti di default viene usato il file Pubmed_Results 
                    prodotto dal primo step di genome_search.py.""")
parser.add_argument("-f", "--from_date", dest="from_date", default=last_update, type=str,
                    help="""Data a partire dalla quale ricercare i nuovi genomi pubblicati. 
                    Formato: mm_yyyy (i.e. 01_2018 per Gennaio 2018). Default: ultimo update del 
                    database.""")

# todo: inserire opzione per prendere automaticamente il file Pubmed_Results dalla cartella dell'ultimo aggiornamento

args = parser.parse_args()

Entrez.email = "hmtdb.update@uniba.it"
Entrez.tool = "HmtDB updating protocol"

logging.basicConfig(filename="log_ml_predict_{}.log".format(args.from_date),
                    format="%(levelname)s %(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)

main_input_dict = {}
main_test_list = []


def prepare_training_file():
    """
    Crea il file main_input.csv che viene poi utilizzato per il resto dell'analisi.
    :return:
    """

    # Reference table
    ref_df = pd.read_sql("SELECT referenceId, pubmedId, genomeId FROM Reference",
                         "sqlite:///{}".format(os.path.abspath("../hmtdb.db")))
    # Genome table
    gen_df = pd.read_sql("SELECT genomeId, genomeType FROM Genome",
                         "sqlite:///{}".format(os.path.abspath("../hmtdb.db")))

    train_df = ref_df.join(gen_df, on = "genomeId", lsuffix = "_ref", rsuffix = "_gen")
    train_df = train_df[["pubmedId", "genomeType"]]
    train_df.drop_duplicates("pubmedId", inplace = True)
    # Remove observations where the pubmedId is Family Tree DNA or is empty
    train_df.drop(train_df[train_df["pubmedId"] == "Family Tree DNA"].index, inplace = True)
    train_df.drop(train_df[train_df["pubmedId"].isna()].index, inplace = True)
    # Rename columns
    train_df.rename({"pubmedId": "PubmedId", "genomeType": "GenomeType"},
                    axis = "columns", inplace = True)

    train_df.to_csv("trainer_ml/main_input.csv", index = False)


def parse_medline(pm_file):
    """Effettua il parsing di un file Medline. Ritorna una lista di elementi, così composta:
    [abstract].
    Se i campi non vengono trovati nel file Medline, il default inserito è None (campo vuoto)."""

    record = Medline.read(pm_file)

    abstract = None

    try:
        abstract = record["AB"]
    except:
        abstract = ""

    data_list = [abstract]

    return data_list


def process_training(main_train_file):
    """Read main_train_csv (PubmedId, GenomeType) and produce model_train (Abstract, Patient)."""

    # main_input_csv: main_input.csv

    logging.info("Parsing trainer file {}...".format(main_train_file))
    with open(main_train_file, "r") as main_input:
        r = csv.reader(main_input)
        for el in r:
            if el[0] == "PubmedId":
                pass
            else:
                main_input_dict[el[0]] = 0 if el[1] == "N" else 1
    logging.info("Done.\n")

    logging.info("Creating file model_train.csv...")
    with open("trainer_ml/model_train.csv", "w") as model_input:
        w = csv.writer(model_input)
        w.writerow(["Abstract", "Patient"])

    for el in main_input_dict:
        time.sleep(1)

        for num_retry in range(5):
            try:
                logging.info("Retrieving PubmedId {}...".format(el))
                entry_pm = Entrez.efetch(db="pubmed", id=el, rettype="medline", retmode="text")
                logging.info("Done. ")
                break
            except Exception:
                logging.error("There was an error. Retry {} in 60 seconds...".format(num_retry))
                time.sleep(60)
        else:
            logging.error("Could not retrieve PubmedId {}. Stopping.".format(el))
            raise

        with open("trainer_ml/model_train.csv", "a") as model_input:
            w = csv.writer(model_input)
            res = parse_medline(entry_pm)
            res.append(main_input_dict[el])
            if len(res[0]) > 1:
                w.writerow(res)
    logging.info("Done.\n")


def process_testing(main_test_file):
    """Read main_test_csv (PubmedId) and produce model_test (PubmedId, Abstract)."""

    # main_test_csv: risultati dell'aggiornamento: Pubmed_Results_MM_YYYY.csv

    logging.info("Parsing tester file {}...".format(main_test_file))
    r = pd.read_csv(main_test_file)
    for el in set(r.PubmedId):
        if el != "Unpublished" and el != "Family Tree DNA":
            main_test_list.append(el)
    logging.info("Done.\n")

    logging.info("Creating file model_test.csv...")
    with open("trainer_ml/model_test.csv", "w") as model_test:
        w = csv.writer(model_test)
        w.writerow(["Abstract", "PubmedId"])

    for el in main_test_list:
        time.sleep(1)

        for num_retry in range(5):
            try:
                logging.info("Retrieving PubmedId {}...".format(el))
                entry_pm = Entrez.efetch(db="pubmed", id=el, rettype="medline", retmode="text")
                logging.info("Done. ")
                break
            except Exception:
                logging.error("There was an error. Retry {} in 60 seconds...".format(num_retry))
                time.sleep(60)
        else:
            logging.error("Could not retrieve PubmedId {}. Stopping.".format(el))
            raise

        with open("trainer_ml/model_test.csv", "a") as model_test:
            w = csv.writer(model_test)
            res = parse_medline(entry_pm)
            res.append(el)
            w.writerow(res)
    logging.info("Done.\n")


def ml_predictor(trainer, tester):
    """Perform the machine learning approach to predict healthy/patient genomes."""

    # expl: training set of publications
    train_df = pd.read_csv(trainer)
    # expl: publications to predict
    test_df = pd.read_csv(tester)

    X_train = train_df.Abstract
    y_train = train_df.Patient
    X_test = test_df.Abstract

    # Models
    clf1 = Pipeline([("vect1", CountVectorizer(stop_words="english", ngram_range=(1, 1))),
                     ("LR", LogisticRegression())])
    clf2 = Pipeline([("vect2", CountVectorizer(stop_words="english", ngram_range=(1, 1))),
                     ("NB", MultinomialNB())])
    svm = LinearSVC()
    clf3 = Pipeline([("vect3", CountVectorizer(stop_words="english", ngram_range=(1, 1))),
                     ("SVM_cal", CalibratedClassifierCV(svm))])
    clf4 = Pipeline([("vect4", CountVectorizer(stop_words="english", ngram_range=(1, 1))),
                     ("ET", ExtraTreesClassifier())])
    clf5 = Pipeline([("vect5", CountVectorizer(stop_words="english", ngram_range=(1, 1))),
                     ("RF", RandomForestClassifier())])

    # Ensemble voting
    logging.info("Learning predictive model...")
    ens = VotingClassifier(estimators=[("LR", clf1), ("NB", clf2), ("SVM", clf3), ("ET", clf4),
                                       ("RF", clf5)],
                           voting="soft")
    ens = ens.fit(X_train, y_train)
    logging.info("Done.")

    # Prediction
    logging.info("Predicting classification and probabilities...")
    y_pred = ens.predict(X_test)
    y_prob = ens.predict_proba(X_test)
    logging.info("Done.")

    pred_list = []
    for n, el in enumerate(test_df.PubmedId):
        pred_list.append((el, y_pred[n], y_prob[n][int(y_pred[n])]))

    # Save results
    logging.info("Saving prediction results in trainer_ml/model_output.csv...")
    pred_df = pd.DataFrame.from_records(pred_list, columns=["PubmedId", "PredClass", "PredProba"])
    pred_df.reset_index(drop=True).to_csv("trainer_ml/model_output.csv")
    logging.info("Done.\n")


def compile_pubmed_results(main_test_file):
    """Compile the last column of main_test_file (Pubmed_Results csv) based on prediction results,
    reporting for each genome whether it is a healthy or patient genome."""

    # todo: non serve per ora finché tanto ci sono anche i genomi Unpublished

    logging.info("Parsing publications counter from {}...".format(main_test_file))
    in_dict = {}
    in_df = pd.read_csv(main_test_file)
    pred_df = pd.read_csv("trainer_ml/model_output.csv", index_col=0)
    pred_df.index += 1  # index fa le veci di Counter negli altri csv
    for el in in_df.itertuples():
        print el.Counter
        in_dict[el.Counter] = [el.PubmedId, pred_df.iloc[el.Counter].PredClass]
    logging.info("Done.")

    print(in_dict)

    # logging.info("Loading ML predictions...")
    #
    # for counter in in_dict:
    #     in_dict[counter].append(pred_df[counter].PredClass)
    # logging.info("Done.")

    logging.info("Compiling Gen_Healthy column in GenBank Results csv...")
    gb_df = pd.read_csv("HmtDB_update_{}/GenBank_Results_{}.csv".format(last_update, last_update),
                        index_col=0)
    for counter in in_dict:
        gb_df[counter].Gen_Healthy = "N" if in_dict[counter][1] == 0 else "P"
    logging.info("Done.")

    logging.info("Overwriting GenBank Results file...")
    gb_df.to_csv("HmtDB_update_{}/GenBank_Results_{}.csv".format(last_update, last_update))
    logging.info("Done.\n")


if __name__ == '__main__':
    if args.update_training:
        prepare_training_file()
    process_training(args.main_train_file)
    process_testing(args.main_test_file)
    ml_predictor("trainer_ml/model_train.csv", "trainer_ml/model_test.csv")
    # compile_pubmed_results(args.main_test_file)
