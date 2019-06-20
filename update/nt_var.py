#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import sys
from os import path
import os
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import logging
import argparse
import pandas as pd
# pd.set_option("display.precision", 5)
from app import db
from app.site.models import NtVariability
from app.static.dbdata import last_update


parser = argparse.ArgumentParser(description="""Modulo per effettuare il parsing dei sitevar 
    generati dallo script calc_sitevar.py e caricare i nuovi dati sul db. Crea anche i file csv 
    con gli NtVar cumulativi che possono essere usati per aggiornare HmtVar.""")
parser.add_argument("-M", "--complete", action="store_true", dest="on_complete", 
                    help="Perform calculations only on complete genomes (default: False).", 
                    default=False)
parser.add_argument("-l", "--local", action="store_true", dest="local",
                    help="Do not work on db, but rather create local files (default: False).",
                    default=False)
parser.add_argument("-c", "--calculate", action="store_true", dest="calculate",
                    help="Calculate sitevars before creating the cumulative ntvar file (default: "
                         "False).",
                    default=False)

args = parser.parse_args()


logging.basicConfig(filename="log_ntvar_{}.log".format(last_update),
                    format="%(levelname)s\t%(asctime)s: %(message)s",
                    datefmt="%d/%m/%Y %I:%M:%S", filemode="w", level=logging.DEBUG)


def delete_table():
    logging.info("Resetting the NtVariability table...")
    db.session.query(NtVariability).delete()
    db.session.commit()
    logging.info("Done.\n")

    return


def load_var(ntvar):
    """
    Load the variability values into the NtVariability table.
    :param ntvar: pandas dataframe with the NtVar (healthy or patient)
    """
    to_add = []
    logging.info("Collecting values to update...")
    for el in ntvar.itertuples():
        to_add.append(NtVariability(nucleotidePosition=el.nucleotidePosition,
                                    insertionPosition=el.insertionPosition,
                                    var_tot=el.var_tot,
                                    var_eu=el.var_eu,
                                    var_am=el.var_am,
                                    var_af=el.var_af,
                                    var_as=el.var_as,
                                    var_oc=el.var_oc,
                                    compVar_tot=el.compVar_tot,
                                    compVar_eu=el.compVar_eu,
                                    compVar_am=el.compVar_am,
                                    compVar_af=el.compVar_af,
                                    compVar_as=el.compVar_as,
                                    compVar_oc=el.compVar_oc,
                                    genomeType=el.genomeType))
    logging.info("Updating the NtVariability table...")
    db.session.add_all(to_add)
    db.session.commit()
    logging.info("Done.")

    return


def parse_tables():
    logging.info("Loading all the sitevars...")
    
    if args.on_complete: 
        vars_df = ["sitevar_healthy_complete", "sitevar_patient_complete", 
                   "sitevar_healthy_complete_AF", "sitevar_healthy_complete_AM", 
                   "sitevar_healthy_complete_AS", "sitevar_healthy_complete_EU", 
                   "sitevar_healthy_complete_OC", "sitevar_patient_complete_AF", 
                   "sitevar_patient_complete_AM", "sitevar_patient_complete_AS", 
                   "sitevar_patient_complete_EU", "sitevar_patient_complete_OC"]
    else: 
        vars_df = ["sitevar_healthy", "sitevar_patient", 
                   "sitevar_healthy_AF", "sitevar_healthy_AM", 
                   "sitevar_healthy_AS", "sitevar_healthy_EU", 
                   "sitevar_healthy_OC", "sitevar_patient_AF", 
                   "sitevar_patient_AM", "sitevar_patient_AS", 
                   "sitevar_patient_EU", "sitevar_patient_OC"]
    
    # Load the sitevars
    sv_h = pd.read_csv("sitevars/{}".format(vars_df[0]), header=0,
                       names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                       dtype={"site": str})
    sv_h = sv_h[sv_h["site"].astype(float) < 16570]

    sv_p = pd.read_csv("sitevars/{}".format(vars_df[1]), header=0,
                       names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                       dtype={"site": str})
    sv_p = sv_p[sv_p["site"].astype(float) < 16570]

    sv_h_af = pd.read_csv("sitevars/{}".format(vars_df[2]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_h_af = sv_h_af[sv_h_af["site"].astype(float) < 16570]

    sv_h_am = pd.read_csv("sitevars/{}".format(vars_df[3]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_h_am = sv_h_am[sv_h_am["site"].astype(float) < 16570]

    sv_h_as = pd.read_csv("sitevars/{}".format(vars_df[4]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_h_as = sv_h_as[sv_h_as["site"].astype(float) < 16570]

    sv_h_eu = pd.read_csv("sitevars/{}".format(vars_df[5]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_h_eu = sv_h_eu[sv_h_eu["site"].astype(float) < 16570]

    sv_h_oc = pd.read_csv("sitevars/{}".format(vars_df[6]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_h_oc = sv_h_oc[sv_h_oc["site"].astype(float) < 16570]

    sv_p_af = pd.read_csv("sitevars/{}".format(vars_df[7]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_p_af = sv_p_af[sv_p_af["site"].astype(float) < 16570]

    sv_p_am = pd.read_csv("sitevars/{}".format(vars_df[8]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_p_am = sv_p_am[sv_p_am["site"].astype(float) < 16570]

    sv_p_as = pd.read_csv("sitevars/{}".format(vars_df[9]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_p_as = sv_p_as[sv_p_as["site"].astype(float) < 16570]

    sv_p_eu = pd.read_csv("sitevars/{}".format(vars_df[10]), header=0,
                          names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                          dtype={"site": str})
    sv_p_eu = sv_p_eu[sv_p_eu["site"].astype(float) < 16570]

    try:  # expl: difficile che questo file venga creato
        sv_p_oc = pd.read_csv("sitevars/{}".format(vars_df[11]), header=0,
                              names=["site", "variability", "nucl_type", "A", "C", "G", "T", "Gap", "Others"],
                              dtype={"site": str})
        sv_p_oc = sv_p_oc[sv_p_oc["site"].astype(float) < 16570]

    except IOError:
        sv_p_oc = sv_p.copy()
        sv_p_oc["variability"] = 0.0
        sv_p_oc["A"] = 0.0
        sv_p_oc["C"] = 0.0
        sv_p_oc["G"] = 0.0
        sv_p_oc["T"] = 0.0
        sv_p_oc["Gap"] = 0.0
        sv_p_oc["Others"] = 0.0

    logging.info("Done. ")

    dfs = [sv_h, sv_h_af, sv_h_am, sv_h_as, sv_h_eu, sv_h_oc,
           sv_p, sv_p_af, sv_p_am, sv_p_as, sv_p_eu, sv_p_oc]

    # Data manipulation
    for df in dfs:
        # df.site = df.site.astype(str)
        # df["variation"] = df.site + "_" + df.nucl_type
        # df.drop_duplicates("variation", inplace=True)
        df["ntPos"], df["insPos"] = df["site"].str.split(".").str
        # df.set_index("variation", inplace=True)
        df["ntPos"].fillna(0, inplace=True)
        df["insPos"].fillna(0, inplace=True)
        df["ntPos"] = df["ntPos"].astype(int)
        df["insPos"] = df["insPos"].astype(int)
        df.set_index("site", inplace=True)

    healthy_list = []
    patient_list = []

    logging.info("Creating the complete healthy NtVar...")
    # Healthy complete NtVar
    for el in sv_h.index:
        new_dict = {"nucleotidePosition": sv_h.loc[el]["ntPos"],
                    "insertionPosition": sv_h.loc[el]["insPos"],
                    "var_tot": sv_h.loc[el]["variability"],
                    "var_eu": sv_h_eu["variability"].get(el, 0.0),
                    "var_am": sv_h_am["variability"].get(el, 0.0),
                    "var_af": sv_h_af["variability"].get(el, 0.0),
                    "var_as": sv_h_as["variability"].get(el, 0.0),
                    "var_oc": sv_h_oc["variability"].get(el, 0.0),
                    "compVar_tot": ";".join([str(sv_h.loc[el]["A"]), str(sv_h.loc[el]["C"]),
                                             str(sv_h.loc[el]["G"]), str(sv_h.loc[el]["T"]),
                                             str(sv_h.loc[el]["Gap"] + sv_h.loc[el]["Others"])]),
                    "compVar_eu": ";".join([str(sv_h_eu["A"].get(el, 0.0)), str(sv_h_eu["C"].get(el, 0.0)),
                                            str(sv_h_eu["G"].get(el, 0.0)), str(sv_h_eu["T"].get(el, 0.0)),
                                            str(sv_h_eu["Gap"].get(el, 100.0) + sv_h_eu["Others"].get(el, 0.0))]),
                    "compVar_am": ";".join([str(sv_h_am["A"].get(el, 0.0)), str(sv_h_am["C"].get(el, 0.0)),
                                            str(sv_h_am["G"].get(el, 0.0)), str(sv_h_am["T"].get(el, 0.0)),
                                            str(sv_h_am["Gap"].get(el, 100.0) + sv_h_am["Others"].get(el, 0.0))]),
                    "compVar_af": ";".join([str(sv_h_af["A"].get(el, 0.0)), str(sv_h_af["C"].get(el, 0.0)),
                                            str(sv_h_af["G"].get(el, 0.0)), str(sv_h_af["T"].get(el, 0.0)),
                                            str(sv_h_af["Gap"].get(el, 100.0) + sv_h_af["Others"].get(el, 0.0))]),
                    "compVar_as": ";".join([str(sv_h_as["A"].get(el, 0.0)), str(sv_h_as["C"].get(el, 0.0)),
                                            str(sv_h_as["G"].get(el, 0.0)), str(sv_h_as["T"].get(el, 0.0)),
                                            str(sv_h_as["Gap"].get(el, 100.0) + sv_h_as["Others"].get(el, 0.0))]),
                    "compVar_oc": ";".join([str(sv_h_oc["A"].get(el, 0.0)), str(sv_h_oc["C"].get(el, 0.0)),
                                            str(sv_h_oc["G"].get(el, 0.0)), str(sv_h_oc["T"].get(el, 0.0)),
                                            str(sv_h_oc["Gap"].get(el, 100.0) + sv_h_oc["Others"].get(el, 0.0))]),
                    "genomeType": "N"}
        healthy_list.append(new_dict)
    ntvar_h = pd.DataFrame(healthy_list)
    logging.info("Done.")

    logging.info("Creating the complete patient NtVar...")
    # Patient complete NtVar
    for el in sv_p.index:
        new_dict = {"nucleotidePosition": sv_p.loc[el]["ntPos"],
                    "insertionPosition": sv_p.loc[el]["insPos"],
                    "var_tot": sv_p.loc[el]["variability"],
                    "var_eu": sv_p_eu["variability"].get(el, 0.0),
                    "var_am": sv_p_am["variability"].get(el, 0.0),
                    "var_af": sv_p_af["variability"].get(el, 0.0),
                    "var_as": sv_p_as["variability"].get(el, 0.0),
                    "var_oc": sv_p_oc["variability"].get(el, 0.0),
                    "compVar_tot": ";".join([str(sv_p.loc[el]["A"]), str(sv_p.loc[el]["C"]),
                                             str(sv_p.loc[el]["G"]), str(sv_p.loc[el]["T"]),
                                             str(sv_p.loc[el]["Gap"] + sv_p.loc[el]["Others"])]),
                    "compVar_eu": ";".join([str(sv_p_eu["A"].get(el, 0.0)), str(sv_p_eu["C"].get(el, 0.0)),
                                            str(sv_p_eu["G"].get(el, 0.0)), str(sv_p_eu["T"].get(el, 0.0)),
                                            str(sv_p_eu["Gap"].get(el, 100.0) + sv_p_eu["Others"].get(el, 0.0))]),
                    "compVar_am": ";".join([str(sv_p_am["A"].get(el, 0.0)), str(sv_p_am["C"].get(el, 0.0)),
                                            str(sv_p_am["G"].get(el, 0.0)), str(sv_p_am["T"].get(el, 0.0)),
                                            str(sv_p_am["Gap"].get(el, 100.0) + sv_p_am["Others"].get(el, 0.0))]),
                    "compVar_af": ";".join([str(sv_p_af["A"].get(el, 0.0)), str(sv_p_af["C"].get(el, 0.0)),
                                            str(sv_p_af["G"].get(el, 0.0)), str(sv_p_af["T"].get(el, 0.0)),
                                            str(sv_p_af["Gap"].get(el, 100.0) + sv_p_af["Others"].get(el, 0.0))]),
                    "compVar_as": ";".join([str(sv_p_as["A"].get(el, 0.0)), str(sv_p_as["C"].get(el, 0.0)),
                                            str(sv_p_as["G"].get(el, 0.0)), str(sv_p_as["T"].get(el, 0.0)),
                                            str(sv_p_as["Gap"].get(el, 100.0) + sv_p_as["Others"].get(el, 0.0))]),
                    "compVar_oc": ";".join([str(sv_p_oc["A"].get(el, 0.0)), str(sv_p_oc["C"].get(el, 0.0)),
                                            str(sv_p_oc["G"].get(el, 0.0)), str(sv_p_oc["T"].get(el, 0.0)),
                                            str(sv_p_oc["Gap"].get(el, 100.0) + sv_p_oc["Others"].get(el, 0.0))]),
                    "genomeType": "P"}
        patient_list.append(new_dict)
    ntvar_p = pd.DataFrame(patient_list)
    logging.info("Done. ")

    if args.on_complete:
        ntvar_h.to_csv("sitevars/ntvar_healthy_complete.csv", index=None)
        ntvar_p.to_csv("sitevars/ntvar_patient_complete.csv", index=None)
    else:
        ntvar_h.to_csv("sitevars/ntvar_healthy.csv", index=None)
        ntvar_p.to_csv("sitevars/ntvar_patient.csv", index=None)

    if args.local:
        logging.info("New ntvar csv files saved.")
    else:
        logging.info("Loading healthy genomes...")
        load_var(ntvar_h)
        logging.info("Done.")
        logging.info("Loading patient genomes...")
        load_var(ntvar_p)
        logging.info("Done.")

    return


def launch_sitevar(genome_type, continent=None):
    """
    Launch the calc_sitevar.py script in order to obtain the sitevar files.
    :param genome_type: ["N", "P"]
    :param continent: ["AF", "AM", "AS", "EU", "OC"]
    """
    logging.info("Calculating SiteVar on {} genomes...".format("healthy" if genome_type == "N" else "patient"))
    command = "python calc_sitevar.py "
    if genome_type == "N":
        command += "-H "
        if continent:
            logging.info("Continent {}".format(continent))
            command += " -C {}".format(continent)
    elif genome_type == "P":
        if continent:
            logging.info("Continent {}".format(continent))
            command += " -C {}".format(continent)

    if args.on_complete:
        command += " -M"

    os.system(command)

    return


if __name__ == '__main__':
    if args.calculate:
        launch_sitevar("N")
        launch_sitevar("P")

        for cont in ["AF", "AM", "AS", "EU", "OC"]:
            launch_sitevar("N", cont)
            launch_sitevar("P", cont)

    if args.local:
        parse_tables()
    else:
        delete_table()
        parse_tables()



