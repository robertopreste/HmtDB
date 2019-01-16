#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import json
import requests
from .models import Country, Sources, Disease, Reference, Genome, Stats, NtVariability


def retrieveHmtVar(pos, alt):
    """Retrieve information from HmtVar so that the Variant Card can be shown."""

    url = "https://www.hmtvar.uniba.it/rdconnect?variant_start={}&variant_alternateBases={}".format(
        pos, alt)

    mainApi = requests.get(url)
    json_res = json.loads(mainApi.text)

    if isinstance(json_res, dict):
        if json_res["success"] == "false":
            resp = 0
        else:
            hmtvar_id = json_res["url"].split("/")[-1]
            resp = hmtvar_id
    else:
        resp = 0

    return resp


def getCountries(cont):
    lista = []
    for el in Country.query.filter_by(continentCode=cont).all():
        if (el.countryId, el.countryName) not in lista:
            lista.append((el.countryId, el.countryName))
    return lista


def getSources():
    lista = set()
    for el in Sources.query.all():
        lista.add((el.sourceId, el.sourceName))
    return list(lista)


def getDiseases():
    lista = set()
    for el in Disease.query.all():
        lista.add((el.diseaseId, el.diseaseName))
    return list(lista)


def getJournals():
    lista = set()
    for el in Reference.query.all():
        lista.add((el.paper, el.paper))
    return list(lista)


def getHaplogroups(letter):
    lista = []
    for el in Genome.query.order_by(Genome.haplogroupHmdb).filter(
            Genome.haplogroupHmdb.like(letter + "%")).all():
        if (el.haplogroupHmdb, el.haplogroupHmdb) not in lista:
            lista.append((el.haplogroupHmdb, el.haplogroupHmdb))
    return lista


def getUserHaplogroups(letter):
    lista = []
    for el in Genome.query.order_by(Genome.haplogroupUser).filter(
            Genome.haplogroupUser.like(letter + "%")).all():
        if (el.haplogroupUser, el.haplogroupUser) not in lista:
            lista.append((el.haplogroupUser, el.haplogroupUser))
    return lista


def getStats():
    worldDict = {
        "Africa": {
            "N": [],
            "P": [],
            "T": 0
        },
        "America": {
            "N": [],
            "P": [],
            "T": 0
        },
        "Asia": {
            "N": [],
            "P": [],
            "T": 0
        },
        "Europe": {
            "N": [],
            "P": [],
            "T": 0
        },
        "Oceania": {
            "N": [],
            "P": [],
            "T": 0
        },
        "Undefined Continent": {
            "N": [],
            "P": [],
            "T": 0
        }
    }

    all_dict = {
        "N": [0, 0, 0],
        "P": [0, 0, 0],
        "T": 0
    }

    for continent in worldDict:
        norm_compl_q = Stats.query.filter(Stats.continentName == continent) \
            .filter(Stats.genomeType == "N", Stats.completeGenome == "Y").first()
        norm_coding_q = Stats.query.filter(Stats.continentName == continent) \
            .filter(Stats.genomeType == "N", Stats.completeGenome == "N").first()
        pat_compl_q = Stats.query.filter(Stats.continentName == continent) \
            .filter(Stats.genomeType == "P", Stats.completeGenome == "Y").first()
        pat_coding_q = Stats.query.filter(Stats.continentName == continent) \
            .filter(Stats.genomeType == "P", Stats.completeGenome == "N").first()

        norm_compl = int(norm_compl_q.total) if norm_compl_q else 0
        norm_coding = int(norm_coding_q.total) if norm_coding_q else 0
        pat_compl = int(pat_compl_q.total) if pat_compl_q else 0
        pat_coding = int(pat_coding_q.total) if pat_coding_q else 0
        norm_tot = norm_compl + norm_coding
        pat_tot = pat_compl + pat_coding

        worldDict[continent]["N"] = [norm_tot, norm_compl, norm_coding]
        worldDict[continent]["P"] = [pat_tot, pat_compl, pat_coding]
        worldDict[continent]["T"] = norm_tot + pat_tot

        all_dict["N"][0] += norm_tot
        all_dict["N"][1] += norm_compl
        all_dict["N"][2] += norm_coding
        all_dict["P"][0] += pat_tot
        all_dict["P"][1] += pat_compl
        all_dict["P"][2] += pat_coding
        all_dict["T"] += (norm_tot + pat_tot)

    worldDict["All Continents"] = all_dict

    return worldDict


def getVariants():
    africaQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                         NtVariability.var_af > 0,
                                         NtVariability.genomeType == "N").count()
    americaQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                          NtVariability.var_am > 0,
                                          NtVariability.genomeType == "N").count()
    asiaQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                       NtVariability.var_as > 0,
                                       NtVariability.genomeType == "N").count()
    europeQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                         NtVariability.var_eu > 0,
                                         NtVariability.genomeType == "N").count()
    oceaniaQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                          NtVariability.var_oc > 0,
                                          NtVariability.genomeType == "N").count()
    totalQ = NtVariability.query.filter(NtVariability.insertionPosition == 0,
                                        NtVariability.var_tot > 0,
                                        NtVariability.genomeType == "N").count()

    varDict = {
        "Africa": africaQ,
        "America": americaQ,
        "Asia": asiaQ,
        "Europe": europeQ,
        "Oceania": oceaniaQ,
        "All Continents": totalQ
    }

    return varDict


def populateHaploScript():
    stringa = """
function populateHaplogroups(s1, s2) {
    // Inserts haplogroups values in the HmtDB-assigned list based on the selected macrohaplogroup
    s1 = document.getElementById(s1);
    s2 = document.getElementById(s2);
    var optionArray;

    s2.innerHTML = "--All Haplogroups--";
    """

    stringa += """
    if (s1.value == "A") {
        optionArray = ["|--All Haplogroups--", """
    haplos = getHaplogroups("A")
    for tup in haplos:
        stringa += """
        "%s|%s", """ % (tup[0], tup[1])
    stringa += """]; """

    for lett in ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q",
                 "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]:
        stringa += """
    } else if (s1.value == "%s") {
        optionArray = ["|--All Haplogroups--", """ % (lett,)
        haplos = getHaplogroups(lett)
        for tup in haplos:
            stringa += """
            "{t[0]}|{t[1]}", """.format(t=tup)
        stringa += """]; """

    stringa += """
    }

    for (var i = s2.options.length - 1; i >= 0; i--) {
        s2.remove(i); 
    } 

    for (var option in optionArray) {
        console.log(option);
        console.log(optionArray[option]);
        var pair = optionArray[option].split("|"); 
        var newOption = document.createElement("option"); 
        newOption.value = pair[0]; 
        newOption.innerHTML = pair[1]; 

        s2.options.add(newOption); 
    }
}
    """

    return stringa


def populateUserHaploScript():
    stringa = """
function populateUserHaplogroups(s1, s2) {
    // Inserts haplogroups values in the User-assigned list based on the selected macrohaplogroup
    s1 = document.getElementById(s1);
    s2 = document.getElementById(s2);
    var optionArray;

    s2.innerHTML = "--All Haplogroups--";
    """

    stringa += """
    if (s1.value == "A") {
        optionArray = ["|--All Haplogroups--", """
    haplos = getUserHaplogroups("A")
    for tup in haplos:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)
    stringa += """]; """

    for lett in ["B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q",
                 "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]:
        stringa += """
    } else if (s1.value == "%s") {
        optionArray = ["|--All Haplogroups--", """ % lett
        haplos = getUserHaplogroups(lett)
        for tup in haplos:
            stringa += """
            "{t[0]}|{t[1]}", """.format(t=tup)
        stringa += """]; """

    stringa += """
    }

    for (var i = s2.options.length - 1; i >= 0; i--) {
        s2.remove(i); 
    } 

    for (var option in optionArray) {
        console.log(option);
        console.log(optionArray[option]);
        var pair = optionArray[option].split("|"); 
        var newOption = document.createElement("option"); 
        newOption.value = pair[0]; 
        newOption.innerHTML = pair[1]; 

        s2.options.add(newOption); 
    }
}
    """

    return stringa


def populateCountriesScript():
    stringa = """
function populateCountries(s1, s2) {
    // Inserts countries values in the country list based on the selected continent
    s1 = document.getElementById(s1);
    s2 = document.getElementById(s2);
    var optionArray;

    s2.innerHTML = "--All Countries--";
    """

    # africa
    stringa += """
    if (s1.value == "AF") {
        optionArray = ["|--All Countries--", """
    countrAF = getCountries("AF")
    for tup in countrAF:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)

    # america
    stringa += """];
    } else if (s1.value == "AM") {
        optionArray = ["|--All Countries--", """
    countrAM = getCountries("AM")
    for tup in countrAM:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)

    # asia
    stringa += """];
    } else if (s1.value == "AS") {
        optionArray = ["|--All Countries--", """
    countrAF = getCountries("AS")
    for tup in countrAF:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)

    # europe
    stringa += """];
    } else if (s1.value == "EU") {
        optionArray = ["|--All Countries--", """
    countrAF = getCountries("EU")
    for tup in countrAF:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)

    # oceania
    stringa += """];
    } else if (s1.value == "OC") {
        optionArray = ["|--All Countries--", """
    countrAF = getCountries("OC")
    for tup in countrAF:
        stringa += """
        "{t[0]}|{t[1]}", """.format(t=tup)

    stringa += """]; 
    }

    for (var option in optionArray) {
        var pair = optionArray[option].split("|");
        var newOption = document.createElement("option");

        newOption.value = pair[0];
        newOption.innerHTML = pair[1];
        s2.options.add(newOption);
    }
}
"""

    return stringa


def otherFunctions():
    stringa = """
function openCenteredPopup(target, w, h) {
    // Opens centered popups for elements like Statistics, People, ecc.
    var t = Math.floor((screen.height-h)/2);
    var l = Math.floor((screen.width-w)/2);
    var stile = "top=" + t + ", left=" + l + ", width=" + w + ", height=" + h + ", status=no, menubar=no, toolbar=no, scrollbar=no";
    window.open(target, "", stile);
}

"""
    return stringa

