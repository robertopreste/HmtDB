#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import json
import requests
import operator
from .models import Country, Sources, Disease, Reference, Genome, Stats, NtVariability, GenAlignment


def retrieveHmtVar(pos, alt):
    """Retrieve information from HmtVar so that the Variant Card
    can be shown.

    :param int pos: mitochondrial input position

    :param str alt: alternative allele

    :return: str
    """
    url = "https://www.hmtvar.uniba.it/rdconnect?variant_start={}&variant_alternateBases={}".format(
        pos, alt)

    resp = requests.get(url)
    json_res = resp.json()

    if resp.ok and json_res["success"] == "true":
        hmtvar_id = json_res["url"].split("/")[-1]
    else:
        hmtvar_id = "0"

    return hmtvar_id


def getCountries(cont):
    """Return the list of countryId and countryName for the given
    continent code.

    :param str cont: input continent code

    :return: List[Tuple[int,str]]
    """
    res = set()
    countries = Country.query.filter_by(continentCode=cont).all()
    for el in countries:
        res.add((el.countryId, el.countryName))
    res_list = list(res)
    res_list.sort(key=operator.itemgetter(1))
    return res_list


def getSources():
    """Return the list of sourceId and sourceName for each entry
    in Sources.

    :return: List[Tuple[int,str]]
    """
    res = set()
    sources = Sources.query.all()
    for el in sources:
        res.add((el.sourceId, el.sourceName))
    res_list = list(res)
    res_list.sort(key=operator.itemgetter(1))
    return res_list


def getDiseases():
    """Return the list of diseaseId and diseaseName for each entry
    in Disease.

    :return: List[Tuple[int,str]]
    """
    res = set()
    diseases = Disease.query.all()
    for el in diseases:
        res.add((el.diseaseId, el.diseaseName))
    res_list = list(res)
    res_list.sort(key=operator.itemgetter(1))
    return res_list


def getJournals():
    """Return the list of papers for each entry in Reference.

    :return: List[Tuple[str,str]]
    """
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


def get_alignments(id_list, with_ref=True):
    """Return a list of tuples with (genomeId, alignment) from query results.

    :param id_list: list of genomeIds
    :param with_ref: include reference genome as first element (default: True)
    :return: [(gnomeId, alignment), ...]
    """
    if with_ref:
        id_list = [0] + id_list
    q = GenAlignment.query.filter(GenAlignment.genomeId.in_(id_list)).all()
    algs = [(el.genomeId, el.alignment) for el in q]

    return algs


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

