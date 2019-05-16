#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import requests
import json
import os
from flask import Blueprint, render_template, flash, redirect, session, \
    url_for, request, g, jsonify, send_file, after_this_request
from .scripts import retrieveHmtVar, get_alignments

www = Blueprint("site", __name__)

from sqlalchemy import or_, and_
from app import app, db
from config import ADMINS
# from flask_login import current_user, login_user, logout_user, login_required
from .forms import QueryForm
from .models import Genome, Country, IndividualsData, GenomeSnp, Insertion, \
    Sources, Methods, Deletion, Reference, EthnicGroups, Disease, NtVariability, \
    AaVariability
from .query import queryLocus, queryNtVar_N, queryNtVar_P, queryMitomapDna, \
    queryMitomapAa, queryAaVar_N, queryAaVar_P, queryDisease, queryDeletion, \
    queryInsertion, getAltCodon, aa_dict, getAa, queryMitomapDnaDiseases
from app.static import dbdata


# RD-Connect common API
@www.route("/rdconnect", methods=["GET"])
def rdconnect():
    """Return information about one specific variant in the format
    recommended by RD-Connect.
    The request gets forwarded to HmtVar, then the response is parsed to
    retrieve the variant Id, which is then applied to the variantCard
    function below."""

    gene_id = request.args.get("gene_id")
    variant_assembly = request.args.get("variant_assembly")
    variant_chromosome = request.args.get("variant_chromosome")
    variant_start = request.args.get("variant_start")
    variant_end = request.args.get("variant_end")
    variant_refBases = request.args.get("variant_referenceBases")
    variant_altBases = request.args.get("variant_alternateBases")
    variant_name = request.args.get("variant_name")

    apiString = "https://www.hmtvar.uniba.it/rdconnect?"

    if gene_id:
        apiString += "gene_id={}&".format(gene_id.upper())

    if variant_start:
        apiString += "variant_start={}&".format(variant_start)

    if variant_end:
        apiString += "variant_end={}&".format(variant_end)

    if variant_refBases:
        apiString += "variant_referenceBases={}&".format(variant_refBases.upper())

    if variant_altBases:
        apiString += "variant_alternateBases={}&".format(variant_altBases.upper())

    if variant_name:
        pass

    if variant_assembly:
        if variant_assembly != "hg19" and variant_assembly != "GRCh37":
            return jsonify({"success": "false"})

    if variant_chromosome:
        if variant_chromosome != "M" and variant_chromosome != "MT":
            return jsonify({"success": "false"})

    mainApi = requests.get(apiString)
    json_res = json.loads(mainApi.text)

    if isinstance(json_res, list):
        resp = []
        for el in json_res:
            hmtvar_id = el["url"].split("/")[-1]
            resp.append({"url": url_for("site.variantCard", idVar=hmtvar_id,
                                        _external=True, _scheme="https"),
                         "success": "true"})
    elif isinstance(json_res, dict):
        if json_res["success"] == "false":
            resp = {"success": "false"}
        else:
            hmtvar_id = json_res["url"].split("/")[-1]
            resp = {"url": url_for("site.variantCard", idVar=hmtvar_id,
                                   _external=True, _scheme="https"),
                    "success": "true"}
    else:
        resp = {"success": "false"}

    return jsonify(resp)


# Variant Card to show data required by RD-Connect
@www.route("/variantCard/<int:idVar>", methods=["GET"])
def variantCard(idVar):

    hmtvar_res = requests.get("https://www.hmtvar.uniba.it/api/main/{}".format(idVar))
    json_res = json.loads(hmtvar_res.text)

    nt_start = json_res["nt_start"]
    ref_nt = json_res["ref_rCRS"]
    alt_nt = json_res["alt"]
    locus = json_res["locus"]
    dis_score = json_res["disease_score"]
    pathogen = json_res["pathogenicity"]
    nt_var_h = json_res["Variab"]["nt_var"]
    nt_var_p = json_res["Variab"]["nt_var_patients"]
    aa_var_h = json_res["Variab"]["aa_var"]
    aa_var_p = json_res["Variab"]["aa_var_patients"]

    return render_template("variantCard.html",
                           title="Variant Card",
                           varId=idVar,
                           position=nt_start,
                           from_nt=ref_nt,
                           to_nt=alt_nt,
                           locus=locus,
                           dis_score=dis_score,
                           pathogen=pathogen,
                           nt_var_h=nt_var_h,
                           nt_var_p=nt_var_p,
                           aa_var_h=aa_var_h,
                           aa_var_p=aa_var_p)


# Home Page
@www.route("/index", methods=["GET"])
@www.route("/home", methods=["GET"])
@www.route("/", methods=["GET"])
def index():
    return render_template("index.html",
                           title="Human Mitochondrial Database",
                           totalGens=dbdata.worldDict["All Continents"]["T"],
                           latest_update=dbdata.latest_update)


# Old Home Page link
@www.route("/hmdb", methods=["GET"])
@www.route("/hmdb/", methods=["GET"])
def index2():
    flash("Please note that the HmtDB web address has changed!")
    flash("The correct address is https://www.hmtdb.uniba.it, without the trailing /hmdb/.")
    return render_template("index.html",
                           title="Human Mitochondrial Database",
                           totalGens=dbdata.worldDict["All Continents"]["T"],
                           latest_update=dbdata.latest_update)


# About Page
@www.route("/hmdb/about", methods=["GET"])
@www.route("/about", methods=["GET"])
def about():
    return render_template("about.html",
                           title="About HmtDB")


# API Page
@www.route("/hmdb/apis", methods=["GET"])
@www.route("/apis", methods=["GET"])
def apis():
    return render_template("apis.html",
                           title="HmtDB API")


# DB Functions Page
@www.route("/hmdb/dbfunctions", methods=["GET"])
@www.route("/dbfunctions", methods=["GET"])
def dbfunctions():
    return render_template("dbfunctions.html",
                           title="Database Functions")


# DB Citations Page
@www.route("/hmdb/dbcitations", methods=["GET"])
@www.route("/dbcitations", methods=["GET"])
def dbcitations():
    return render_template("dbcitations.html",
                           title="Citing HmtDB")


# What's New Page
@www.route("/hmdb/whatsnew", methods=["GET"])
@www.route("/whatsnew", methods=["GET"])
def whatsnew():
    return render_template("whatsnew.html",
                           title="What's New",
                           latest_update=dbdata.latest_update)


# People Page
@www.route("/hmdb/people", methods=["GET"])
@www.route("/people", methods=["GET"])
def people():
    return render_template("people.html",
                           title="People at HmtDB")


# Contacts Page
@www.route("/hmdb/contacts", methods=["GET"])
@www.route("/contacts", methods=["GET"])
def contacts():
    return render_template("contacts.html",
                           title="Contacts")


# Main Page
@www.route("/hmdb/body", methods=["GET"])
@www.route("/body", methods=["GET"])
def body():
    return render_template("body.html",
                           title="Menu")


# Downloads Page
@www.route("/hmdb/downloads", methods=["GET"])
@www.route("/downloads", methods=["GET"])
def downloads():
    return render_template("downloads.html",
                           title="Downloads")


@www.route("/download_file", methods=["GET"])
def download_file():
    dataset = request.args.get("dataset", "", type=str)

    return send_file("static/zips/" + dataset, mimetype="text/plain",
                     as_attachment=True, attachment_filename=dataset)


@www.route("/download_algs", methods=["GET"])
def download_algs():
    genome_ids = request.args.get("genome_ids", "", type=str)
    id_list = list(map(int, genome_ids.split("_")))
    algs = get_alignments(id_list)
    with open("app/site/alg_dl/hmtdb_algs.fa", "w") as f:
        for el in algs:
            f.write(">{}\n{}\n".format(el[0], el[1]))

    @after_this_request
    def remove_file(response):
        os.remove("app/site/alg_dl/hmtdb_algs.fa")
        return response

    return send_file("site/alg_dl/hmtdb_algs.fa", mimetype="text/plain",
                     as_attachment=True, attachment_filename="hmtdb_algs.fa")


# Site Variability Page
@www.route("/hmdb/siteVariability", methods=["GET"])
@www.route("/siteVariability", methods=["GET"])
def siteVariability():
    return render_template("siteVariability.html",
                           title="Site Variability")


# MToolBox Page
@www.route("/hmdb/mtoolbox", methods=["GET"])
@www.route("/mtoolbox", methods=["GET"])
def mtoolbox():
    return render_template("mtoolbox.html",
                           title="MToolBox")


# Statistics Page
@www.route("/hmdb/stats", methods=["GET"])
@www.route("/stats", methods=["GET"])
def stats():

    return render_template("stats.html",
                           title="Statistics",
                           worlddict=dbdata.worldDict,
                           vardict=dbdata.varDict)


# Query Page
@www.route("/hmdb/query", methods=["GET", "POST"])
@www.route("/query", methods=["GET", "POST"])
def queryCriteria():
    form = QueryForm()
    if request.method == "GET":

        return render_template("queryCriteria.html",
                               title="Query",
                               sources=dbdata.sources,
                               diseases=dbdata.diseases,
                               journals=dbdata.journals)

    elif request.method == "POST":

        return redirect(url_for("site.queryResults",
                                haplotypeHmdb=form.haplotypeHmdb.data,
                                referenceDbId=form.referenceDbId.data,
                                continent=form.continent.data,
                                country=form.country.data,
                                macrohap=form.macrohap.data,
                                # haplogroup=form.haplogroup.data,  # temporarily disabled
                                # haplogroup_user=form.haplogroup_user.data,  # temporarily disabled
                                complete_genome=form.complete_genome.data,
                                snp_position=form.snp_position.data,
                                transit=form.transit.data,
                                transv=form.transv.data,
                                insertion=form.insertion.data,
                                insertion_position=form.insertion_position.data,
                                insertion_sequence=form.insertion_sequence.data,
                                #insertion_length=form.insertion_length.data,
                                deletion=form.deletion.data,
                                start_deletion=form.start_deletion.data,
                                end_deletion=form.end_deletion.data,
                                age=form.age.data,
                                sex=form.sex.data,
                                source=form.source.data,
                                genome_type=form.genome_type.data,
                                disease=form.disease.data,
                                pubmedId=form.pubmedId.data,
                                journal=form.journal.data,
                                author=form.author.data))


# Query Results Page
@www.route("/hmdb/results", methods=["GET", "POST"])
@www.route("/results", methods=["GET", "POST"])
def queryResults():
    if request.method == "GET":

        haplotypeHmdb = request.args.get("haplotypeHmdb", "", type=str)
        referenceDbId = request.args.get("referenceDbId", "", type=str)
        continent = request.args.get("continent", "", type=str)
        country = request.args.get("country", "", type=str)
        macrohap = request.args.get("macrohap", "", type=str)
        # haplogroup = request.args.get("haplogroup", "", type=str)  # temporarily disabled
        # haplogroup_user = request.args.get("haplogroup_user", "", type=str)  # temporarily disabled
        complete_genome = request.args.get("complete_genome", "", type=str)
        snp_position = request.args.get("snp_position", "", type=str)
        transit = request.args.getlist("transit")
        transv = request.args.getlist("transv")
        insertion = request.args.get("insertion")
        insertion_position = request.args.get("insertion_position", "", type=str)
        insertion_sequence = request.args.get("insertion_sequence", "", type=str)
        #insertion_length = request.args.get("insertion_length", 0, type=int)
        deletion = request.args.get("deletion")
        start_deletion = request.args.get("start_deletion", 0, type=int)
        end_deletion = request.args.get("end_deletion", 0, type=int)
        age = request.args.get("age", 0, type=int)
        sex = request.args.get("sex", "", type=str)
        source = request.args.get("source", 0, type=int)
        genome_type = request.args.get("genome_type", "", type=str)
        disease = request.args.getlist("disease")
        pubmedId = request.args.get("pubmedId", "", type=str)
        journal = request.args.get("journal", "", type=str)
        author = request.args.get("author", "", type=str)

        qString = "Genome.query"

        if haplotypeHmdb:
            q = Genome.query.filter(Genome.haplotypeHmdb == haplotypeHmdb.upper()).first()
            return redirect(url_for("site.genomeCard", idGenome=q.genomeId))

        if referenceDbId:
            q = Genome.query.filter(Genome.referenceDbId == referenceDbId.upper()).first()
            return redirect(url_for("site.genomeCard", idGenome=q.genomeId))

        if continent:
            if country:
                countryQ = Genome.query.join(IndividualsData).filter(IndividualsData.countryId == country).subquery()
            else:
                countryQ = Genome.query.join(IndividualsData.query.join(Country, IndividualsData.countryId == Country.countryId).filter(Country.continentCode == continent)).subquery()
            qString += ".join(countryQ, Genome.genomeId == countryQ.c.genomeId)"

        if macrohap:
            # temporarily disabled to fix issue
            # if haplogroup:
            #     haplogroupQ = Genome.query.filter(Genome.haplogroupHmdb == haplogroup).subquery()
            # elif haplogroup_user:
            #     haplogroupQ = Genome.query.filter(Genome.haplogroupUser == haplogroup_user).subquery()
            # else:
            #     haplogroupQ = Genome.query.filter(Genome.haplogroupHmdb.like(macrohap+"%")).subquery()
            haplogroupQ = Genome.query.filter(Genome.haplogroupHmdb.like(macrohap+"%")).subquery()
            qString += ".join(haplogroupQ, Genome.genomeId == haplogroupQ.c.genomeId)"

        if complete_genome and complete_genome != "0":
            completeGenomeQ = Genome.query.filter(Genome.completeGenome == complete_genome.upper()).subquery()
            qString += ".join(completeGenomeQ, Genome.genomeId == completeGenomeQ.c.genomeId)"

        if snp_position:
            if "_" in snp_position:
                positions = snp_position.split("_")
                query_snp = "Genome.query.join(GenomeSnp).filter(or_(GenomeSnp.snpPosition == {}".format(int(positions[0]))
                for n in range(1, len(positions)):
                    query_snp += ", GenomeSnp.snpPosition == {}".format(int(positions[n]))
                query_snp += ")).subquery()"
                snpPositionQ = eval(query_snp)
            elif "-" in snp_position:
                start_pos, end_pos = snp_position.split("-")
                snpPositionQ = Genome.query.join(GenomeSnp).filter(GenomeSnp.snpPosition >= int(start_pos), GenomeSnp.snpPosition <= int(end_pos)).subquery()
            else:
                snpPositionQ = Genome.query.join(GenomeSnp).filter(GenomeSnp.snpPosition == int(snp_position)).subquery()
            qString += ".join(snpPositionQ, Genome.genomeId == snpPositionQ.c.genomeId)"

        if transit:
            if "0" in transit:
                pass
                # print "zero"
                # query all transitions?
            else:
                if len(transit) == 1:
                    query_transit = "Genome.query.join(GenomeSnp).filter(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}').subquery()".format(t=transit[0])

                else:
                    query_transit = "Genome.query.join(GenomeSnp).filter(or_(and_(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}')"
                    temp_str = ", and_(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}')"

                    query_transit = query_transit.format(t=transit[0])
                    for n in range(1, len(transit)):
                        query_transit += temp_str.format(t=transit[n])

                    query_transit += ")).subquery()"
                transitQ = eval(query_transit)
                qString += ".join(transitQ, Genome.genomeId == transitQ.c.genomeId)"

        if transv:
            if "0" in transv:
                pass
                # query all transversions?
            else:
                if len(transv) == 1:
                    query_transv = "Genome.query.join(GenomeSnp).filter(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}').subquery()".format(t=transv[0])

                else:
                    query_transv = "Genome.query.join(GenomeSnp).filter(or_(and_(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}')"
                    temp_str = ", and_(GenomeSnp.rcrsType == '{t[0]}', GenomeSnp.snpType == '{t[1]}')"

                    query_transv = query_transv.format(t=transv[0])
                    for n in range(1, len(transv)):
                        query_transv += temp_str.format(t=transv[n])

                    query_transv += ")).subquery()"
                transvQ = eval(query_transv)
                qString += ".join(transvQ, Genome.genomeId == transvQ.c.genomeId)"

        # todo: to be continued
        # if insertion:
        #    pass

        if insertion_position:
            insertionPositionQ = Genome.query.join(Insertion).filter(Insertion.position5P == int(insertion_position)).subquery()
            qString += ".join(insertionPositionQ, Genome.genomeId == insertionPositionQ.c.genomeId)"

        if insertion_sequence:
            insertionSequenceQ = Genome.query.join(Insertion).filter(Insertion.sequence == insertion_sequence.upper()).subquery()
            qString += ".join(insertionSequenceQ, Genome.genomeId == insertionSequenceQ.c.genomeId)"

        # todo: to be continued
        # if insertion_length:
        #     insertionLengthQ = Genome.query.join(Insertion).filter(Insertion.sequence == insertion_length).subquery()
        #     qString += ".join(insertionLengthQ, Genome.genomeId == insertionLengthQ.c.genomeId)"

        # todo: to be continued
        # if deletion:
        #     if deletion == "False":
        #         #query_del = Deletion.query.filter(Deletion.genomeId != 0).subquery()
        #         #qString += ".join(query_del, Genome.genomeId != query_del.c.genomeId)"
        #         qString += ".except_all(Deletion.query.filter(Deletion.genomeId != 0).subquery())"

        if start_deletion:
            startDeletionQ = Genome.query.join(Deletion).filter(Deletion.fromPosition == int(start_deletion)).subquery()
            qString += ".join(startDeletionQ, Genome.genomeId == startDeletionQ.c.genomeId)"

        if end_deletion:
            endDeletionQ = Genome.query.join(Deletion).filter(Deletion.toPosition == int(end_deletion)).subquery()
            qString += ".join(endDeletionQ, Genome.genomeId == endDeletionQ.c.genomeId)"

        if age:
            ageQ = Genome.query.join(IndividualsData).filter(IndividualsData.age == int(age)).subquery()
            qString += ".join(ageQ, Genome.genomeId == ageQ.c.genomeId)"

        if sex:
            sexQ = Genome.query.join(IndividualsData).filter(IndividualsData.sex == sex).subquery()
            qString += ".join(sexQ, Genome.genomeId == sexQ.c.genomeId)"

        if source:
            sourceQ = Genome.query.filter(Genome.sourceId == int(source)).subquery()
            qString += ".join(sourceQ, Genome.genomeId == sourceQ.c.genomeId)"

        if genome_type and genome_type != "0":
            if genome_type == "P" and disease:
                genomeTypeQ = Genome.query.join(IndividualsData).filter(Genome.genomeType == "P").filter(IndividualsData.diseases == int(disease[0])).subquery()
            else:
                genomeTypeQ = Genome.query.filter(Genome.genomeType == genome_type).subquery()
            qString += ".join(genomeTypeQ, Genome.genomeId == genomeTypeQ.c.genomeId)"

        if pubmedId:
            pubmedQ = Genome.query.join(Reference).filter(Reference.pubmedId == pubmedId).subquery()
            qString += ".join(pubmedQ, Genome.genomeId == pubmedQ.c.genomeId)"

        if journal:
            journalQ = Genome.query.join(Reference).filter(Reference.paper == journal).subquery()
            qString += ".join(journalQ, Genome.genomeId == journalQ.c.genomeId)"

        if author:
            authorQ = Genome.query.join(Reference).filter(Reference.author.like("%"+author+"%")).subquery()
            qString += ".join(authorQ, Genome.genomeId == authorQ.c.genomeId)"


        qString += ".order_by(Genome.genomeId).all()"
        if qString == "Genome.query.order_by(Genome.genomeId).all()":
            return redirect(url_for("site.queryCriteria"))
        mainQuery = eval(qString)

        return render_template("queryResults.html",
                               title="Results",
                               numResults=len(mainQuery),
                               results=mainQuery,
                               genids="_".join([str(el.genomeId)
                                                for el in mainQuery]))


@www.route("/hmdb/genomeCard/<int:idGenome>", methods=["GET", "POST"])
@www.route("/genomeCard/<int:idGenome>", methods=["GET", "POST"])
def genomeCard(idGenome):
    if request.method == "GET":

        genome = Genome.query.filter(Genome.genomeId == idGenome).first()
        source = Sources.query.filter(Sources.sourceId == genome.sourceId).first()
        method = Methods.query.filter(Methods.methodId == genome.methodId).first()
        refs = Reference.query.filter(Reference.genomeId == idGenome).all()
        indiv = IndividualsData.query.filter(IndividualsData.individualId == genome.individualId).first()

        try:
            country = Country.query.filter(Country.countryId == indiv.countryId).first()
            ethnic = EthnicGroups.query.filter(EthnicGroups.groupId == indiv.groupId).first()
            disease = Disease.query.filter(Disease.diseaseId == indiv.diseases).first()
        except:
            country = ""
            ethnic = ""
            disease = ""
        snps = GenomeSnp.query.filter(GenomeSnp.genomeId == idGenome).all()

        return render_template("genomeCard.html",
                               title="Genome Card",
                               genome=genome,
                               source=source,
                               method=method,
                               num_refs=len(refs),
                               refs=refs,
                               indiv=indiv,
                               country=country,
                               ethnic=ethnic,
                               disease=disease,
                               num_snps=len(snps),
                               snps=snps,
                               queryLocus=queryLocus,
                               queryNtVar_N=queryNtVar_N,
                               queryNtVar_P=queryNtVar_P,
                               queryMitomapDna=queryMitomapDna,
                               queryMitomapAa=queryMitomapAa,
                               queryMitomapDnaDiseases=queryMitomapDnaDiseases,
                               queryAaVar_N=queryAaVar_N,
                               queryAaVar_P=queryAaVar_P,
                               queryDisease=queryDisease,
                               queryDeletion=queryDeletion,
                               queryInsertion=queryInsertion,
                               getAltCodon=getAltCodon,
                               aa_dict=aa_dict,
                               getAa=getAa,
                               getId=retrieveHmtVar)

    elif request.method == "POST":
        # TODO: I can't remember why we also need the POST method...maybe I can remove it
        pass


@www.route("/hmdb/ntSitevar/<int:ntPos>", methods=["GET"])
@www.route("/ntSitevar/<int:ntPos>", methods=["GET"])
def ntSitevar(ntPos):

    posN = NtVariability.query.filter(NtVariability.nucleotidePosition == ntPos,
                                      NtVariability.genomeType == "N").first()
    posP = NtVariability.query.filter(NtVariability.nucleotidePosition == ntPos,
                                      NtVariability.genomeType == "P").first()

    return render_template("ntSitevar.html",
                           title="Nt Site Variability",
                           ntPos=ntPos,
                           posN=posN,
                           posP=posP)


@www.route("/hmdb/aaSitevar/<int:aaPos>/<gene>", methods=["GET"])
@www.route("/aaSitevar/<int:aaPos>/<gene>", methods=["GET"])
def aaSitevar(aaPos, gene):

    aaVarN = AaVariability.query.filter(AaVariability.aaPos == aaPos,
                                        AaVariability.geneName == gene,
                                        AaVariability.genomeType == "N").first()
    aaVarP = AaVariability.query.filter(AaVariability.aaPos == aaPos,
                                        AaVariability.geneName == gene,
                                        AaVariability.genomeType == "P").first()

    return render_template("aaSitevar.html",
                           title="Aa Site Variability",
                           aaPos=aaPos,
                           gene=gene,
                           aaVarN=aaVarN,
                           aaVarP=aaVarP)


@www.errorhandler(404)
def page_not_found(e):
    return render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
def internal_server_error(e):
    return render_template("500.html", title="Error 500"), 500

