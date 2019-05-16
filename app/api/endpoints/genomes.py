#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from sqlalchemy import or_
from app.api.views import api
from app.site.models import Genome, Country, IndividualsData, GenomeSnp
from app.api.models import genome_schema, genome_schema_many, genome_schema_comp

ns = api.namespace("genome", description="Retrieve data from the Genome table.")


@ns.deprecated
@ns.route("/")
class GenomeList(Resource):
    def get(self):
        """Get all the entries in Genome table.

        It is not recommended to run this query, as it may load a huge
        number of entries and consequently slow down your browser.
        Will return a list of entries.
        """
        q = Genome.query.all()
        return genome_schema_many.jsonify(q)


@ns.route("/<int:genome_id>")
class GenomeId(Resource):
    def get(self, genome_id):
        """Find the genome associated with the specified id.

        Will return a single entry.
        """
        q = Genome.query.filter(Genome.genomeId == genome_id).first()
        return genome_schema_comp.jsonify(q)


@ns.route("/continent/<cont_name>")
class GenomeContinent(Resource):
    def get(self, cont_name):
        """Find genomes associated with the specified continent name.

        Allowed continent names: [Africa, Asia, America, Europe, Oceania]
        Will return a list of entries.
        """
        cont_q = Genome.query.join(IndividualsData.query.join(Country, IndividualsData.countryId == Country.countryId).filter(Country.continentCode == cont_name)).subquery()
        q = Genome.query.join(cont_q, Genome.genomeId == cont_q.c.genomeId).order_by(Genome.genomeId).all()
        return genome_schema_many.jsonify(q)


@ns.route("/country/<country_name>")
class GenomeCountry(Resource):
    def get(self, country_name):
        """Find genomes associated with the specified country name.

        Will return a list of entries.
        """
        country_q = Genome.query.join(IndividualsData.query.join(Country, IndividualsData.countryId == Country.countryId).filter(Country.countryName == country_name)).subquery()
        q = Genome.query.join(country_q, Genome.genomeId == country_q.c.genomeId).order_by(Genome.genomeId).all()
        return genome_schema_many.jsonify(q)


@ns.route("/haplo/<haplotype>")
class GenomeHaplotype(Resource):
    def get(self, haplotype):
        """Find genomes associated with the specified macrohaplogroup and haplotype.

        Will return a list of entries.
        """
        haplotype_q = Genome.query.filter(Genome.haplogroupHmdb == haplotype).subquery()
        q = Genome.query.join(haplotype_q, Genome.genomeId == haplotype_q.c.genomeId).order_by(Genome.genomeId).all()
        return genome_schema_many.jsonify(q)


@ns.route("/pos/<snp_pos>")
class GenomeSnpPosition(Resource):
    def get(self, snp_pos):
        """Find genomes with the snp(s) in the specified position(s).

        Snps can be in one position (e.g. 263), several positions
        (e.g. 245_2145_11789) or in a specific region (1120-2780).
        Will return a list of entries.
        """
        if "_" in snp_pos:
            positions = snp_pos.split("_")
            query_snp = "Genome.query.join(GenomeSnp).filter(or_(GenomeSnp.snpPosition == %d" % int(positions[0])
            for n in range(1, len(positions)):
                query_snp += ", GenomeSnp.snpPosition == %d" % int(positions[n])
            query_snp += ")).subquery()"
            snpPositionQ = eval(query_snp)
        elif "-" in snp_pos:
            start_pos, end_pos = snp_pos.split("-")
            snpPositionQ = Genome.query.join(GenomeSnp).filter(GenomeSnp.snpPosition >= int(start_pos), GenomeSnp.snpPosition <= int(end_pos)).subquery()
        else:
            snpPositionQ = Genome.query.join(GenomeSnp).filter(GenomeSnp.snpPosition == int(snp_pos)).subquery()
        q = Genome.query.join(snpPositionQ, Genome.genomeId == snpPositionQ.c.genomeId).order_by(Genome.genomeId).all()
        return genome_schema_many.jsonify(q)


@ns.route("/genbank/<gb_id>")
class GenomeGenbankId(Resource):
    def get(self, gb_id):
        """Find the genome with the given Genbank ID.

        Will return a single entry.
        """
        q = Genome.query.filter(Genome.referenceDbId == gb_id).first()
        return genome_schema_comp.jsonify(q)
