#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import GenomeSnp
from app.api.models import snp_schema, snp_schema_many

ns = api.namespace("snps", description="Retrieve data from the GenomeSnp table.")


@ns.deprecated
@ns.route("/")
class GenomeSnpList(Resource):
    def get(self):
        """
        Get all the entries in GenomeSnp table.
        It is not recommended to run this query, as it may load a huge number of entries and consequently slow down your browser. Will return a list of entries.
        """
        q = GenomeSnp.query.all()
        return snp_schema_many.jsonify(q)


@ns.route("/<int:snp_id>")
class GenomeSnpId(Resource):
    def get(self, snp_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = GenomeSnp.query.filter(GenomeSnp.snpId == snp_id).first()
        return snp_schema.jsonify(q)