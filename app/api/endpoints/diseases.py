#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Disease
from app.api.models import disease_schema, disease_schema_many

ns = api.namespace("disease", description="Retrieve data from the Disease table.")


@ns.deprecated
@ns.route("/")
class DiseaseList(Resource):
    def get(self):
        """
        Get all the entries in Disease table.
        It is not recommended to run this query, as it may load a huge number of entries and consequently slow down your browser. Will return a list of entries.
        """
        q = Disease.query.all()
        return disease_schema_many.jsonify(q)


@ns.route("/<int:disease_id>")
class DiseaseId(Resource):
    def get(self, disease_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = Disease.query.filter(Disease.diseaseId == disease_id).first()
        return disease_schema.jsonify(q)