#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Insertion
from app.api.models import insertion_schema, insertion_schema_many

ns = api.namespace("insertion", description="Retrieve data from the Insertion table.")


@ns.deprecated
@ns.route("/")
class InsertionList(Resource):
    def get(self):
        """
        Get all the entries in Insertion table.
        It is not recommended to run this query, as it may load a huge number of entries and consequently slow down your browser. Will return a list of entries.
        """
        q = Insertion.query.all()
        return insertion_schema_many.jsonify(q)


@ns.route("/<int:insertion_id>")
class InsertionId(Resource):
    def get(self, insertion_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = Insertion.query.filter(Insertion.insertionId == insertion_id).first()
        return insertion_schema.jsonify(q)