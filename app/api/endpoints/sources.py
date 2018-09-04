#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Sources
from app.api.models import source_schema, source_schema_many

ns = api.namespace("sources", description="Retrieve data from the Sources table.")


@ns.deprecated
@ns.route("/")
class SourceList(Resource):
    def get(self):
        """
        Get all the entries in Sources table.
        It is not recommended to run this query, as it may load a huge number of entries and consequently slow down your browser. Will return a list of entries.
        """
        q = Sources.query.all()
        return source_schema_many.jsonify(q)


@ns.route("/<int:source_id>")
class SourceId(Resource):
    def get(self, source_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = Sources.query.filter(Sources.sourceId == source_id).first()
        return source_schema.jsonify(q)