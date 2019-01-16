#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import GenAnnotation
from app.api.models import genAnnotation_schema, genAnnotation_schema_many

ns = api.namespace("genAnnotation", description="Retrieve data from the GenAnnotation table.")


@ns.deprecated
@ns.route("/")
class GenAnnotationList(Resource):
    def get(self):
        """
        Get all the entries in GenAnnotation table.
        It is not recommended to run this query, as it may load a huge number of entries and
        consequently slow down your browser.
        Will return a list of entries.
        """
        q = GenAnnotation.query.all()
        return genAnnotation_schema_many.jsonify(q)


@ns.route("/<int:annotation_id>")
class GenAnnotationId(Resource):
    def get(self, annotation_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = GenAnnotation.query.filter(GenAnnotation.annotationId == annotation_id).first()
        return genAnnotation_schema.jsonify(q)
