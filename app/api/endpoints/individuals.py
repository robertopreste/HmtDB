#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import IndividualsData
from app.api.models import individual_schema, individual_schema_many

ns = api.namespace("individuals",
                   description="Retrieve data from the IndividualsData table.")


@ns.deprecated
@ns.route("/")
class IndividualsList(Resource):
    def get(self):
        """Get all the entries in IndividualsData table.

        It is not recommended to run this query, as it may load a huge
        number of entries and consequently slow down your browser.
        Will return a list of entries.
        """
        q = IndividualsData.query.all()
        return individual_schema_many.jsonify(q)


@ns.route("/<int:individual_id>")
class IndividualsId(Resource):
    def get(self, individual_id):
        """Find the information associated with the specified id.

        Will return a single entry.
        """
        q = IndividualsData.query.filter(IndividualsData.individualId == individual_id).first()
        return individual_schema.jsonify(q)
