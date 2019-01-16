#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import EthnicGroups
from app.api.models import ethnic_schema, ethnic_schema_many

ns = api.namespace("ethnics", description="Retrieve data from the EthnicGroups table.")


@ns.deprecated
@ns.route("/")
class EthnicList(Resource):
    def get(self):
        """
        Get all the entries in EthnicGroups table.
        It is not recommended to run this query, as it may load a huge number of entries and
        consequently slow down your browser.
        Will return a list of entries.
        """
        q = EthnicGroups.query.all()
        return ethnic_schema_many.jsonify(q)


@ns.route("/<int:group_id>")
class EthnicId(Resource):
    def get(self, group_id):
        """
        Find the information associated with the specified id.
        Will return a single entry.
        """
        q = EthnicGroups.query.filter(EthnicGroups.groupId == group_id).first()
        return ethnic_schema.jsonify(q)
