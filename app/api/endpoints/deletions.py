#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Deletion
from app.api.models import deletion_schema, deletion_schema_many

ns = api.namespace("deletion",
                   description="Retrieve data from the Deletion table.")


@ns.deprecated
@ns.route("/")
class DeletionList(Resource):
    def get(self):
        """Get all the entries in Deletion table.

        It is not recommended to run this query, as it may load a huge
        number of entries and consequently slow down your browser.
        Will return a list of entries.
        """
        q = Deletion.query.all()
        return deletion_schema_many.jsonify(q)


@ns.route("/<int:deletion_id>")
class DeletionId(Resource):
    def get(self, deletion_id):
        """Find the information associated with the specified id.

        Will return a single entry.
        """
        q = Deletion.query.filter(Deletion.deletionId == deletion_id).first()
        return deletion_schema.jsonify(q)
