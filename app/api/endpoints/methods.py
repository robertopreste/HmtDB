#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Methods
from app.api.models import method_schema, method_schema_many

ns = api.namespace("methods",
                   description="Retrieve data from the Methods table.")


@ns.deprecated
@ns.route("/")
class MethodsList(Resource):
    def get(self):
        """Get all the entries in Methods table.

        It is not recommended to run this query, as it may load a huge
        number of entries and consequently slow down your browser.
        Will return a list of entries.
        """
        q = Methods.query.all()
        return method_schema_many.jsonify(q)


@ns.route("/<int:method_id>")
class MethodsId(Resource):
    def get(self, method_id):
        """Find the information associated with the specified id.

        Will return a single entry.
        """
        q = Methods.query.filter(Methods.methodId == method_id).first()
        return method_schema.jsonify(q)

