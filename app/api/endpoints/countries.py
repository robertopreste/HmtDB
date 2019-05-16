#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_restplus import Resource
from app.api.views import api
from app.site.models import Country
from app.api.models import country_schema, country_schema_many

ns = api.namespace("country",
                   description="Retrieve data from the Country table.")


@ns.deprecated
@ns.route("/")
class CountryList(Resource):
    def get(self):
        """Get all the entries in Country table.

        It is not recommended to run this query, as it may load a huge
        number of entries and consequently slow down your browser.
        Will return a list of entries.
        """
        q = Country.query.all()
        return country_schema_many.jsonify(q)


@ns.route("/<int:country_id>")
class CountryId(Resource):
    def get(self, country_id):
        """Find the information associated with the specified id.

        Will return a single entry.
        """
        q = Country.query.filter(Country.countryId == country_id).first()
        return country_schema.jsonify(q)

