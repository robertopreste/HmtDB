#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask import Blueprint, url_for
from flask_restplus import Api


class MyApi(Api):
    @property
    def specs_url(self):
        """Monkey patch for HTTPS"""
        scheme = "http" if "5000" in self.base_url else "https"
        return url_for(self.endpoint("specs"), _external=True, _scheme=scheme)


res = Blueprint("api", __name__)
api = MyApi(res, version="1.0", title="HmtDB API",
            description="A simple API for data hosted on HmtDB.")

