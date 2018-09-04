#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask import Blueprint
from flask_restplus import Api


res = Blueprint("api", __name__)
api = Api(res, version="1.0", title="HmtDB API", description="A simple API for data hosted on HmtDB.")
