#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_marshmallow import Marshmallow
from flask_cors import CORS
from flask_login import LoginManager
from flask_mail import Mail

app = Flask(__name__)
CORS(app)
app.config.from_object("config")
Bootstrap(app)
db = SQLAlchemy(app)
ma = Marshmallow(app)
login = LoginManager(app)
login.login_view = "site.login"
mail = Mail(app)

from api.views import res
from site.views import mod

from api.views import api
from api.endpoints.genomes import ns as genome_namespace
from api.endpoints.individuals import ns as individual_namespace
from api.endpoints.methods import ns as method_namespace
from api.endpoints.countries import ns as country_namespace
from api.endpoints.diseases import ns as disease_namespace
from api.endpoints.ethnics import ns as ethnic_namespace
from api.endpoints.sources import ns as source_namespace
from api.endpoints.genAnnotations import ns as genAnnotation_namespace
from api.endpoints.deletions import ns as deletion_namespace
from api.endpoints.insertions import ns as insertion_namespace
from api.endpoints.genomeSnps import ns as snp_namespace

api.add_namespace(genome_namespace)
api.add_namespace(individual_namespace)
api.add_namespace(method_namespace)
api.add_namespace(country_namespace)
api.add_namespace(disease_namespace)
api.add_namespace(ethnic_namespace)
api.add_namespace(source_namespace)
api.add_namespace(genAnnotation_namespace)
api.add_namespace(deletion_namespace)
api.add_namespace(insertion_namespace)
api.add_namespace(snp_namespace)

app.register_blueprint(mod, static_folder="site/static")
app.register_blueprint(res, url_prefix="/api")
