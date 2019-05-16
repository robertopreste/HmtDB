#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_marshmallow import Marshmallow
from flask_cors import CORS
import click
import datetime
import pandas as pd

app = Flask(__name__)
CORS(app)
app.config.from_object("config")
Bootstrap(app)
db = SQLAlchemy(app)
ma = Marshmallow(app)

from .api.views import res
from .site.views import www

from .api.views import api
from .api.endpoints.genomes import ns as genome_namespace
from .api.endpoints.individuals import ns as individual_namespace
from .api.endpoints.methods import ns as method_namespace
from .api.endpoints.countries import ns as country_namespace
from .api.endpoints.diseases import ns as disease_namespace
from .api.endpoints.ethnics import ns as ethnic_namespace
from .api.endpoints.sources import ns as source_namespace
from .api.endpoints.genAnnotations import ns as genAnnotation_namespace
from .api.endpoints.deletions import ns as deletion_namespace
from .api.endpoints.insertions import ns as insertion_namespace
from .api.endpoints.genomeSnps import ns as snp_namespace

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

app.register_blueprint(www, static_folder="site/static")
app.register_blueprint(res, url_prefix="/api")


@app.cli.command()
def create_db():
    click.echo("Creating new database... ", nl=False)
    db.drop_all()
    db.create_all()
    click.echo("Done.")


@app.cli.command()
def migrate_db():
    from app.site.scripts import getSources, getDiseases, getJournals, populateCountriesScript, \
        populateHaploScript, populateUserHaploScript, otherFunctions, getStats, getVariants
    click.echo("Migrating db and saving data... ", nl=False)
    with open("app/static/dbdata.py", "w") as d:
        sources = getSources()
        d.write("sources = " + repr(sources) + "\n")

        diseases = getDiseases()
        d.write("diseases = " + repr(diseases) + "\n")

        journals = getJournals()
        d.write("journals = " + repr(journals) + "\n")

        world_dict = getStats()
        d.write("worldDict = " + repr(world_dict) + "\n")

        var_dict = getVariants()
        d.write("varDict = " + repr(var_dict) + "\n")

        latest_update = "%s %s" % (
        datetime.date.today().strftime("%B"), str(datetime.date.today().year))
        d.write("latest_update = " + repr(latest_update) + "\n")
        last_update = "%s_%s" % (
        datetime.date.today().strftime("%m"), str(datetime.date.today().year))
        d.write("last_update = " + repr(last_update) + "\n")

    # create the file script.js with data to populate countries and haplogroup menus
    with open("app/static/js/script.js", "w") as s:
        # countries
        s.write(populateCountriesScript())

        s.write("\n\n")

        # haplogroups
        # temporarily disabled to fix issue
        # s.write(populateHaploScript())
        # s.write(populateUserHaploScript())

        # other functions
        s.write(otherFunctions())

    click.echo("Done.")


@app.cli.command()
def update_db():
    click.echo("Updating database tables...")
    db.engine.execute("SET FOREIGN_KEY_CHECKS=0")
    sources = ("Blosum", "Country", "Disease", "EthnicGroups", "Locus",
               "Methods", "Sources", "Stats", "AaVariability", "Deletion",
               "GenAlignment", "GenAnnotation", "GenomeSnp", "IndividualsData",
               "Insertion", "MitomapAa", "MitomapDna", "NtVariability",
               "Reference")
    for el in sources:
        click.echo("\tUpdating {} table... ".format(el), nl=False)
        df = pd.read_csv("update/data/tables/{}.csv".format(el),
                         keep_default_na=False)
        # df.reset_index(drop=True, inplace=True)
        # df.reset_index(inplace=True)
        # df.rename(columns={"index": "id"}, inplace=True)
        df.to_sql(name=el, con=db.engine, index=False, if_exists="append")
        # df.to_sql(name=el, con=db.engine, index=False, if_exists="replace",
        #           index_label="id")
        click.echo("Complete.")
    # Genome table is a bit complicated
    click.echo("\tUpdating Genome table... ", nl=False)
    df = pd.read_csv("update/data/tables/Genome.csv", keep_default_na=False)
    db.engine.execute("SET SESSION sql_mode='NO_AUTO_VALUE_ON_ZERO'")
    df.to_sql(name="Genome", con=db.engine, index=False, if_exists="append")
    click.echo("Complete.")

    db.engine.execute("SET FOREIGN_KEY_CHECKS=1")
    # click.echo("Database correctly updated. Please migrate it before use.")
    click.echo("Done.")
