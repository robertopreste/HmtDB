#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import imp
import datetime
from migrate.versioning import api
from app import db
from config import SQLALCHEMY_DATABASE_URI
from config import SQLALCHEMY_MIGRATE_REPO
from app.site.scripts import getSources, getDiseases, getJournals, populateCountriesScript, populateHaploScript, populateUserHaploScript, otherFunctions, getStats, getVariants

v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
migration = SQLALCHEMY_MIGRATE_REPO + ("/versions/%03d_migration.py" % (v + 1))
tmp_module = imp.new_module("old_model")
old_model = api.create_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

exec(old_model, tmp_module.__dict__)

script = api.make_update_script_for_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, tmp_module.meta, db.metadata)
open(migration, "wt").write(script)
api.upgrade(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

# questo crea il file dbdata.py con i dati per popolare i menu a tendina
with open("db_repo/dbdata.py", "w") as d:

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

    # il seguente serve per visualizzare nella home page quando è stato aggiornato il db
    latest_update = "%s %s" % (datetime.date.today().strftime("%B"), str(datetime.date.today().year))
    d.write("latest_update = " + repr(latest_update) + "\n")
    # il seguente serve per l'update del db
    last_update = "%s_%s" % (datetime.date.today().strftime("%m"), str(datetime.date.today().year))
    d.write("last_update = " + repr(last_update) + "\n")


# questo crea il file script.js con i dati per popolare i menu countries e haplogroup
with open("app/static/js/script.js", "w") as s:
    # countries
    s.write(populateCountriesScript())

    s.write("\n\n")

    # haplogroups
    # debug: temporarily disabled to fix issue
    # s.write(populateHaploScript())
    # s.write(populateUserHaploScript())

    # other functions
    s.write(otherFunctions())


print("New migration saved as " + migration)
print("Current database version: " + str(v))
