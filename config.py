#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os

basedir = os.path.abspath(os.path.dirname(__file__))

WTF_CSRF_ENABLED = True
SECRET_KEY = "secret_key"

SQLALCHEMY_TRACK_MODIFICATIONS = False

# database
SQLALCHEMY_DATABASE_URI = "sqlite:///" + os.path.join(basedir, "hmtdb.db")
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, "db_repo")

# administrators
ADMINS = ["admin@mail.com"]
