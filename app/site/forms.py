#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_wtf import Form
from wtforms import StringField, SelectField, SelectMultipleField, RadioField, \
    BooleanField, widgets
from wtforms.validators import NumberRange


class MultiCheckboxField(SelectMultipleField):
    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()


class QueryForm(Form):
    haplotypeHmdb = StringField("haplotypeHmdb", default=None)
    referenceDbId = StringField("referenceDbId", default=None)
    continent = SelectField("continent", default=None)
    country = SelectField("country", default=None)
    macrohap = SelectField("macrohap", default=None)
    haplogroup = SelectField("haplogroup", default=None)
    haplogroup_user = SelectField("haplogroup_user", default=None)
    complete_genome = RadioField("complete_genome",
                                 choices=[("0", "Whole Database"),
                                          ("Y", "Complete Genomes"),
                                          ("N", "Coding Regions")],
                                 default="0")
    snp_position = StringField("snp_position",
                               validators=[NumberRange(min=1, max=17000)],
                               default=None)
    transit = MultiCheckboxField("transit",
                                 choices=[("0", "All"), ("AG", "AG"),
                                          ("GA", "GA"), ("CT", "CT"),
                                          ("TC", "TC")],
                                 default="0")
    transv = MultiCheckboxField("transv",
                                choices=[("0", "All"), ("AC", "AC"),
                                         ("AT", "AT"), ("CA", "CA"),
                                         ("CG", "CG"), ("GC", "GC"),
                                         ("GT", "GT"), ("TA", "TA"),
                                         ("TG", "TG")],
                                default="0")
    insertion = BooleanField("insertion", default=True)
    insertion_position = StringField("insertion_position",
                                     validators=[NumberRange(min=1, max=17000)],
                                     default=None)
    insertion_sequence = StringField("insertion_sequence", default=None)
    #insertion_length = StringField("insertion_length", default=None)
    deletion = BooleanField("deletion", default=True)
    start_deletion = StringField("start_deletion",
                                 validators=[NumberRange(min=1, max=17000)],
                                 default=None)
    end_deletion = StringField("end_deletion",
                               validators=[NumberRange(min=1, max=17000)],
                               default=None)
    age = StringField("age", validators=[NumberRange(min=1, max=100)],
                      default=None)
    sex = RadioField("sex",
                     choices=[("0", "All"), ("M", "Male"), ("F", "Female")],
                     default="0")
    source = SelectField("source", default=None)
    genome_type = RadioField("genome_type",
                             choices=[("0", "All"), ("N", "Healthy"),
                                      ("P", "Pathologic")],
                             default="0")
    disease = SelectMultipleField("disease", default=None)
    pubmedId = StringField("pubmedId", default=None)
    journal = SelectField("journal", default=None)
    author = StringField("author", default=None)

    def __init__(self, *args, **kwargs):
        Form.__init__(self, *args, **kwargs)

