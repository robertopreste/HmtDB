#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
from app.site.scripts import retrieveHmtVar, getCountries, getSources, getDiseases


def test_retrieveHmtVar():
    expect = "1255"
    result = retrieveHmtVar(3308, "C")
    assert result == expect


def test_getCountries():
    expect = [(235, 'American Samoa'), (10, 'Australia'), (251, 'Christmas Island'), (252, 'Cocos (Keeling) Islands'), (44, 'Cook Islands'), (118, 'Federal States of Micronesia'), (65, 'Fiji'), (160, 'French Polynesia'), (270, 'Guam'), (274, 'Heard Island and McDonald Islands'), (98, 'Kiribati'), (294, 'Marshall Islands'), (302, 'Nauru'), (304, 'New Caledonia'), (136, 'New Zealand'), (141, 'Niue'), (305, 'Norfold Island'), (306, 'Northern Mariana Islands'), (308, 'Palau'), (154, 'Papua New Guinea'), (309, 'Pitcairn'), (171, 'Samoa'), (182, 'Solomon Islands'), (329, 'Tokelau'), (201, 'Tonga'), (204, 'Tuvalu'), (212, 'Undefined Oceanian Country'), (333, 'United States Minor Outlying Islands'), (218, 'Vanuatu'), (225, 'Wallis and Futuna')]
    result = getCountries("OC")
    assert result == expect


def test_getSources():
    expect = [(24, 'B-Lymphocyte'), (6, 'Benign tissue')]
    result = getSources()[0:2]
    assert result == expect


def test_getDiseases():
    expect = [(1, 'Adenocarcinoma'), (2, 'Adenoma')]
    result = getDiseases()[0:2]
    assert result == expect



