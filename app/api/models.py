#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app.site.models import Country, Genome, IndividualsData, Methods, Disease, EthnicGroups, Sources, GenAnnotation, Deletion, Insertion, GenomeSnp
from app import ma


class CountrySchema(ma.ModelSchema):
    class Meta:
        model = Country


class DiseaseSchema(ma.ModelSchema):
    class Meta:
        model = Disease


class EthnicSchema(ma.ModelSchema):
    class Meta:
        model = EthnicGroups


class SourceSchema(ma.ModelSchema):
    class Meta:
        model = Sources


class GenAnnotationSchema(ma.ModelSchema):
    class Meta:
        model = GenAnnotation


class DeletionSchema(ma.ModelSchema):
    class Meta:
        model = Deletion


class InsertionSchema(ma.ModelSchema):
    class Meta:
        model = Insertion


class GenomeSnpSchema(ma.ModelSchema):
    class Meta:
        model = GenomeSnp


class GenomeSchema(ma.ModelSchema):
    class Meta:
        model = Genome


class IndividualSchema(ma.ModelSchema):
    class Meta:
        model = IndividualsData

    Country = ma.Nested(CountrySchema, exclude=("countryId", "countryCode",
                                                "continentCode", "individuals"))
    Disease = ma.Nested(DiseaseSchema, exclude=("diseaseId", "individuals",
                                                "mitoAa", "mitoDna"))
    EthnicGroups = ma.Nested(EthnicSchema, exclude=("groupId",
                                                    "groupDescription",
                                                    "individuals"))


class MethodSchema(ma.ModelSchema):
    class Meta:
        model = Methods


class GenomeSchemaComplete(ma.ModelSchema):
    class Meta:
        model = Genome

    IndividualsData = ma.Nested(IndividualSchema, exclude=("genomeId", ))
    Methods = ma.Nested(MethodSchema, exclude=("methodId", "genomes"))
    Sources = ma.Nested(SourceSchema, exclude=("sourceId", "genomes"))
    annotationId = ma.Nested(GenAnnotationSchema, exclude=("Genome", ))
    deletions = ma.Nested(DeletionSchema, many=True, exclude=("Genome", ))
    insertions = ma.Nested(InsertionSchema, many=True, exclude="Genome", )
    snps = ma.Nested(GenomeSnpSchema, many=True, exclude=("Genome", ))


genome_schema = GenomeSchema()
genome_schema_many = GenomeSchema(many=True)
genome_schema_comp = GenomeSchemaComplete()
genome_schema_comp_many = GenomeSchemaComplete(many=True)
country_schema = CountrySchema()
country_schema_many = CountrySchema(many=True)
individual_schema = IndividualSchema()
individual_schema_many = IndividualSchema(many=True)
method_schema = MethodSchema()
method_schema_many = MethodSchema(many=True)
disease_schema = DiseaseSchema()
disease_schema_many = DiseaseSchema(many=True)
ethnic_schema = EthnicSchema()
ethnic_schema_many = EthnicSchema(many=True)
source_schema = SourceSchema()
source_schema_many = SourceSchema(many=True)
genAnnotation_schema = GenAnnotationSchema()
genAnnotation_schema_many = GenAnnotationSchema(many=True)
deletion_schema = DeletionSchema()
deletion_schema_many = DeletionSchema(many=True)
insertion_schema = InsertionSchema()
insertion_schema_many = InsertionSchema(many=True)
snp_schema = GenomeSnpSchema()
snp_schema_many = GenomeSnpSchema(many=True)
