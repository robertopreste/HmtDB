#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from datetime import datetime
from threading import Thread
from app import db, login, mail, app
from flask import render_template
from flask_login import UserMixin
from flask_mail import Message
from werkzeug.security import generate_password_hash, check_password_hash
from config import ADMINS


deletion_genome = db.Table("deletion_genome",
                           db.Column("deletion_id", db.Integer, db.ForeignKey("Deletion.deletionId")),
                           db.Column("genome_id", db.Integer, db.ForeignKey("Genome.genomeId")))

disease_individuals = db.Table("disease_individuals",
                               db.Column("disease_id", db.Integer, db.ForeignKey("Disease.diseaseId")),
                               db.Column("individual_id", db.Integer, db.ForeignKey("IndividualsData.individualId")))

genome_reference = db.Table("genome_reference",
                            db.Column("genome_id", db.Integer, db.ForeignKey("Genome.genomeId")),
                            db.Column("reference_id", db.Integer, db.ForeignKey("Reference.referenceId")))

genome_snps = db.Table("genome_snps",
                       db.Column("genome_id", db.Integer, db.ForeignKey("Genome.genomeId")),
                       db.Column("snps_id", db.Integer, db.ForeignKey("GenomeSnp.snpId")))


class AaVariability(db.Model):
    __tablename__ = "AaVariability"

    aaId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    aaPos = db.Column(db.Integer, nullable=False)
    rcrsAa = db.Column(db.String, nullable=True, default=None)
    varAa_intrahs = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_intermam = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_eu = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_am = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_af = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_as = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_oc = db.Column(db.Float, nullable=True, default=0.00000)
    genomeType = db.Column(db.String, nullable=False, default="N")
    # relationship with Locus - "many" side
    geneName = db.Column(db.String, db.ForeignKey("Locus.geneName"), nullable=False)

    def __repr__(self):
        return """AaVariability(aaId: {self.aaId}, aaPos: {self.aaPos}, rcrsAa: {self.rcrsAa}, varAa_intrahs: {self.varAa_intrahs}, varAa_intermam: {self.varAa_intermam}, varAa_eu: {self.varAa_eu}, varAa_am: {self.varAa_am}, varAa_af: {self.varAa_af}, varAa_as: {self.varAa_as}, varAa_oc: {self.varAa_oc}, genomeType: {self.genomeType}, geneName: {self.geneName})""".format(self=self)


class Blosum(db.Model):
    __tablename__ = "Blosum"

    blosumId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    aaChange = db.Column(db.String, nullable=False)
    blosumIndex = db.Column(db.Numeric(6, 3), nullable=False)

    def __repr__(self):
        return """Blosum(blosumId: {self.blosumId}, aaChange: {self.aaChange}, blosumIndex: {self.blosumIndex})""".format(self=self)


class Country(db.Model):
    __tablename__ = "Country"

    countryId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    countryName = db.Column(db.String, nullable=False)
    countryCode = db.Column(db.String, nullable=False, default="XX")
    continentName = db.Column(db.String, nullable=False)
    continentCode = db.Column(db.String, nullable=False, default="XX")
    # relationship with IndividualsData - "one" side
    individuals = db.relationship("IndividualsData", backref="Country", lazy="dynamic")

    def __repr__(self):
        return """Country(countryId: {self.countryId}, countryName: {self.countryName}, countryCode: {self.countryCode}, continentName: {self.continentName}, continentCode: {self.continentCode}, individuals: {self.individuals})""".format(self=self)


class Deletion(db.Model):
    __tablename__ = "Deletion"

    deletionId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    fromPosition = db.Column(db.Integer, nullable=False)
    toPosition = db.Column(db.Integer, nullable=False)
    # relationship with Genome (M -- M) (ora M)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """Deletion(deletionId: {self.deletionId}, fromPosition: {self.fromPosition}, toPosition: {self.toPosition}, genomeId: {self.genomeId})""".format(self=self)


class Disease(db.Model):
    __tablename__ = "Disease"

    diseaseId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    diseaseName = db.Column(db.String, nullable=False)
    diseaseAcronym = db.Column(db.String, nullable=True, default=None)
    url = db.Column(db.String, nullable=True, default=None)
    # relationship with IndividualsData (M -- M) (ora 1)
    individuals = db.relationship("IndividualsData", backref="Disease", lazy="dynamic")
    # relationships with MitomapAa and MitomapDna - one side
    mitoAa = db.relationship("MitomapAa", backref="Disease", lazy="dynamic")
    mitoDna = db.relationship("MitomapDna", backref="Disease", lazy="dynamic")

    def __repr__(self):
        return """Disease(diseaseId: {self.diseaseId}, diseaseName: {self.diseaseName}, diseaseAcronym: {self.diseaseAcronym}, url: {self.url}, individuals: {self.individuals}, mitoAa: {self.mitoAa}, mitoDna: {self.mitoDna})""".format(self=self)


class EthnicGroups(db.Model):
    __tablename__ = "EthnicGroups"

    groupId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    groupName = db.Column(db.String, nullable=True, default=None)
    groupDescription = db.Column(db.String, nullable=True, default=None)
    # relationship with IndividualsData - "one" side
    individuals = db.relationship("IndividualsData", backref="EthnicGroups", lazy="dynamic")

    def __repr__(self):
        return """EthnicGroups(groupId: {self.groupId}, groupName: {self.groupName}, groupDescription: {self.groupDescription}, individuals: {self.individuals})""".format(self=self)


class GenAlignment(db.Model):
    __tablename__ = "GenAlignment"

    alignmentId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    alignment = db.Column(db.String, nullable=False)
    # relationship with Genome - (one to one) convertita a Many
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """GenAlignment(alignmentId: {self.alignmentId}, alignment: {self.alignment}, genomeId: {self.genomeId})""".format(self=self)


class GenAnnotation(db.Model):
    __tablename__ = "GenAnnotation"

    annotationId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    annotations = db.Column(db.String, nullable=False)
    # relationship with Genome - (one to one) convertita a Many
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """GenAnnotation(annotationId: {self.annotationId}, annotations: {self.annotations}, genomeId: {self.genomeId})""".format(self=self)


class Genome(db.Model):
    __tablename__ = "Genome"

    genomeId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    genomeSequence = db.Column(db.String, nullable=False)
    completeGenome = db.Column(db.String, nullable=False)
    startPosition = db.Column(db.Integer, nullable=True, default=None)
    endPosition = db.Column(db.Integer, nullable=True, default=None)
    haplotypeUser = db.Column(db.String, nullable=True, default=None)
    haplotypeHmdb = db.Column(db.String, nullable=True, default=None)
    haplogroupUser = db.Column(db.String, nullable=True, default=None)
    haplogroupHmdb = db.Column(db.String, nullable=True, default=None)
    referenceDb = db.Column(db.String, nullable=True, default=None)
    referenceDbId = db.Column(db.String, nullable=True, default=None)
    genomeType = db.Column(db.String, nullable=False, default="N")
    sourceDbId = db.Column(db.String, nullable=True, default=None)
    # relationship with IndividualsData - (one to one) convertita a Many
    individualId = db.Column(db.Integer, db.ForeignKey("IndividualsData.individualId"), nullable=True, default=None)
    # relationship with GenAlignment - (one to one) convertita a One
    alignmentId = db.relationship("GenAlignment", uselist=False, backref="Genome")
    # relationship with GenAnnotation - (one to one) convertita a One
    annotationId = db.relationship("GenAnnotation", uselist=False, backref="Genome")
    # relationship with Methods - "many" side
    methodId = db.Column(db.Integer, db.ForeignKey("Methods.methodId"), nullable=True, default=None)
    # relationship with Sources - "many" side
    sourceId = db.Column(db.Integer, db.ForeignKey("Sources.sourceId"), nullable=True, default=None)
    # relationship with References - "one" side
    references = db.relationship("Reference", backref="Genome", lazy="dynamic")
    # relationship with GenomeSnp - "one" side (andrebbe fatto M -- M)
    snps = db.relationship("GenomeSnp", backref="Genome", lazy="dynamic")
    # relationship with Deletion (M -- M) (ora 1)
    deletions = db.relationship("Deletion", backref="Genome", lazy="dynamic")
    # relationship with Insertion - ora 1 (andrebbe fatto M -- M)
    insertions = db.relationship("Insertion", backref="Genome", lazy="dynamic")

    def __repr__(self):
        return """Genome(genomeId: {self.genomeId}, genomeSequence: {self.genomeSequence}, completeGenome: {self.completeGenome}, startPosition: {self.startPosition}, endPosition: {self.endPosition}, haplotypeUser: {self.haplotypeUser}, haplotypeHmdb: {self.haplotypeHmdb}, haplogroupUser: {self.haplogroupUser}, haplogroupHmdb: {self.haplogroupUser}, referenceDb: {self.referenceDb}, referenceDbId: {self.referenceDbId}, genomeType: {self.genomeType}, sourceDbId: {self.sourceDbId}, individualId: {self.individualId}, alignmentId: {self.alignmentId}, annotationId: {self.annotationId}, methodId: {self.methodId}, sourceId: {self.sourceId}, references: {self.references}, snps: {self.snps}, deletions: {self.deletions}, insertions: {self.insertions})""".format(self=self)


class GenomeSnp(db.Model):
    __tablename__ = "GenomeSnp"

    snpId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    snpPosition = db.Column(db.Integer, nullable=False)
    snpType = db.Column(db.String, nullable=False)
    rcrsType = db.Column(db.String, nullable=False, default="")
    # relationship with Genome - "many" side (andrebbe fatto M -- M)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """GenomeSnp(snpId: {self.snpId}, snpPosition: {self.snpPosition}, snpType: {self.snpType}, rcrsType: {self.rcrsType}, genomeId: {self.genomeId})""".format(self=self)


class IndividualsData(db.Model):
    __tablename__ = "IndividualsData"

    individualId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    age = db.Column(db.Integer, nullable=True, default=None)
    sex = db.Column(db.String, nullable=True, default=None)
    individualType = db.Column(db.String, nullable=True, default=None)
    # relationship with Genome - (one to one) convertita a One
    genomeId = db.relationship("Genome", backref="IndividualsData", lazy="dynamic")
    # relationship with Country - "many" side
    countryId = db.Column(db.Integer, db.ForeignKey("Country.countryId"), nullable=True, default=None)
    # relationship with EthnicGroups - "many" side
    groupId = db.Column(db.Integer, db.ForeignKey("EthnicGroups.groupId"), nullable=True, default=None)
    # relationship with Disease (M -- M) (ora M)
    diseases = db.Column(db.String, db.ForeignKey("Disease.diseaseId"), nullable=True, default=None)

    def __repr__(self):
        return """IndividualsData(individualId: {self.individualId}, age: {self.age}, sex: {self.sex}, individualType: {self.individualType}, genomeId: {self.genomeId}, countryId: {self.countryId}, groupId: {self.groupId}, diseases: {self.diseases})""".format(self=self)


class Insertion(db.Model):
    __tablename__ = "Insertion"

    insertionId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    position5P = db.Column(db.Integer, nullable=False)
    sequence = db.Column(db.String, nullable=True)
    # relationship with Genome - 1 side
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """Insertion(insertionId: {self.insertionId}, position5P: {self.position5P}, sequence: {self.sequence}, genomeId: {self.genomeId})""".format(self=self)


class Locus(db.Model):
    __tablename__ = "Locus"

    locusId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    geneName = db.Column(db.String, nullable=False, default="")
    locusType = db.Column(db.String, nullable=False)
    startPosition = db.Column(db.Integer, nullable=False)
    endPosition = db.Column(db.Integer, nullable=False)
    startPattern = db.Column(db.String, nullable=True, default=None)
    endPattern = db.Column(db.String, nullable=True, default=None)
    rcrsAaSeq = db.Column(db.String, nullable=True, default=None)
    rcrsDnaSeq = db.Column(db.String, nullable=True, default=None)
    # relationship with MitomapAa - "one" side
    mitomapAa = db.relationship("MitomapAa", backref="Locus", lazy="dynamic")
    # relationship with MitomapDna - "one" side
    mitomapDna = db.relationship("MitomapDna", backref="Locus", lazy="dynamic")
    # relationship with AaVariability - "one" side
    aaVariability = db.relationship("AaVariability", backref="Locus", lazy="dynamic")
    # relationship with NtVariability - "one" side
    ntVariability = db.relationship("NtVariability", backref="Locus", lazy="dynamic")


    def __repr__(self):
        return """Locus(locusId: {self.locusId}, geneName: {self.geneName}, locusType: {self.locusType}, startPosition: {self.startPosition}, endPosition: {self.endPosition}, startPattern: {self.startPattern}, endPattern: {self.endPattern}, rcrsAaSeq: {self.rcrsAaSeq}, mitomapAa: {self.mitomapAa}, mitomapDna: {self.mitomapDna}, aaVariability: {self.aaVariability}, ntVariability: {self.ntVariability})""".format(self=self)


class Methods(db.Model):
    __tablename__ = "Methods"

    methodId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    methodName = db.Column(db.String, nullable=False)
    methodDescription = db.Column(db.String)
    # relationship with Genome - "one" side
    genomes = db.relationship("Genome", backref="Methods", lazy="dynamic")

    def __repr__(self):
        return """Methods(methodId: {self.methodId}, methodName: {self.methodName}, methodDescription: {self.methodDescription}, genomes: {self.genomes})""".format(self=self)


class MitomapAa(db.Model):
    __tablename__ = "MitomapAa"

    varAaId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    aaChange = db.Column(db.String, nullable=False, default="SYNON")
    aaPosition = db.Column(db.Integer, nullable=False)
    variationType = db.Column(db.Integer, nullable=True, default=None)
    # relationship with Locus - "many" side
    geneName = db.Column(db.String, db.ForeignKey("Locus.geneName"), nullable=False)
    # relationship with MitomapDna - "many" side
    varDnaId1 = db.Column(db.Integer, nullable=False)
    varDnaId2 = db.Column(db.Integer, nullable=True, default=None)
    varDnaId3 = db.Column(db.Integer, nullable=True, default=None)
    # relationship with Disease - association table mitomapAaVSdisease
    diseases = db.Column(db.String, db.ForeignKey("Disease.diseaseId"), nullable=True, default=None)

    def __repr__(self):
        return """MitomapAa(varAaId: {self.varAaId}, aaChange: {self.aaChange}, aaPosition: {self.aaPosition}, variationType: {self.variationType}, geneName: {self.geneName}, varDnaId1: {self.varDnaId1}, varDnaId2: {self.varDnaId2}, varDnaId3: {self.varDnaId3}, diseases: {self.diseases})""".format(self=self)


class MitomapDna(db.Model):
    __tablename__ = "MitomapDna"

    varDnaId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    nucleotidePosition = db.Column(db.Integer, nullable=False)
    snpType = db.Column(db.String, nullable=False)
    rcrsType = db.Column(db.String, nullable=False, default="")
    # relationship with Locus - "many" side
    geneName = db.Column(db.Integer, db.ForeignKey("Locus.geneName"), nullable=True, default=None)
    # relationship with Disease - association table mitomapDnaVSdisease
    diseases = db.Column(db.String, db.ForeignKey("Disease.diseaseId"), nullable=True, default=None)

    def __repr__(self):
        return """MitomapDna(varDnaId: {self.varDnaId}, nucleotidePosition: {self.nucleotidePosition}, snpType: {self.snpType}, rcrsType: {self.rcrsType}, geneName: {self.geneName}, diseases: {self.diseases})""".format(self=self)


class NtVariability(db.Model):
    __tablename__ = "NtVariability"

    ntId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    nucleotidePosition = db.Column(db.Integer, nullable=False)
    insertionPosition = db.Column(db.Integer, nullable=False)
    var_tot = db.Column(db.Float, nullable=True, default=0.00000)
    var_eu = db.Column(db.Float, nullable=True, default=0.00000)
    var_am = db.Column(db.Float, nullable=True, default=0.00000)
    var_af = db.Column(db.Float, nullable=True, default=0.00000)
    var_as = db.Column(db.Float, nullable=True, default=0.00000)
    var_oc = db.Column(db.Float, nullable=True, default=0.00000)
    compVar_tot = db.Column(db.String, nullable=True, default=None)
    compVar_eu = db.Column(db.String, nullable=True, default=None)
    compVar_am = db.Column(db.String, nullable=True, default=None)
    compVar_af = db.Column(db.String, nullable=True, default=None)
    compVar_as = db.Column(db.String, nullable=True, default=None)
    compVar_oc = db.Column(db.String, nullable=True, default=None)
    genomeType = db.Column(db.String, nullable=False, default="N")
    # relationship with Locus - "many" side
    # questo va inserito
    geneName = db.Column(db.Integer, db.ForeignKey("Locus.geneName"), nullable=True, default=None)

    def __repr__(self):
        return """NtVariability(ntId: {self.ntId}, nucleotidePosition: {self.nucleotidePosition}, insertionPosition: {self.insertionPosition}, var_tot: {self.var_tot}, var_eu: {self.var_eu}, var_am: {self.var_am}, var_af: {self.var_af}, var_as: {self.var_as}, var_oc: {self.var_oc}, compVar_tot: {self.compVar_tot}, compVar_eu: {self.compVar_eu}, compVar_am: {self.compVar_am}, compVar_af: {self.compVar_af}, compVar_as: {self.compVar_as}, compVar_oc: {self.compVar_oc}, genomeType: {self.genomeType}, geneName: {self.geneName})""".format(self=self)


class Reference(db.Model):
    __tablename__ = "Reference"

    referenceId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    pubmedId = db.Column(db.String, nullable=True, default=None)
    author = db.Column(db.String, nullable=True, default=None)
    title = db.Column(db.String, nullable=True, default=None)
    volume = db.Column(db.String, nullable=True, default=None)
    year = db.Column(db.String, nullable=True, default=None)
    firstPage = db.Column(db.String, nullable=True, default=None)
    paper = db.Column(db.String, nullable=True, default=None)
    issue = db.Column(db.String, nullable=True, default=None)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"), nullable=False)

    def __repr__(self):
        return """Reference(referenceId: {self.referenceId}, pubmedId: {self.pubmedId}, author: {self.author}, title: {self.title}, volume: {self.volume}, year: {self.year}, firstPage: {self.firstPage}, paper: {self.paper}, issue: {self.issue}, genomeId: {self.genomeId})""".format(self=self)


class Sources(db.Model):
    __tablename__ = "Sources"

    sourceId = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    sourceName = db.Column(db.String, nullable=False)
    # relationship with Genome - "one" side
    genomes = db.relationship("Genome", backref="Sources", lazy="dynamic")

    def __repr__(self):
        return """Sources(sourceId: {self.sourceId}, sourceName: {self.sourceName}, genomes: {self.genomes})""".format(self=self)


class Stats(db.Model):
    __tablename__ = "Stats"

    id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    continentName = db.Column(db.String, nullable=False)
    completeGenome = db.Column(db.String, nullable=False)
    genomeType = db.Column(db.String, nullable=False, default="N")
    total = db.Column(db.Integer, nullable=False, default="0")

    def __repr__(self):
        return """Stats(id: {self.id}, continentName: {self.continentName}, completeGenome: {self.completeGenome}, genomeType: {self.genomeType}, total: {self.total})""".format(self=self)


class User(UserMixin, db.Model):
    __tablename__ = "User"

    id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    username = db.Column(db.String, index=True, unique=True, nullable=False)
    email = db.Column(db.String, index=True, unique=True, nullable=False)
    first_name = db.Column(db.String, nullable=True, default=None)
    last_name = db.Column(db.String, nullable=True, default=None)
    password_hash = db.Column(db.String)
    registr_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    last_access = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    approved = db.Column(db.Boolean, default=False)
    downloads = db.relationship("Downloads", backref="User", lazy="dynamic")

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def set_approval(self):
        self.approved = True
        db.session.commit()

    def unset_approval(self):
        self.approved = False
        db.session.commit()

    def update_last_access(self):
        self.last_access = datetime.utcnow()
        db.session.commit()

    def __repr__(self):
        return """User(id: {self.id}, username: {self.username})""".format(self=self)


class Downloads(db.Model):
    __tablename__ = "Downloads"

    id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
    dataset = db.Column(db.String, nullable=False)
    dl_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    user_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)

    def __repr__(self):
        return """Downloads(id: {self.id}, dataset: {self.dataset}, dl_date: {self.dl_date}, user_id: {self.user_id})""".format(self=self)


@login.user_loader
def load_user(id):
    return User.query.get(int(id))

