#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app import db


"""
 - Locus 1 --- M AaVariability
 - Country 1 -- M IndividualsData
 - Genome 1 -- M Deletion
 - Disease 1 -- M IndividualsData
 - Genome 1 -- M GenAlignment
 - Genome 1 -- M GenAnnotation
 - IndividualsData 1 -- M Genome
 - Methods 1 -- M Genome
 - Sources 1 -- M Genome
 - Genome 1 -- M References
 - Genome 1 -- M GenomeSnp
 - Genome 1 -- M Insertion
 - Disease 1 -- M MitomapAa and Disease 1 -- M MitomapDna
"""

deletion_genome = db.Table("deletion_genome",
                           db.Column("deletion_id", db.Integer,
                                     db.ForeignKey("Deletion.deletionId")),
                           db.Column("genome_id", db.Integer,
                                     db.ForeignKey("Genome.genomeId")))

disease_individuals = db.Table("disease_individuals",
                               db.Column("disease_id", db.Integer,
                                         db.ForeignKey("Disease.diseaseId")),
                               db.Column("individual_id", db.Integer,
                                         db.ForeignKey("IndividualsData.individualId")))

genome_reference = db.Table("genome_reference",
                            db.Column("genome_id", db.Integer,
                                      db.ForeignKey("Genome.genomeId")),
                            db.Column("reference_id", db.Integer,
                                      db.ForeignKey("Reference.referenceId")))

genome_snps = db.Table("genome_snps",
                       db.Column("genome_id", db.Integer,
                                 db.ForeignKey("Genome.genomeId")),
                       db.Column("snps_id", db.Integer,
                                 db.ForeignKey("GenomeSnp.snpId")))


class AaVariability(db.Model):
    __tablename__ = "AaVariability"

    aaId = db.Column(db.Integer, nullable=False, autoincrement=True,
                     primary_key=True)
    aaPos = db.Column(db.Integer, nullable=False)
    rcrsAa = db.Column(db.String(16), nullable=True, default=None)
    varAa_intrahs = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_intermam = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_eu = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_am = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_af = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_as = db.Column(db.Float, nullable=True, default=0.00000)
    varAa_oc = db.Column(db.Float, nullable=True, default=0.00000)
    genomeType = db.Column(db.String(16), nullable=False, default="N")
    # relationship with Locus - "many" side
    geneName = db.Column(db.String(16), db.ForeignKey("Locus.geneName"),
                         nullable=False)

    def __repr__(self):
        return """AaVariability(aaId: {self.aaId}, aaPos: {self.aaPos}, rcrsAa: {self.rcrsAa}, 
        varAa_intrahs: {self.varAa_intrahs}, varAa_intermam: {self.varAa_intermam}, 
        varAa_eu: {self.varAa_eu}, varAa_am: {self.varAa_am}, varAa_af: {self.varAa_af}, 
        varAa_as: {self.varAa_as}, varAa_oc: {self.varAa_oc}, genomeType: {self.genomeType}, 
        geneName: {self.geneName})""".format(self=self)


class Blosum(db.Model):
    __tablename__ = "Blosum"

    blosumId = db.Column(db.Integer, nullable=False, autoincrement=True,
                         primary_key=True)
    aaChange = db.Column(db.String(16), nullable=False)
    blosumIndex = db.Column(db.Numeric(6, 3), nullable=False)

    def __repr__(self):
        return """Blosum(blosumId: {self.blosumId}, aaChange: {self.aaChange}, 
        blosumIndex: {self.blosumIndex})""".format(self=self)


class Country(db.Model):
    __tablename__ = "Country"

    countryId = db.Column(db.Integer, nullable=False, autoincrement=True,
                          primary_key=True)
    countryName = db.Column(db.String(64), nullable=False)
    countryCode = db.Column(db.String(16), nullable=False, default="XX")
    continentName = db.Column(db.String(32), nullable=False)
    continentCode = db.Column(db.String(16), nullable=False, default="XX")
    # relationship with IndividualsData - "one" side
    individuals = db.relationship("IndividualsData", backref="Country",
                                  lazy="dynamic")

    def __repr__(self):
        return """Country(countryId: {self.countryId}, countryName: {self.countryName}, 
        countryCode: {self.countryCode}, continentName: {self.continentName}, 
        continentCode: {self.continentCode}, individuals: {self.individuals})""".format(self=self)


class Deletion(db.Model):
    __tablename__ = "Deletion"

    deletionId = db.Column(db.Integer, nullable=False, autoincrement=True,
                           primary_key=True)
    fromPosition = db.Column(db.Integer, nullable=False)
    toPosition = db.Column(db.Integer, nullable=False)
    # relationship with Genome (M -- M)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """Deletion(deletionId: {self.deletionId}, fromPosition: {self.fromPosition}, 
        toPosition: {self.toPosition}, genomeId: {self.genomeId})""".format(self=self)


class Disease(db.Model):
    __tablename__ = "Disease"

    diseaseId = db.Column(db.Integer, nullable=False, autoincrement=True,
                          primary_key=True)
    diseaseName = db.Column(db.String(200), nullable=False)
    diseaseAcronym = db.Column(db.String(64), nullable=True, default=None)
    url = db.Column(db.String(200), nullable=True, default=None)
    # relationship with IndividualsData (M -- M)
    individuals = db.relationship("IndividualsData", backref="Disease",
                                  lazy="dynamic")
    # relationships with MitomapAa and MitomapDna - one side
    # mitoAa = db.relationship("MitomapAa", backref="Disease", lazy="dynamic")
    # mitoDna = db.relationship("MitomapDna", backref="Disease", lazy="dynamic")
    # TODO: fix diseases related to MitomapDna

    def __repr__(self):
        return """Disease(diseaseId: {self.diseaseId}, diseaseName: {self.diseaseName}, 
        diseaseAcronym: {self.diseaseAcronym}, url: {self.url}, individuals: {self.individuals})""".format(self=self)


class EthnicGroups(db.Model):
    __tablename__ = "EthnicGroups"

    groupId = db.Column(db.Integer, nullable=False, autoincrement=True,
                        primary_key=True)
    groupName = db.Column(db.String(200), nullable=True, default=None)
    groupDescription = db.Column(db.String(100), nullable=True, default=None)
    # relationship with IndividualsData - "one" side
    individuals = db.relationship("IndividualsData", backref="EthnicGroups",
                                  lazy="dynamic")

    def __repr__(self):
        return """EthnicGroups(groupId: {self.groupId}, groupName: {self.groupName}, 
        groupDescription: {self.groupDescription}, individuals: {self.individuals})""".format(self=self)


class GenAlignment(db.Model):
    __tablename__ = "GenAlignment"

    alignmentId = db.Column(db.Integer, nullable=False, autoincrement=True,
                            primary_key=True)
    alignment = db.Column(db.Text, nullable=False)
    # relationship with Genome - (one to one)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """GenAlignment(alignmentId: {self.alignmentId}, alignment: {self.alignment}, 
        genomeId: {self.genomeId})""".format(self=self)


class GenAnnotation(db.Model):
    __tablename__ = "GenAnnotation"

    annotationId = db.Column(db.Integer, nullable=False, autoincrement=True,
                             primary_key=True)
    annotations = db.Column(db.String(64), nullable=False)
    # relationship with Genome - (one to one)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """GenAnnotation(annotationId: {self.annotationId}, annotations: {self.annotations}, 
        genomeId: {self.genomeId})""".format(self=self)


class Genome(db.Model):
    __tablename__ = "Genome"

    genomeId = db.Column(db.Integer, nullable=False, autoincrement=True,
                         primary_key=True)
    genomeSequence = db.Column(db.Text, nullable=False)
    completeGenome = db.Column(db.String(16), nullable=False)
    startPosition = db.Column(db.Integer, nullable=True, default=None)
    endPosition = db.Column(db.Integer, nullable=True, default=None)
    haplotypeUser = db.Column(db.String(200), nullable=True, default=None)
    haplotypeHmdb = db.Column(db.String(200), nullable=True, default=None)
    haplogroupUser = db.Column(db.String(200), nullable=True, default=None)
    haplogroupHmdb = db.Column(db.String(200), nullable=True, default=None)
    referenceDb = db.Column(db.String(100), nullable=True, default=None)
    referenceDbId = db.Column(db.String(100), nullable=True, default=None)
    genomeType = db.Column(db.String(16), nullable=False, default="N")
    sourceDbId = db.Column(db.String(100), nullable=True, default=None)
    # relationship with IndividualsData - (one to one)
    individualId = db.Column(db.Integer,
                             db.ForeignKey("IndividualsData.individualId"),
                             nullable=True, default=None)
    # relationship with GenAlignment - (one to one)
    alignmentId = db.relationship("GenAlignment", uselist=False,
                                  backref="Genome")
    # relationship with GenAnnotation - (one to one)
    annotationId = db.relationship("GenAnnotation", uselist=False,
                                   backref="Genome")
    # relationship with Methods - "many" side
    methodId = db.Column(db.Integer, db.ForeignKey("Methods.methodId"),
                         nullable=True, default=None)
    # relationship with Sources - "many" side
    sourceId = db.Column(db.Integer, db.ForeignKey("Sources.sourceId"),
                         nullable=True, default=None)
    # relationship with References - "one" side
    references = db.relationship("Reference", backref="Genome", lazy="dynamic")
    # relationship with GenomeSnp - "one" side
    snps = db.relationship("GenomeSnp", backref="Genome", lazy="dynamic")
    # relationship with Deletion (M -- M) (ora 1)
    deletions = db.relationship("Deletion", backref="Genome", lazy="dynamic")
    # relationship with Insertion - ora 1
    insertions = db.relationship("Insertion", backref="Genome", lazy="dynamic")

    def __repr__(self):
        return """Genome(genomeId: {self.genomeId}, genomeSequence: {self.genomeSequence}, 
        completeGenome: {self.completeGenome}, startPosition: {self.startPosition}, 
        endPosition: {self.endPosition}, haplotypeUser: {self.haplotypeUser}, 
        haplotypeHmdb: {self.haplotypeHmdb}, haplogroupUser: {self.haplogroupUser}, 
        haplogroupHmdb: {self.haplogroupUser}, referenceDb: {self.referenceDb}, 
        referenceDbId: {self.referenceDbId}, genomeType: {self.genomeType}, 
        sourceDbId: {self.sourceDbId}, individualId: {self.individualId}, 
        alignmentId: {self.alignmentId}, annotationId: {self.annotationId}, 
        methodId: {self.methodId}, sourceId: {self.sourceId}, references: {self.references}, 
        snps: {self.snps}, deletions: {self.deletions}, insertions: {self.insertions})""".format(self=self)


class GenomeSnp(db.Model):
    __tablename__ = "GenomeSnp"

    snpId = db.Column(db.Integer, nullable=False, autoincrement=True,
                      primary_key=True)
    snpPosition = db.Column(db.Integer, nullable=False)
    snpType = db.Column(db.String(64), nullable=False)
    rcrsType = db.Column(db.String(64), nullable=False, default="")
    # relationship with Genome - "many" side
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """GenomeSnp(snpId: {self.snpId}, snpPosition: {self.snpPosition}, 
        snpType: {self.snpType}, rcrsType: {self.rcrsType}, genomeId: {self.genomeId})""".format(self=self)


class IndividualsData(db.Model):
    __tablename__ = "IndividualsData"

    individualId = db.Column(db.Integer, nullable=False, autoincrement=True,
                             primary_key=True)
    age = db.Column(db.Integer, nullable=True, default=None)
    sex = db.Column(db.String(16), nullable=True, default=None)
    individualType = db.Column(db.String(16), nullable=True, default=None)
    # relationship with Genome - (one to one)
    genomeId = db.relationship("Genome", backref="IndividualsData",
                               lazy="dynamic")
    # relationship with Country - "many" side
    countryId = db.Column(db.Integer, db.ForeignKey("Country.countryId"),
                          nullable=True, default=None)
    # relationship with EthnicGroups - "many" side
    groupId = db.Column(db.Integer, db.ForeignKey("EthnicGroups.groupId"),
                        nullable=True, default=None)
    # relationship with Disease (M -- M)
    diseases = db.Column(db.Integer, db.ForeignKey("Disease.diseaseId"),
                         nullable=True, default=None)

    def __repr__(self):
        return """IndividualsData(individualId: {self.individualId}, age: {self.age}, 
        sex: {self.sex}, individualType: {self.individualType}, genomeId: {self.genomeId}, 
        countryId: {self.countryId}, groupId: {self.groupId}, diseases: {self.diseases})""".format(self=self)


class Insertion(db.Model):
    __tablename__ = "Insertion"

    insertionId = db.Column(db.Integer, nullable=False, autoincrement=True,
                            primary_key=True)
    position5P = db.Column(db.Integer, nullable=False)
    sequence = db.Column(db.String(64), nullable=True)
    # relationship with Genome - 1 side
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """Insertion(insertionId: {self.insertionId}, position5P: {self.position5P}, 
        sequence: {self.sequence}, genomeId: {self.genomeId})""".format(self=self)


class Locus(db.Model):
    __tablename__ = "Locus"

    locusId = db.Column(db.Integer, nullable=False, autoincrement=True,
                        primary_key=True)
    geneName = db.Column(db.String(64), nullable=False, default="", index=True)
    locusType = db.Column(db.String(64), nullable=False)
    startPosition = db.Column(db.Integer, nullable=False)
    endPosition = db.Column(db.Integer, nullable=False)
    startPattern = db.Column(db.String(2000), nullable=True, default=None)
    endPattern = db.Column(db.String(2000), nullable=True, default=None)
    rcrsAaSeq = db.Column(db.String(4000), nullable=True, default=None)
    rcrsDnaSeq = db.Column(db.String(4000), nullable=True, default=None)
    # relationship with MitomapAa - "one" side
    mitomapAa = db.relationship("MitomapAa", backref="Locus", lazy="dynamic")
    # relationship with MitomapDna - "one" side
    mitomapDna = db.relationship("MitomapDna", backref="Locus", lazy="dynamic")
    # relationship with AaVariability - "one" side
    aaVariability = db.relationship("AaVariability", backref="Locus",
                                    lazy="dynamic")
    # relationship with NtVariability - "one" side
    ntVariability = db.relationship("NtVariability", backref="Locus",
                                    lazy="dynamic")

    def __repr__(self):
        return """Locus(locusId: {self.locusId}, geneName: {self.geneName}, 
        locusType: {self.locusType}, startPosition: {self.startPosition}, 
        endPosition: {self.endPosition}, startPattern: {self.startPattern}, 
        endPattern: {self.endPattern}, rcrsAaSeq: {self.rcrsAaSeq}, mitomapAa: {self.mitomapAa}, 
        mitomapDna: {self.mitomapDna}, aaVariability: {self.aaVariability}, 
        ntVariability: {self.ntVariability})""".format(self=self)


class Methods(db.Model):
    __tablename__ = "Methods"

    methodId = db.Column(db.Integer, nullable=False, autoincrement=True,
                         primary_key=True)
    methodName = db.Column(db.String(200), nullable=False)
    methodDescription = db.Column(db.String(2000))
    # relationship with Genome - "one" side
    genomes = db.relationship("Genome", backref="Methods", lazy="dynamic")

    def __repr__(self):
        return """Methods(methodId: {self.methodId}, methodName: {self.methodName}, 
        methodDescription: {self.methodDescription}, genomes: {self.genomes})""".format(self=self)


class MitomapAa(db.Model):
    __tablename__ = "MitomapAa"

    varAaId = db.Column(db.Integer, nullable=False, autoincrement=True,
                        primary_key=True)
    aaChange = db.Column(db.String(16), nullable=False, default="SYNON")
    aaPosition = db.Column(db.Integer, nullable=False)
    variationType = db.Column(db.Integer, nullable=True, default=None)
    # relationship with Locus - "many" side
    geneName = db.Column(db.String(32), db.ForeignKey("Locus.geneName"),
                         nullable=False)
    # relationship with MitomapDna - "many" side
    varDnaId1 = db.Column(db.Integer, nullable=False)
    varDnaId2 = db.Column(db.Integer, nullable=True, default=None)
    varDnaId3 = db.Column(db.Integer, nullable=True, default=None)
    diseases = db.Column(db.String(200), nullable=True, default=None)
    # relationship with Disease
    # diseases = db.Column(db.Integer, db.ForeignKey("Disease.diseaseId"),
    #                      nullable=True, default=None)

    def __repr__(self):
        return """MitomapAa(varAaId: {self.varAaId}, aaChange: {self.aaChange}, 
        aaPosition: {self.aaPosition}, variationType: {self.variationType}, 
        geneName: {self.geneName}, varDnaId1: {self.varDnaId1}, varDnaId2: {self.varDnaId2}, 
        varDnaId3: {self.varDnaId3})""".format(self=self)


class MitomapDna(db.Model):
    __tablename__ = "MitomapDna"

    varDnaId = db.Column(db.Integer, nullable=False, autoincrement=True,
                         primary_key=True)
    nucleotidePosition = db.Column(db.Integer, nullable=False)
    snpType = db.Column(db.String(64), nullable=False)
    rcrsType = db.Column(db.String(64), nullable=False, default="")
    # relationship with Locus - "many" side
    geneName = db.Column(db.String(16), db.ForeignKey("Locus.geneName"),
                         nullable=True, default=None)
    diseases = db.Column(db.String(200), nullable=True, default=None)
    # relationship with Disease
    # diseases = db.Column(db.String(200), db.ForeignKey("Disease.diseaseId"),
    #                      nullable=True, default=None)

    def __repr__(self):
        return """MitomapDna(varDnaId: {self.varDnaId}, 
        nucleotidePosition: {self.nucleotidePosition}, snpType: {self.snpType}, 
        rcrsType: {self.rcrsType}, geneName: {self.geneName})""".format(self=self)


class NtVariability(db.Model):
    __tablename__ = "NtVariability"

    ntId = db.Column(db.Integer, nullable=False, autoincrement=True,
                     primary_key=True)
    nucleotidePosition = db.Column(db.Integer, nullable=False)
    insertionPosition = db.Column(db.Integer, nullable=False)
    var_tot = db.Column(db.Float, nullable=True, default=0.00000)
    var_eu = db.Column(db.Float, nullable=True, default=0.00000)
    var_am = db.Column(db.Float, nullable=True, default=0.00000)
    var_af = db.Column(db.Float, nullable=True, default=0.00000)
    var_as = db.Column(db.Float, nullable=True, default=0.00000)
    var_oc = db.Column(db.Float, nullable=True, default=0.00000)
    compVar_tot = db.Column(db.String(200), nullable=True, default=None)
    compVar_eu = db.Column(db.String(200), nullable=True, default=None)
    compVar_am = db.Column(db.String(200), nullable=True, default=None)
    compVar_af = db.Column(db.String(200), nullable=True, default=None)
    compVar_as = db.Column(db.String(200), nullable=True, default=None)
    compVar_oc = db.Column(db.String(200), nullable=True, default=None)
    genomeType = db.Column(db.String(100), nullable=False, default="N")
    # relationship with Locus - "many" side
    geneName = db.Column(db.String(16), db.ForeignKey("Locus.geneName"),
                         nullable=True, default=None)

    def __repr__(self):
        return """NtVariability(ntId: {self.ntId}, nucleotidePosition: {self.nucleotidePosition}, 
        insertionPosition: {self.insertionPosition}, var_tot: {self.var_tot}, var_eu: {self.var_eu}, 
        var_am: {self.var_am}, var_af: {self.var_af}, var_as: {self.var_as}, var_oc: {self.var_oc}, 
        compVar_tot: {self.compVar_tot}, compVar_eu: {self.compVar_eu}, 
        compVar_am: {self.compVar_am}, compVar_af: {self.compVar_af}, compVar_as: {self.compVar_as}, 
        compVar_oc: {self.compVar_oc}, genomeType: {self.genomeType}, 
        geneName: {self.geneName})""".format(self=self)


class Reference(db.Model):
    __tablename__ = "Reference"

    referenceId = db.Column(db.Integer, nullable=False, autoincrement=True,
                            primary_key=True)
    pubmedId = db.Column(db.String(32), nullable=True, default=None)
    author = db.Column(db.String(500), nullable=True, default=None)
    title = db.Column(db.String(2000), nullable=True, default=None)
    volume = db.Column(db.String(64), nullable=True, default=None)
    year = db.Column(db.String(32), nullable=True, default=None)
    firstPage = db.Column(db.String(32), nullable=True, default=None)
    paper = db.Column(db.String(64), nullable=True, default=None)
    issue = db.Column(db.String(32), nullable=True, default=None)
    genomeId = db.Column(db.Integer, db.ForeignKey("Genome.genomeId"),
                         nullable=False)

    def __repr__(self):
        return """Reference(referenceId: {self.referenceId}, pubmedId: {self.pubmedId}, 
        author: {self.author}, title: {self.title}, volume: {self.volume}, year: {self.year}, 
        firstPage: {self.firstPage}, paper: {self.paper}, issue: {self.issue}, 
        genomeId: {self.genomeId})""".format(self=self)


class Sources(db.Model):
    __tablename__ = "Sources"

    sourceId = db.Column(db.Integer, nullable=False, autoincrement=True,
                         primary_key=True)
    sourceName = db.Column(db.String(100), nullable=False)
    # relationship with Genome - "one" side
    genomes = db.relationship("Genome", backref="Sources", lazy="dynamic")

    def __repr__(self):
        return """Sources(sourceId: {self.sourceId}, sourceName: {self.sourceName}, 
        genomes: {self.genomes})""".format(self=self)


class Stats(db.Model):
    __tablename__ = "Stats"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    continentName = db.Column(db.String(64), nullable=False)
    completeGenome = db.Column(db.String(16), nullable=False)
    genomeType = db.Column(db.String(16), nullable=False, default="N")
    total = db.Column(db.Integer, nullable=False, default="0")

    def __repr__(self):
        return """Stats(id: {self.id}, continentName: {self.continentName}, 
        completeGenome: {self.completeGenome}, genomeType: {self.genomeType}, 
        total: {self.total})""".format(self=self)

