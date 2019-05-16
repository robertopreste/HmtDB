#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from .models import Locus, NtVariability, MitomapAa, MitomapDna, AaVariability, \
    Disease, Deletion, Insertion


def queryLocus(snpPosition):
    """Return the locus object corresponding to the specified position."""
    if snpPosition >= 16024 or snpPosition <= 576:
        return Locus.query.filter(Locus.locusType == "dloop").first()
    return Locus.query.filter(snpPosition >= Locus.startPosition,
                              snpPosition <= Locus.endPosition).first()


def queryNtVar_N(snp):
    """Return the nt variability object of the snp, for normal genome."""
    return NtVariability.query.filter(
        NtVariability.nucleotidePosition == snp.snpPosition,
        NtVariability.genomeType == "N").first()


def queryNtVar_P(snp):
    """Return the nt variability object of the snp, for pathologic genome."""
    return NtVariability.query.filter(
        NtVariability.nucleotidePosition == snp.snpPosition,
        NtVariability.genomeType == "P").first()


def queryMitomapAa(snp_aa_pos, locus):
    """Return the aa mitomap object for the corresponding locus and aa position."""
    return MitomapAa.query.filter(
        MitomapAa.geneName == locus.geneName,
        MitomapAa.aaPosition == int(snp_aa_pos)).first()


def queryMitomapDna(snp):
    """Return the dna mitomap object for the corresponding nucleotide position and snp type."""
    return MitomapDna.query.filter(
        MitomapDna.nucleotidePosition == snp.snpPosition,
        MitomapDna.rcrsType == snp.rcrsType,
        MitomapDna.snpType == snp.snpType).first()


def queryAaVar_N(snp_aa_pos, locus):
    """Return the aa variability object of the snp, for normal genome."""
    return AaVariability.query.filter(
        AaVariability.aaPos == snp_aa_pos,
        AaVariability.geneName == locus.geneName,
        AaVariability.genomeType == "N").first()


def queryAaVar_P(snp_aa_pos, locus):
    """Return the aa variability object of the snp, for pathologic genome."""
    return AaVariability.query.filter(
        AaVariability.aaPos == snp_aa_pos,
        AaVariability.geneName == locus.geneName,
        AaVariability.genomeType == "P").first()


def queryDisease(dId):
    """Return the disease object corresponding to the specified dId (diseaseId)."""
    return Disease.query.filter(Disease.diseaseId == dId).first()


def queryMitomapDnaDiseases(mito_dna):
    """Return the list of disease objects corresponding to the given MitomapDna entry."""
    if mito_dna is None or mito_dna.diseases is None:
        return None
    if ";" in mito_dna.diseases:
        dis = mito_dna.diseases.split(";")
        return [Disease.query.filter(Disease.diseaseId == int(el)).first()
                for el in dis]
    return [Disease.query.filter(Disease.diseaseId == int(mito_dna.diseases)).first()]


def queryDeletion(gId):
    """Return the list of deletions of the specified genome."""
    return Deletion.query.filter(Deletion.genomeId == gId).all()


def queryInsertion(gId):
    """Return the list of insertion of the specified genome."""
    return Insertion.query.filter(Insertion.genomeId == gId).all()


def getCodon(locus, snp_pos):
    if snp_pos % 3 == 0:  # codon position 3
        res = locus.rcrsDnaSeq[snp_pos - 3] \
              + locus.rcrsDnaSeq[snp_pos - 2] \
              + locus.rcrsDnaSeq[snp_pos - 1]
    elif snp_pos % 3 == 1:  # codon position 1
        res = locus.rcrsDnaSeq[snp_pos - 1] \
              + locus.rcrsDnaSeq[snp_pos] \
              + locus.rcrsDnaSeq[snp_pos + 1]
    else:  # codon position 2
        res = locus.rcrsDnaSeq[snp_pos - 2] \
              + locus.rcrsDnaSeq[snp_pos - 1] \
              + locus.rcrsDnaSeq[snp_pos]
    return res


def getAltCodon(locus, snp_pos, snp):
    # everything is shifted of 1 position due to Python's 0-based numbering
    if snp_pos % 3 == 0:
        res = locus.rcrsDnaSeq[snp_pos - 3] \
              + locus.rcrsDnaSeq[snp_pos - 2] \
              + snp.snpType
    elif snp_pos % 3 == 1:
        res = snp.snpType \
              + locus.rcrsDnaSeq[snp_pos] \
              + locus.rcrsDnaSeq[snp_pos + 1]
    else:
        res = locus.rcrsDnaSeq[snp_pos - 2] \
              + snp.snpType \
              + locus.rcrsDnaSeq[snp_pos]
    return res


aa_dict = {
    'AAA': 'K', 'CAA': 'Q', 'GAA': 'E', 'TAA': '*',
    'AAC': 'N', 'CAC': 'H', 'GAC': 'D', 'TAC': 'Y',
    'AAG': 'K', 'CAG': 'Q', 'GAG': 'E', 'TAG': '*',
    'AAT': 'N', 'CAT': 'H', 'GAT': 'D', 'TAT': 'Y',

    'ACA': 'T', 'CCA': 'P', 'GCA': 'A', 'TCA': 'S',
    'ACC': 'T', 'CCC': 'P', 'GCC': 'A', 'TCC': 'S',
    'ACG': 'T', 'CCG': 'P', 'GCG': 'A', 'TCG': 'S',
    'ACT': 'T', 'CCT': 'P', 'GCT': 'A', 'TCT': 'S',

    'AGA': '*', 'CGA': 'R', 'GGA': 'G', 'TGA': 'W',
    'AGC': 'S', 'CGC': 'R', 'GGC': 'G', 'TGC': 'C',
    'AGG': '*', 'CGG': 'R', 'GGG': 'G', 'TGG': 'W',
    'AGT': 'S', 'CGT': 'R', 'GGT': 'G', 'TGT': 'C',

    'ATA': 'M', 'CTA': 'L', 'GTA': 'V', 'TTA': 'L',
    'ATC': 'I', 'CTC': 'L', 'GTC': 'V', 'TTC': 'F',
    'ATG': 'M', 'CTG': 'L', 'GTG': 'V', 'TTG': 'L',
    'ATT': 'I', 'CTT': 'L', 'GTT': 'V', 'TTT': 'F',

    'AAR': 'K', 'CAR': 'Q', 'GAR': 'E', 'TAR': '*',
    'AAY': 'N', 'CAY': 'H', 'GAY': 'D', 'TAY': 'Y',
    'AAS': 'K/N', 'CAS': 'H/Q', 'GAS': 'D/E', 'TAS': 'Y/*',
    'AAW': 'K/N', 'CAW': 'H/Q', 'GAW': 'D/E', 'TAW': 'Y/*',
    'AAK': 'K/N', 'CAK': 'H/Q', 'GAK': 'D/E', 'TAK': 'Y/*',
    'AAM': 'K/N', 'CAM': 'H/Q', 'GAM': 'D/E', 'TAM': 'Y/*',
    'AAB': 'K/N', 'CAB': 'H/Q', 'GAB': 'D/E', 'TAB': 'Y/*',
    'AAD': 'K/N', 'CAD': 'H/Q', 'GAD': 'D/E', 'TAD': 'Y/*',
    'AAH': 'K/N', 'CAH': 'H/Q', 'GAH': 'D/E', 'TAH': 'Y/*',
    'AAV': 'K/N', 'CAV': 'H/Q', 'GAV': 'D/E', 'TAV': 'Y/*',
    'AAN': 'K/N', 'CAN': 'H/Q', 'GAN': 'D/E', 'TAN': 'Y/*',

    'ACR': 'T', 'CCR': 'P', 'GCR': 'A', 'TCR': 'S',
    'ACY': 'T', 'CCY': 'P', 'GCY': 'A', 'TCY': 'S',
    'ACS': 'T', 'CCS': 'P', 'GCS': 'A', 'TCS': 'S',
    'ACW': 'T', 'CCW': 'P', 'GCW': 'A', 'TCW': 'S',
    'ACK': 'T', 'CCK': 'P', 'GCK': 'A', 'TCK': 'S',
    'ACM': 'T', 'CCM': 'P', 'GCM': 'A', 'TCM': 'S',
    'ACB': 'T', 'CCB': 'P', 'GCB': 'A', 'TCB': 'S',
    'ACD': 'T', 'CCD': 'P', 'GCD': 'A', 'TCD': 'S',
    'ACH': 'T', 'CCH': 'P', 'GCH': 'A', 'TCH': 'S',
    'ACV': 'T', 'CCV': 'P', 'GCV': 'A', 'TCV': 'S',
    'ACN': 'T', 'CCN': 'P', 'GCN': 'A', 'TCN': 'S',

    'AGR': '*', 'CGR': 'R', 'GGR': 'G', 'TGR': 'W',
    'AGY': 'S', 'CGY': 'R', 'GGY': 'G', 'TGY': 'C',
    'AGS': 'S/*', 'CGS': 'R', 'GGS': 'G', 'TGS': 'C/W',
    'AGW': 'S/*', 'CGW': 'R', 'GGW': 'G', 'TGW': 'C/W',
    'AGK': 'S/*', 'CGK': 'R', 'GGK': 'G', 'TGK': 'C/W',
    'AGM': 'S/*', 'CGM': 'R', 'GGM': 'G', 'TGM': 'C/W',
    'AGB': 'S/*', 'CGB': 'R', 'GGB': 'G', 'TGB': 'C/W',
    'AGD': 'S/*', 'CGD': 'R', 'GGD': 'G', 'TGD': 'C/W',
    'AGH': 'S/*', 'CGH': 'R', 'GGH': 'G', 'TGH': 'C/W',
    'AGV': 'S/*', 'CGV': 'R', 'GGV': 'G', 'TGV': 'C/W',
    'AGN': 'S/*', 'CGN': 'R', 'GGN': 'G', 'TGN': 'C/W',

    'ATR': 'M', 'CTR': 'L', 'GTR': 'V', 'TTR': 'L',
    'ATY': 'I', 'CTY': 'L', 'GTY': 'V', 'TTY': 'F',
    'ATS': 'I/M', 'CTS': 'L', 'GTS': 'V', 'TTS': 'F/L',
    'ATW': 'I/M', 'CTW': 'L', 'GTW': 'V', 'TTW': 'F/L',
    'ATK': 'I/M', 'CTK': 'L', 'GTK': 'V', 'TTK': 'F/L',
    'ATM': 'I/M', 'CTM': 'L', 'GTM': 'V', 'TTM': 'F/L',
    'ATB': 'I/M', 'CTB': 'L', 'GTB': 'V', 'TTB': 'F/L',
    'ATD': 'I/M', 'CTD': 'L', 'GTD': 'V', 'TTD': 'F/L',
    'ATH': 'I/M', 'CTH': 'L', 'GTH': 'V', 'TTH': 'F/L',
    'ATV': 'I/M', 'CTV': 'L', 'GTV': 'V', 'TTV': 'F/L',
    'ATN': 'I/M', 'CTN': 'L', 'GTN': 'V', 'TTN': 'F/L',

    'ARA': 'K/*', 'CRA': 'Q/R', 'GRA': 'E/G', 'TRA': 'W/*',
    'AYA': 'M/T', 'CYA': 'L/P', 'GYA': 'A/V', 'TYA': 'L/S',
    'ASA': 'T/*', 'CSA': 'P/R', 'GSA': 'A/G', 'TSA': 'S/W',
    'AWA': 'K/M', 'CWA': 'L/Q', 'GWA': 'E/V', 'TWA': 'L/*',
    'AKA': 'M/*', 'CKA': 'L/R', 'GKA': 'G/V', 'TKA': 'L/W',
    'AMA': 'K/T', 'CMA': 'P/Q', 'GMA': 'A/E', 'TMA': 'S/*',

    'ARC': 'N/S', 'CRC': 'H/R', 'GRC': 'D/G', 'TRC': 'C/Y',
    'AYC': 'I/T', 'CYC': 'L/P', 'GYC': 'A/V', 'TYC': 'F/S',
    'ASC': 'S/T', 'CSC': 'P/R', 'GSC': 'A/G', 'TSC': 'C/S',
    'AWC': 'I/N', 'CWC': 'H/L', 'GWC': 'D/V', 'TWC': 'F/Y',
    'AKC': 'I/S', 'CKC': 'L/R', 'GKC': 'G/V', 'TKC': 'C/F',
    'AMC': 'N/T', 'CMC': 'H/P', 'GMC': 'A/D', 'TMC': 'S/Y',

    'ARG': 'K/*', 'CRG': 'Q/R', 'GRG': 'E/G', 'TRG': 'W/*',
    'AYG': 'M/T', 'CYG': 'L/P', 'GYG': 'A/V', 'TYG': 'L/S',
    'ASG': 'T/*', 'CSG': 'P/R', 'GSG': 'A/G', 'TSG': 'S/W',
    'AWG': 'K/M', 'CWG': 'L/Q', 'GWG': 'E/V', 'TWG': 'L/*',
    'AKG': 'M/*', 'CKG': 'L/R', 'GKG': 'G/V', 'TKG': 'L/W',
    'AMG': 'K/T', 'CMG': 'P/Q', 'GMG': 'A/E', 'TMG': 'S/*',

    'ART': 'N/S', 'CRT': 'H/R', 'GRT': 'D/G', 'TRT': 'C/Y',
    'AYT': 'I/T', 'CYT': 'L/P', 'GYT': 'A/V', 'TYT': 'F/S',
    'AST': 'S/T', 'CST': 'P/R', 'GST': 'A/G', 'TST': 'C/S',
    'AWT': 'I/N', 'CWT': 'H/L', 'GWT': 'D/V', 'TWT': 'F/Y',
    'AKT': 'I/S', 'CKT': 'L/R', 'GKT': 'G/V', 'TKT': 'C/F',
    'AMT': 'N/T', 'CMT': 'H/P', 'GMT': 'A/D', 'TMT': 'S/Y',

    'ARR': 'K/*', 'CRR': 'Q/R', 'GRR': 'E/G', 'TRR': 'W/*',
    'AYR': 'M/T', 'CYR': 'L/P', 'GYR': 'A/V', 'TYR': 'L/S',
    'ASR': 'T/*', 'CSR': 'P/R', 'GSR': 'A/G', 'TSR': 'S/W',
    'AWR': 'K/M', 'CWR': 'L/Q', 'GWR': 'E/V', 'TWR': 'L/*',
    'AKR': 'M/*', 'CKR': 'L/R', 'GKR': 'G/V', 'TKR': 'L/W',
    'AMR': 'K/T', 'CMR': 'P/Q', 'GMR': 'A/E', 'TMR': 'S/*',
    'ARY': 'N/S', 'CRY': 'H/R', 'GRY': 'D/G', 'TRY': 'C/Y',
    'AYY': 'I/T', 'CYY': 'L/P', 'GYY': 'A/V', 'TYY': 'F/S',
    'ASY': 'S/T', 'CSY': 'P/R', 'GSY': 'A/G', 'TSY': 'C/S',
    'AWY': 'I/N', 'CWY': 'H/L', 'GWY': 'D/V', 'TWY': 'F/Y',
    'AKY': 'I/S', 'CKY': 'L/R', 'GKY': 'G/V', 'TKY': 'C/F',
    'AMY': 'N/T', 'CMY': 'H/P', 'GMY': 'A/D', 'TMY': 'S/Y',

    'CYS': 'L/P', 'GYS': 'A/V',
    'CYW': 'L/P', 'GYW': 'A/V',
    'CYK': 'L/P', 'GYK': 'A/V',
    'CYM': 'L/P', 'GYM': 'A/V',
    'CYB': 'L/P', 'GYB': 'A/V',
    'CYD': 'L/P', 'GYD': 'A/V',
    'CYH': 'L/P', 'GYH': 'A/V',
    'CYV': 'L/P', 'GYV': 'A/V',
    'CYN': 'L/P', 'GYN': 'A/V',

    'CSS': 'P/R', 'GSS': 'A/G',
    'CSW': 'P/R', 'GSW': 'A/G',
    'CSK': 'P/R', 'GSK': 'A/G',
    'CSM': 'P/R', 'GSM': 'A/G',
    'CSB': 'P/R', 'GSB': 'A/G',
    'CSD': 'P/R', 'GSD': 'A/G',
    'CSH': 'P/R', 'GSH': 'A/G',
    'CSV': 'P/R', 'GSV': 'A/G',
    'CSN': 'P/R', 'GSN': 'A/G',

    'CKS': 'L/R', 'GKS': 'G/V',
    'CKW': 'L/R', 'GKW': 'G/V',
    'CKK': 'L/R', 'GKK': 'G/V',
    'CKM': 'L/R', 'GKM': 'G/V',
    'CKB': 'L/R', 'GKB': 'G/V',
    'CKD': 'L/R', 'GKD': 'G/V',
    'CKH': 'L/R', 'GKH': 'G/V',
    'CKV': 'L/R', 'GKV': 'G/V',
    'CKN': 'L/R', 'GKN': 'G/V',

    'YTR': 'L',
    'YTA': 'L',
    'YTG': 'L',

}


def getAa(codon):
    try:
        return aa_dict[codon]
    except KeyError:
        return 'X'

