import numpy as np

import mummi_ras.online.cg.fast_rdf2d
from mummi_ras.datastructures.rdf import RDF

def RAS_lipid_RDF(syst,bt):
    lipidnames = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    range_max = RDF.RNG_MAX #120  in Angstroms
    bin_size = RDF.SZ_BIN #0.2
    range_min = RDF.RNG_MIN #0.1
    bins = RDF.NBINS #int((range_max-range_min)/bin_size)

    REF = syst.atoms[[syst.select_atoms('name F1')[0].index]]

    values = np.zeros([len(lipidnames)+1,bins])
    for i in range(len(lipidnames)):
        SEL = bt.residues.atoms.select_atoms('resname '+str(lipidnames[i])+' and name C1A D1A T1A R1')
        kras_rdf = mummi_ras.online.cg.fast_rdf2d.InterRDF(REF, SEL, nbins=bins, range=(range_min,range_max))
        kras_rdf.run()
        values[i+1] = kras_rdf.count
        values[0] = kras_rdf.bins
    return values

def RAS4A_lipid_RDF(syst,bt):
    lipidnames = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    range_max = RDF.RNG_MAX #120  in Angstroms
    bin_size = RDF.SZ_BIN #0.2
    range_min = RDF.RNG_MIN #0.1
    bins = RDF.NBINS #int((range_max-range_min)/bin_size)

    F1_resid = syst.select_atoms('name F1')[0].resid
    P1_resid = syst.select_atoms('name P1')[0].resid
    REF = syst.select_atoms('(resid '+str(F1_resid)+' and name F1) or (resid '+str(P1_resid)+' and name P1)')
    #REF = syst.atoms[[syst.select_atoms('name F1')[0].index]]

    values = np.zeros([len(lipidnames)+1,bins])
    for i in range(len(lipidnames)):
        SEL = bt.residues.atoms.select_atoms('resname '+str(lipidnames[i])+' and name C1A D1A T1A R1')
        kras_rdf = mummi_ras.online.cg.fast_rdf2d.InterRDF(REF, SEL, nbins=bins, range=(range_min,range_max))
        kras_rdf.run()
        values[i+1] = kras_rdf.count
        values[0] = kras_rdf.bins
    return values

def RAF_lipid_RDF(syst,bt):
    lipidnames = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    range_max = RDF.RNG_MAX #120  in Angstroms
    bin_size = RDF.SZ_BIN #0.2
    range_min = RDF.RNG_MIN #0.1
    bins = RDF.NBINS #int((range_max-range_min)/bin_size)

    start_resid = syst.select_atoms('name F1')[0].resid
    REF = syst.select_atoms('resid '+str(start_resid+95)+':'+str(start_resid+98)+' '+str(start_resid+108)+':'+str(start_resid+111)+' and name BB')

    values = np.zeros([len(lipidnames)+1,bins])
    for i in range(len(lipidnames)):
        SEL = bt.residues.atoms.select_atoms('resname '+str(lipidnames[i])+' and name C1A D1A T1A R1')
        kras_rdf = mummi_ras.online.cg.fast_rdf2d.InterRDF(REF, SEL, nbins=bins, range=(range_min,range_max))
        kras_rdf.run()
        values[i+1] = kras_rdf.count
        values[0] = kras_rdf.bins
    return values


def KRAS_lipid_RDF_rdf(syst,bt):
    lipidnames = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    range_max = RDF.RNG_MAX #70  in Angstroms
    bin_size = RDF.SZ_BIN #0.2
    range_min = RDF.RNG_MIN #0.1
    bins = RDF.NBINS #int((range_max-range_min)/bin_size)

    REF = syst.atoms[[syst.select_atoms('name F1')[0].index]]

    values = np.zeros([len(lipidnames)+1,bins])
    for i in range(len(lipidnames)):
        SEL = bt.residues.atoms.select_atoms('resname '+str(lipidnames[i])+' and name C1A D1A T1A R1')
        kras_rdf = mummi_ras.online.cg.fast_rdf2d.InterRDF(REF, SEL, nbins=bins, range=(range_min,range_max))
        kras_rdf.run()
        values[i+1] = kras_rdf.count
        values[0] = kras_rdf.bins
    return values
