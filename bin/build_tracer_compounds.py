"""Builds the pinned tracer tables for MetabolicLabelingTracer (tracer_columns.py) from PubChem and RDKit.
Not part of the workflow run (the workflow environment has no RDKit): run by hand when the compound list
changes, and commit the two tables.

  bin/tracer_compounds.tsv       compound, names (';'-joined, lower case), pubchem_cid, smiles, and the number
                                 of atoms a uniform label covers: n_C, n_N, n_O, n_S, n_H (hydrogens on carbon,
                                 all hydrogens when there is none: they exchange elsewhere)
  bin/tracer_isotopologues.tsv   compound, labels ('13C;15N'), smiles: the compound with every such atom labelled

usage (an environment with rdkit and requests): build_tracer_compounds.py bin
"""
import csv, itertools, os, sys, time
import requests
from rdkit import Chem

# compound (PubChem name) -> other names and abbreviations used for it in tracer descriptions. Measured over the
# per-sample factors of Metabolomics Workbench and MetaboLights (2026-11-06): every compound named next to an
# isotope label in at least one study, polymers aside (cellulose, inulin: no single structure).
COMPOUNDS = {
    'D-glucose': ['glucose', 'glc', 'gluc', 'dextrose', 'd-glc'],
    'D-fructose': ['fructose', 'fru'],
    'D-galactose': ['galactose', 'gal'],
    'glycerol': ['glycerin'],
    'L-glutamine': ['glutamine', 'gln', 'l-gln'],
    'L-glutamic acid': ['glutamate', 'glutamic acid', 'l-glutamate', 'glu', 'l-glu'],
    'L-serine': ['serine', 'serin', 'ser', 'l-ser'],
    'glycine': ['gly'],
    'L-alanine': ['alanine', 'ala', 'l-ala'],
    'L-leucine': ['leucine', 'leu', 'l-leu'],
    'L-isoleucine': ['isoleucine', 'ile', 'l-ile'],
    'L-valine': ['valine', 'val', 'l-val'],
    'L-threonine': ['threonine', 'thr', 'l-thr'],
    'L-proline': ['proline', 'pro', 'l-pro'],
    'L-lysine': ['lysine', 'lys', 'l-lys'],
    'L-arginine': ['arginine', 'arg', 'l-arg'],
    'L-methionine': ['methionine', 'met', 'l-met'],
    'L-cysteine': ['cysteine', 'cys', 'l-cys'],
    'L-cystine': ['cystine'],
    'L-asparagine': ['asparagine', 'asn', 'l-asn'],
    'L-aspartic acid': ['aspartate', 'aspartic acid', 'l-aspartate', 'asp'],
    'L-phenylalanine': ['phenylalanine', 'phe', 'l-phe'],
    'L-tyrosine': ['tyrosine', 'tyr', 'l-tyr'],
    'L-tryptophan': ['tryptophan', 'trp', 'l-trp'],
    'L-histidine': ['histidine', 'his', 'l-his'],
    'pyruvic acid': ['pyruvate', 'pyr', 'sodium pyruvate'],
    'L-lactic acid': ['lactate', 'lactic acid', 'l-lactate', 'sodium lactate', 'lac'],
    'acetic acid': ['acetate', 'sodium acetate', 'ac'],
    'butyric acid': ['butyrate', 'sodium butyrate'],
    'propionic acid': ['propionate', 'sodium propionate'],
    'acetoacetic acid': ['acetoacetate'],
    'succinic acid': ['succinate'],
    'citric acid': ['citrate'],
    'palmitic acid': ['palmitate', 'bsa-palmitate'],
    'oleic acid': ['oleate', 'bsa-oleate'],
    'elaidic acid': ['elaidate', 'bsa-elaidate'],
    'choline': [],
    'L-carnitine': ['carnitine'],
    'ethanol': [],
    'methanol': [],
    'hypoxanthine': [],
    'adenine': [],
    'adenosine': [],
    'ammonium': ['ammonia', 'ammonium chloride', 'nh4cl', 'nh4'],
    'sodium bicarbonate': ['bicarbonate', 'nahco3'],
    'water': ['h2o'],
    'oxygen': ['o2', 'air', 'dioxygen'],
}
ISOTOPES = {'13C': ('C', 13), '15N': ('N', 15), '18O': ('O', 18), '34S': ('S', 34), '2H': ('H', 2)}


def pubchem(name):
    for attempt in range(4):
        try:
            r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}'
                             '/property/IsomericSMILES,SMILES/JSON', timeout=60)
            if r.status_code == 404:
                raise SystemExit(f'not in PubChem: {name}')
            r.raise_for_status()
            p = r.json()['PropertyTable']['Properties'][0]
            return p['CID'], p.get('IsomericSMILES') or p.get('SMILES')
        except (requests.RequestException, KeyError) as e:
            time.sleep(2 + 3 * attempt)
            err = e
    raise SystemExit(f'PubChem failed for {name}: {err}')


def label_targets(mol):
    """atoms each isotope labels: every atom of the element; hydrogens: those on carbon, or all when none is"""
    h = Chem.AddHs(mol)
    out = {}
    for iso, (el, _) in ISOTOPES.items():
        if el == 'H':
            on_c = [a.GetIdx() for a in h.GetAtoms() if a.GetSymbol() == 'H' and a.GetNeighbors()[0].GetSymbol() == 'C']
            out[iso] = on_c or [a.GetIdx() for a in h.GetAtoms() if a.GetSymbol() == 'H']
        else:
            out[iso] = [a.GetIdx() for a in h.GetAtoms() if a.GetSymbol() == el]
    return h, out


def labelled(h, targets, isos):
    m = Chem.RWMol(h)
    for iso in isos:
        for i in targets[iso]:
            m.GetAtomWithIdx(i).SetIsotope(ISOTOPES[iso][1])
    return Chem.MolToSmiles(Chem.RemoveHs(m.GetMol()))


def main():
    out = sys.argv[1]
    comp_rows, iso_rows, seen = [], [], {}
    for compound, names in COMPOUNDS.items():
        cid, smiles = pubchem(compound)
        mol = Chem.MolFromSmiles(smiles)
        h, targets = label_targets(mol)
        counts = {iso: len(t) for iso, t in targets.items()}
        all_names = sorted({compound.lower(), *[n.lower() for n in names]})
        for n in all_names:
            if n in seen:
                raise SystemExit(f'name {n} of {compound} also names {seen[n]}')
            seen[n] = compound
        comp_rows.append([compound, ';'.join(all_names), cid, Chem.MolToSmiles(mol),
                          counts['13C'], counts['15N'], counts['18O'], counts['34S'], counts['2H']])
        present = [iso for iso in ISOTOPES if counts[iso]]
        for k in range(1, len(present) + 1):
            for isos in itertools.combinations(present, k):
                iso_rows.append([compound, ';'.join(isos), labelled(h, targets, isos)])
        time.sleep(0.3)   # PubChem asks for at most 5 requests per second
    with open(os.path.join(out, 'tracer_compounds.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['compound', 'names', 'pubchem_cid', 'smiles', 'n_C', 'n_N', 'n_O', 'n_S', 'n_H'])
        w.writerows(comp_rows)
    with open(os.path.join(out, 'tracer_isotopologues.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['compound', 'labels', 'smiles'])
        w.writerows(iso_rows)
    print(len(comp_rows), 'compounds,', len(iso_rows), 'isotopologues')


if __name__ == '__main__':
    main()
