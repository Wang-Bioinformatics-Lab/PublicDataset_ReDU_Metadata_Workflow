"""Fill missing MassSpectrometer values of GNPS/MassIVE rows from the MassIVE dataset
description (PROXI: https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/<MSV>).

Only rows whose MassSpectrometer is 'missing value' are touched, and only when the
dataset lists exactly one instrument that maps unambiguously to an allowed
MassSpectrometer term. Submitter values are never overwritten. The instrument is a
dataset-level statement, so it is not used for anything sample-specific (species are
deliberately not filled this way: blanks and QCs would get the organism).

If anything goes wrong (network, unexpected input), the input is written unchanged so
the workflow is never blocked by this optional enrichment.
"""
import argparse
import json
import re
import shutil
import sys
import time

import pandas as pd
import requests

from ms_instrument_matching import match_mass_spectrometer

PROXI_URL = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{}"
MISSING = "missing value"
# PROXI's generic placeholder when the submitter did not pick a model
PLACEHOLDER_ACCESSIONS = {"MS:1000031"}


def fetch_instruments(session, accession, retries=3, timeout=(10, 60)):
    """Returns the PROXI instrument list of a dataset, or None if it could not be fetched."""
    for attempt in range(1, retries + 1):
        try:
            r = session.get(PROXI_URL.format(accession), timeout=timeout)
            if r.status_code == 200:
                # some dataset descriptions are not valid UTF-8 (e.g. latin-1 degree signs)
                return json.loads(r.content.decode("utf-8", errors="replace")).get("instruments") or []
            if r.status_code < 500 and r.status_code != 429:
                return None
        except (requests.RequestException, ValueError) as e:
            print(f"Attempt {attempt} failed for {accession}: {e}")
        time.sleep(min(2 ** attempt, 30))
    return None


def _is_placeholder(instrument):
    return instrument.get("accession") in PLACEHOLDER_ACCESSIONS \
        or str(instrument.get("name", "")).strip().lower() == "instrument model"


def free_text_instrument(value, allowed_ms_terms):
    """Allowed term for the submitter's free text of a generic 'instrument model' entry,
    e.g. 'Q Exactive|1001911' (MS accession number after a '|'), 'MS:1002634' or a name
    such as 'Exactive (Thermo Scientific)'; None unless it maps to exactly one term."""
    text = str(value or "").strip()
    if not text:
        return None
    number = re.fullmatch(r"(.*?)\|\s*(?:MS:)?(\d{7})", text)
    if number:
        return match_mass_spectrometer(None, allowed_ms_terms, accession="MS:" + number.group(2))
    return match_mass_spectrometer(text, allowed_ms_terms)


def dataset_instrument(instruments, allowed_ms_terms):
    """The allowed MassSpectrometer term for a dataset's instrument list, or None."""
    named = [i for i in instruments if not _is_placeholder(i)]
    if len(named) == 1:
        return match_mass_spectrometer(named[0].get("name"), allowed_ms_terms,
                                       accession=named[0].get("accession"))
    if named:
        return None
    # only generic 'instrument model' entries: use their free text when all of them name
    # the same instrument
    terms = {free_text_instrument(i.get("value"), allowed_ms_terms) for i in instruments if i.get("value")}
    if len(terms) != 1 or None in terms:
        return None
    return next(iter(terms))


def main():
    parser = argparse.ArgumentParser(description="Fill missing MassSpectrometer from MassIVE dataset info")
    parser.add_argument("input_metadata")
    parser.add_argument("output_metadata")
    parser.add_argument("--path_to_allowed_term_json", required=True)
    parser.add_argument("--request_pause", type=float, default=0.2, help="pause between PROXI requests (s)")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_metadata, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        shutil.copyfile(args.input_metadata, args.output_metadata)
        return

    accession_col = "ATTRIBUTE_DatasetAccession" if "ATTRIBUTE_DatasetAccession" in df.columns else "MassiveID"
    if "MassSpectrometer" not in df.columns or accession_col not in df.columns:
        print("No MassSpectrometer/dataset column; nothing to fill")
        df.to_csv(args.output_metadata, sep="\t", index=False)
        return

    with open(args.path_to_allowed_term_json) as f:
        allowed_ms_terms = json.load(f)["MassSpectrometer"]["allowed_values"]

    missing = df["MassSpectrometer"].isin([MISSING, ""])
    datasets = sorted(set(df.loc[missing, accession_col]) - {""})
    datasets = [d for d in datasets if d.startswith("MSV")]
    print(f"{missing.sum()} rows in {len(datasets)} datasets have no MassSpectrometer")

    session = requests.Session()
    filled_by_dataset = {}
    not_fetched = 0
    for accession in datasets:
        instruments = fetch_instruments(session, accession)
        if instruments is None:
            not_fetched += 1
        else:
            term = dataset_instrument(instruments, allowed_ms_terms)
            if term:
                filled_by_dataset[accession] = term
        time.sleep(args.request_pause)

    fill = missing & df[accession_col].isin(filled_by_dataset)
    df.loc[fill, "MassSpectrometer"] = df.loc[fill, accession_col].map(filled_by_dataset)
    print(f"Filled MassSpectrometer for {fill.sum()} rows in {len(filled_by_dataset)} datasets "
          f"({not_fetched} datasets could not be fetched)")
    df.to_csv(args.output_metadata, sep="\t", index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # Optional enrichment: never block the update
        print(f"WARNING: MassSpectrometer enrichment failed ({type(e).__name__}: {e}); passing rows through")
        shutil.copyfile(sys.argv[1], sys.argv[2])
