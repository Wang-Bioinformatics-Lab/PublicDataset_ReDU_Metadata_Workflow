"""Loads the pinned OWL files with owlready2, without network access.

owlready2 follows the owl:imports of an ontology and downloads every import it does not
find in owlready2.onto_path, on every load: their latest version, and a failed download
(e.g. a purl.obolibrary.org 502) stops the step. PO and DOID have imports (DOID's ext.owl
imports further files); they are pinned in data/reference_data.json at the release of the
importing ontology and installed in data/<ontology>_imports/, each file under the last
part of its import IRI, which is where owlready2 looks for it.

    load_ontology(owl_path, imports_dir=None) -> owlready2 Ontology

Imports are read from imports_dir only; an import that is not there is an error, never a
download.
"""
import owlready2
import owlready2.namespace as _owl_ns

_get_onto_file = _owl_ns._get_onto_file


# owlready2 0.46 resolves every file it loads here: a local file, or else the IRI itself,
# which it then downloads
def _local_onto_file(base_iri, name, mode='r', only_local=False):
    path = _get_onto_file(base_iri, name, mode, only_local)
    if mode.startswith('r') and path.startswith(('http://', 'https://')):
        raise FileNotFoundError(f"{path} is imported but not in {owlready2.onto_path or 'any imports folder'}; "
                                f"pin it in data/reference_data.json")
    return path
_owl_ns._get_onto_file = _local_onto_file


def load_ontology(owl_path, imports_dir=None):
    # Each ontology gets its own owlready2 World: in the shared default world, classes an
    # earlier file already declared (e.g. UBERON declares ENVO's 'marine biome') belong to
    # that earlier ontology and silently disappear from this one's classes()
    owlready2.onto_path[:] = [imports_dir] if imports_dir and imports_dir != 'none' else []
    return owlready2.World().get_ontology(owl_path).load()
