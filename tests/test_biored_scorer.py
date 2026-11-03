import csv
import importlib.util
import itertools
import tempfile
import unittest
from pathlib import Path

_spec = importlib.util.spec_from_file_location("biored_evaluate", Path(__file__).parents[1] / "experiments" / "biored" / "evaluate.py")
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)

TEXT = "PAH and GCH1 in PKU and classical PKU and BH4 deficiency"  # one block; spans below index into it


def span(word, nth=0):
    start = -1
    for _ in range(nth + 1):
        start = TEXT.index(word, start + 1)
    return start, start + len(word)


def write(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def outputs(folder, entities, pairs):
    """entities: {id: (type, [(start, end)])}; pairs: [(gene_id, disease_id, status)]."""
    write(folder / "entities.csv", [{"entity_id": k, "entity_type": t} for k, (t, _) in entities.items()], ["entity_id", "entity_type"])
    write(folder / "auxiliary_entities.csv", [], ["entity_id", "entity_type"])
    write(folder / "mentions.csv", [{"entity_id": k, "block_id": "b0001", "start": s, "end": e}
                                    for k, (_, spans) in entities.items() for s, e in spans], ["entity_id", "block_id", "start", "end"])
    fields = ["pair_id", "gene_id", "gene_name", "disease_id", "disease_name", "status", "relation_type", "source_scope", "reason"]
    write(folder / "pairs.csv", [{"pair_id": f"{g}-{d}", "gene_id": g, "gene_name": g, "disease_id": d, "disease_name": d,
                                  "status": st, "relation_type": "association", "source_scope": "本研究结果", "reason": ""}
                                 for g, d, st in pairs], fields)


def gold(concepts, relations):
    """concepts: {id: (type, [(start, end)])}; relations: [(gene concept, disease concept)]."""
    kind = {"gene": "GeneOrGeneProduct", "disease": "DiseaseOrPhenotypicFeature"}
    return {"mentions": [{"identifier": c, "type": kind[t], "text": TEXT[s:e], "block_id": "b0001", "start": s, "end": e}
                         for c, (t, spans) in concepts.items() for s, e in spans],
            "relations": [{"gene": g, "disease": d, "type": "Association", "novel": "Novel"} for g, d in relations]}


class OneToOneScoringTests(unittest.TestCase):
    def score(self, entities, pairs, concepts, relations):
        with tempfile.TemporaryDirectory() as tmp:
            outputs(Path(tmp), entities, pairs)
            counts, rows = E.score_paper("P1", gold(concepts, relations), Path(tmp))
        # Counting identities hold whatever the input.
        self.assertEqual(counts["tp_ours"], counts["tp_gold"])
        self.assertEqual(counts["tp_ours"] + sum(r["outcome"] == "FP" for r in rows), counts["our_supported"])
        self.assertEqual(counts["tp_gold"] + counts["fn"], len(relations))
        return counts, rows

    def test_two_pairs_for_one_gold_relation_count_once(self):
        entities = {"G1": ("gene", [span("PAH")]), "D1": ("disease", [span("PKU")]), "D2": ("disease", [span("classical PKU")])}
        concepts = {"g": ("gene", [span("PAH")]), "d": ("disease", [span("PKU"), span("PKU", 1)])}
        counts, rows = self.score(entities, [("G1", "D1", "supported"), ("G1", "D2", "supported")], concepts, [("g", "d")])
        self.assertEqual((counts["tp_gold"], counts["our_supported"]), (1, 2))
        self.assertEqual([r["detail"] for r in rows if r["outcome"] == "FP"], ["duplicate: its gold relation is credited to another of our pairs"])
        self.assertEqual(counts["tp_ours_many_to_one"], 2)  # the earlier protocol credited both

    def test_one_pair_matching_two_gold_relations_counts_once(self):
        both = (span("PKU")[0], span("deficiency")[1])  # one of our entities spans two gold disease concepts
        entities = {"G1": ("gene", [span("PAH")]), "D1": ("disease", [both])}
        concepts = {"g": ("gene", [span("PAH")]), "d1": ("disease", [span("PKU")]), "d2": ("disease", [span("BH4 deficiency")])}
        counts, rows = self.score(entities, [("G1", "D1", "supported")], concepts, [("g", "d1"), ("g", "d2")])
        self.assertEqual((counts["tp_gold"], counts["fn"]), (1, 1))
        self.assertIn("our pair: supported, credited to another gold relation", [r["detail"] for r in rows])

    def test_maximum_matching_not_greedy_and_order_free(self):
        # Gold r1 can take pair A or B, r2 only A: a greedy pass that gives r1 its first candidate (A) scores 1; the maximum is 2.
        entities = {"G1": ("gene", [span("PAH")]), "D1": ("disease", [span("PKU")]), "D2": ("disease", [span("classical PKU")])}
        concepts = {"g": ("gene", [span("PAH")]), "d1": ("disease", [span("PKU"), span("classical PKU")]), "d2": ("disease", [span("PKU")])}
        pairs = [("G1", "D1", "supported"), ("G1", "D2", "supported")]
        results = set()
        for rel_order in itertools.permutations([("g", "d1"), ("g", "d2")]):
            for pair_order in itertools.permutations(pairs):
                counts, _ = self.score(entities, list(pair_order), concepts, list(rel_order))
                results.add((counts["tp_gold"], counts["fn"], counts["our_supported"]))
        self.assertEqual(results, {(2, 0, 2)})

    def test_empty_gold_relations_make_every_supported_pair_a_false_positive(self):
        entities = {"G1": ("gene", [span("PAH")]), "D1": ("disease", [span("PKU")])}
        concepts = {"g": ("gene", [span("PAH")]), "d": ("disease", [span("PKU")])}
        counts, rows = self.score(entities, [("G1", "D1", "supported")], concepts, [])
        self.assertEqual((counts["tp_gold"], counts["fn"], counts["our_supported"]), (0, 0, 1))
        self.assertEqual([r["detail"] for r in rows], ["gold has both concepts but no relation"])

    def test_no_predictions_make_every_gold_relation_a_miss(self):
        entities = {"G1": ("gene", [span("PAH")]), "D1": ("disease", [span("PKU")])}
        concepts = {"g": ("gene", [span("PAH")]), "d": ("disease", [span("PKU")])}
        counts, rows = self.score(entities, [("G1", "D1", "uncertain")], concepts, [("g", "d")])
        self.assertEqual((counts["tp_gold"], counts["fn"], counts["our_supported"]), (0, 1, 0))
        self.assertEqual(rows[0]["detail"], "our pair: uncertain")
        self.assertEqual(E.prf(counts), (0.0, 0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
