"""Fetch Toxoplasma knowledge-base sources into data/raw/ and record provenance.

Sources (all public APIs, standard library only):
  1. PubMed   - reviews / meta-analyses / systematic reviews about Toxoplasma (NCBI E-utilities)
  2. openFDA  - official drug labels for pyrimethamine and sulfadiazine
  3. UniProt  - reviewed (Swiss-Prot) Toxoplasma protein entries

Outputs (JSONL, one record per line) plus data/raw/MANIFEST.json describing
exactly which URL/query each file came from, when, and how many records.

Usage:  python scripts/fetch_sources.py
"""
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "toxo-chatbot-research/0.1 (student research project)"}


def http_get(url, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read(), dict(r.headers)
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))


def write_jsonl(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


# ---------------------------------------------------------------- PubMed
PUBMED_QUERY = (
    "toxoplasm*[tiab] AND hasabstract AND "
    "(systematic[sb] OR meta-analysis[pt] OR review[pt])"
)


def fetch_pubmed():
    q = urllib.parse.urlencode(
        {"db": "pubmed", "term": PUBMED_QUERY, "retmode": "json",
         "retmax": 10000, "usehistory": "n"})
    search_url = f"{EUTILS}/esearch.fcgi?{q}"
    body, _ = http_get(search_url)
    res = json.loads(body)["esearchresult"]
    ids = res["idlist"]
    print(f"[pubmed] {res['count']} matches, fetching {len(ids)}")

    records = []
    for i in range(0, len(ids), 200):
        batch = ids[i:i + 200]
        url = (f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id="
               + ",".join(batch))
        xml_bytes, _ = http_get(url)
        root = ET.fromstring(xml_bytes)
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//MedlineCitation/PMID")
            title = "".join(art.find(".//ArticleTitle").itertext()) \
                if art.find(".//ArticleTitle") is not None else ""
            parts = []
            for ab in art.findall(".//Abstract/AbstractText"):
                label = ab.get("Label")
                text = "".join(ab.itertext()).strip()
                parts.append(f"{label}: {text}" if label else text)
            year = (art.findtext(".//JournalIssue/PubDate/Year")
                    or (art.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4])
            doi = None
            for aid in art.findall(".//ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = aid.text
            records.append({
                "source": "PubMed",
                "pmid": pmid,
                "title": title,
                "abstract": "\n".join(parts),
                "journal": art.findtext(".//Journal/Title"),
                "year": year,
                "doi": doi,
                "pub_types": [p.text for p in art.findall(".//PublicationType")],
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })
        time.sleep(0.4)  # stay under 3 requests/second without an API key
        print(f"[pubmed] {min(i + 200, len(ids))}/{len(ids)}")

    out = RAW / "pubmed" / "toxoplasma_reviews.jsonl"
    write_jsonl(out, records)
    return {
        "file": str(out.relative_to(ROOT)),
        "records": len(records),
        "query": PUBMED_QUERY,
        "api": "NCBI E-utilities esearch + efetch",
        "search_url": search_url,
        "total_matches_reported": int(res["count"]),
        "licence_note": "Abstracts are publisher-copyrighted; use locally for "
                        "research and cite PMIDs; do not redistribute.",
    }


# --------------------------------------------------------------- openFDA
FDA_DRUGS = ["pyrimethamine", "sulfadiazine"]
FDA_FIELDS = [
    "indications_and_usage", "dosage_and_administration", "contraindications",
    "warnings", "warnings_and_cautions", "precautions", "pregnancy",
    "nursing_mothers", "pediatric_use", "geriatric_use", "adverse_reactions",
    "drug_interactions", "clinical_pharmacology", "description", "overdosage",
]


def fetch_openfda():
    records, urls, dropped = [], [], []
    for drug in FDA_DRUGS:
        url = ("https://api.fda.gov/drug/label.json?search="
               + urllib.parse.quote(f"openfda.generic_name:{drug}") + "&limit=100")
        urls.append(url)
        body, _ = http_get(url)
        data = json.loads(body)
        for r in data.get("results", []):
            of = r.get("openfda", {})
            # The generic_name search for "sulfadiazine" also returns SILVER
            # sulfadiazine (topical burn cream). That is NOT the oral drug used
            # for toxoplasmosis, so exclude it to avoid wrong retrieval.
            generic = " ".join(of.get("generic_name", [])).lower()
            routes = [x.lower() for x in of.get("route", [])]
            if "silver" in generic or (drug == "sulfadiazine" and "oral" not in routes):
                dropped.append({"set_id": r.get("set_id"),
                                "generic_name": of.get("generic_name"),
                                "route": of.get("route")})
                continue
            rec = {
                "source": "openFDA drug label",
                "query_drug": drug,
                "route": of.get("route"),
                "set_id": r.get("set_id"),
                "effective_time": r.get("effective_time"),
                "generic_name": of.get("generic_name"),
                "brand_name": of.get("brand_name"),
                "manufacturer": of.get("manufacturer_name"),
                "url": f"https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid={r.get('set_id')}",
            }
            for field in FDA_FIELDS:
                if field in r:
                    rec[field] = "\n".join(r[field])
            records.append(rec)
        print(f"[openfda] {drug}: {len(data.get('results', []))} labels returned")
        time.sleep(0.5)
    print(f"[openfda] kept {len(records)}, dropped {len(dropped)} (silver/non-oral)")

    out = RAW / "openfda" / "toxoplasma_drug_labels.jsonl"
    write_jsonl(out, records)
    return {
        "file": str(out.relative_to(ROOT)),
        "records": len(records),
        "drugs": FDA_DRUGS,
        "api": "openFDA drug/label endpoint (labels originate from DailyMed)",
        "urls": urls,
        "dropped_labels": dropped,
        "filter_note": "Dropped silver sulfadiazine (topical) and any sulfadiazine "
                       "label without an oral route.",
        "licence_note": "openFDA data is public domain (CC0). Labels are US "
                        "prescribing information only; not every Toxoplasma "
                        "regimen is on-label.",
    }


# --------------------------------------------------------------- UniProt
UNIPROT_QUERY = "taxonomy_id:5811 AND reviewed:true"
UNIPROT_FIELDS = ("accession,id,gene_names,protein_name,organism_name,length,"
                  "cc_function,cc_subcellular_location,go_f,go_p,go_c,xref_pdb")


def fetch_uniprot():
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(UNIPROT_QUERY)
           + "&format=tsv&fields=" + UNIPROT_FIELDS + "&size=500")
    body, _ = http_get(url)
    lines = body.decode("utf-8").splitlines()
    header = lines[0].split("\t")
    records = []
    for line in lines[1:]:
        vals = line.split("\t")
        rec = dict(zip(header, vals))
        rec["source"] = "UniProtKB/Swiss-Prot"
        rec["url"] = f"https://www.uniprot.org/uniprotkb/{rec.get('Entry')}/entry"
        records.append(rec)
    print(f"[uniprot] {len(records)} reviewed entries")

    out = RAW / "uniprot" / "toxoplasma_reviewed_proteins.jsonl"
    write_jsonl(out, records)
    return {
        "file": str(out.relative_to(ROOT)),
        "records": len(records),
        "query": UNIPROT_QUERY,
        "api": "UniProt REST API (uniprotkb/search, TSV)",
        "url": url,
        "licence_note": "UniProt is CC BY 4.0; cite UniProt.",
    }


def main():
    import sys
    available = {"pubmed": fetch_pubmed, "openfda": fetch_openfda,
                 "uniprot": fetch_uniprot}
    wanted = sys.argv[1:] or list(available)  # e.g. `python fetch_sources.py openfda`
    manifest_path = RAW / "MANIFEST.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        manifest = {"sources": {}}
    manifest["generated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for name in wanted:
        try:
            result = available[name]()
            result["fetched_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
            manifest["sources"][name] = result
        except Exception as e:  # noqa: BLE001
            manifest["sources"][name] = {"error": repr(e)}
            print(f"[{name}] FAILED: {e!r}")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print("Wrote data/raw/MANIFEST.json")


if __name__ == "__main__":
    main()

