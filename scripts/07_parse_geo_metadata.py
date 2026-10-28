from pathlib import Path
import csv

INPUT = Path("data/metadata/GSE72056_family.soft")
OUTPUT = Path("data/metadata/GSE72056_sample_metadata.csv")

records = []
current = None

fields_to_keep = {
    "!Sample_title": "title",
    "!Sample_source_name_ch1": "source_name",
    "!Sample_organism_ch1": "organism",
    "!Sample_platform_id": "platform_id",
}

with INPUT.open("r", encoding="utf-8", errors="replace") as handle:
    for raw_line in handle:
        line = raw_line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            if current is not None:
                records.append(current)

            current = {
                "geo_accession": line.split("=", 1)[1].strip(),
                "title": "",
                "source_name": "",
                "organism": "",
                "platform_id": "",
                "characteristics": [],
            }
            continue

        if current is None:
            continue

        for prefix, field in fields_to_keep.items():
            if line.startswith(prefix + " = "):
                value = line.split(" = ", 1)[1].strip()
                if current[field]:
                    current[field] += " | " + value
                else:
                    current[field] = value

        if line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split(" = ", 1)[1].strip()
            current["characteristics"].append(value)

if current is not None:
    records.append(current)

for record in records:
    record["characteristics"] = " | ".join(record["characteristics"])

columns = [
    "geo_accession",
    "title",
    "source_name",
    "organism",
    "platform_id",
    "characteristics",
]

with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(records)

print("GEO sample records parsed:", len(records))
print("Saved:", OUTPUT)
print("\nFirst 20 records:")

for record in records[:20]:
    print(
        record["geo_accession"],
        "|", record["title"],
        "|", record["source_name"],
        "|", record["characteristics"][:180],
    )
