import csv

ob = set()
discontinued = set()
with open("products.txt", encoding="latin-1") as f:
    for row in csv.DictReader(f, delimiter="~"):
        for ing in row["Ingredient"].split(";"):
            ing = ing.strip().upper()
            if row["Type"] == "DISCN":
                discontinued.add(ing)
            else:
                ob.add(ing)

for r in csv.DictReader(open("hits.csv")):
    n = (r["name"] or "").upper()
    if any(n in i or i in n for i in ob):
        status = "FDA"
    elif any(n in i or i in n for i in discontinued):
        status = "FDA-discontinued"
    else:
        status = "NOT-FDA"
    print(r["chembl_id"], r["name"], r["score"], status)