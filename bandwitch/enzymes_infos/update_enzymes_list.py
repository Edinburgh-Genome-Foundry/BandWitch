from lxml import html
import requests
from collections import OrderedDict, Counter

if __name__ == "__main__":
    website_url = "http://rebase.neb.com"

    # GET THE CONTENT FROM ALL PROVIDER PAGES
    all_providers_page = requests.get(website_url + "/cgi-bin/lapsupplist")
    parsed_providers_page = html.fromstring(all_providers_page.text)
    providers_pages = OrderedDict(sorted([
        (e.text, requests.get((website_url + e.values()[0])))
        for e in parsed_providers_page.xpath('//a')
        if e.values()[0].startswith("/cgi-bin/laplist")
    ]))

    # CREATE A TABLE WITH ALL ENZYMES FROM ALL PAGES
    full_table = []
    for provider, page in providers_pages.items():
        parsed_page = html.fromstring(page.text)
        tables = parsed_page.xpath('//table[@border=2]')
        for table in tables:
            parsed_table = [
                [v.text_content() for v in row]
                for row in table[2:]
            ]
            if (len(parsed_table) == 0) or (parsed_table[0][0] != 'overlaps:'):
                continue
            header, parsed_table = parsed_table[0], parsed_table[1:]
            full_table += [tuple(row) for row in parsed_table]

    # REMOVE (AND COUNT) DUPLICATES BETWEEN PROVIDERS. WRITE TO A CSV
    counter = Counter(full_table)
    enzymes = [k[0] for k in counter.keys()]
    assert len(enzymes) == len(set(enzymes))
    with open("enzymes_data.csv", "w+") as f:
        f.write("\n".join(
            [";".join(header + ['suppliers'])] + sorted([
                ";".join(line + (str(n),))
                for line, n in counter.items()
            ])))
