import time
import tqdm
import pandas
import requests

if __name__ == "__main__":
    website_url = "http://rebase.neb.com"
    methylations = ['Dam', 'Dcm', 'CpG', 'EcoBI', 'EcoKI']
    columns = ['enzyme', 'req_seq', 'suppliers'] + methylations

    print("Getting enzyme list and suppliers...")

    html = requests.get(website_url + '/cgi-bin/azlist?re2+cy').text
    tables = pandas.read_html(html, header=0)
    table = [t for t in tables if (t.values.shape[1] == 4)][0]
    suppliers_dict = {
        enzyme: len(suppliers)
        for enzyme, suppliers in zip(table.Enzymes, table.Suppliers)
    }

    print("Getting enzymes methylation data...")

    def get_enzyme_data(enzyme_name):
        url = website_url + "/cgi-bin/damlist?e" + enzyme_name
        html = requests.get(url).text
        table = [
            t for t in pandas.read_html(html)
            if t.values.shape == (3, 6)
        ][0]
        rec_seq = table.values[0, 0].replace(enzyme_name, '')
        suppliers = suppliers_dict[enzyme_name]
        return dict([('enzyme', enzyme_name),
                    ('req_seq', rec_seq),
                    ('suppliers', suppliers)] +
                    list(zip(methylations, table.values[2, 1:6])))
    data = []
    errored = []
    for enzyme_name in tqdm.tqdm(list(suppliers_dict)):
        time.sleep(0.5)
        try:
            data.append(get_enzyme_data(enzyme_name))
        except Exception as e:
            errored.append(enzyme_name)
    if len(errored):
        print ("%d errors:" % len(errored), " ".join(errored))
    
    print("Generating the spreadsheet...")

    dataframe = pandas.DataFrame(data, columns=columns)
    dataframe.to_csv("./enzymes_infos.csv", index=False)

    print ('Done !')