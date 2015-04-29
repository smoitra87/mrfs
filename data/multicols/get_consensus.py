from lxml import html
import requests
import os, sys

def fetch_consensus(pfamid):
    page = requests.get('http://gremlin.bakerlab.org/pfam.php?id={}&year=2013'.format(pfamid))
    tree = html.fromstring(page.text)
    consensus = tree.xpath('//*[@id="textarea"]/text()')
    print "{} found sequence {}".format(pfamid, consensus)

    baker_file  = os.path.join(pfamid, "baker.fasta")
    with open(baker_file,'w') as fout:
        print >>fout, ">baker"
        print >>fout, consensus[0]

    return consensus
