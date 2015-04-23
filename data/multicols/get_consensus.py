from lxml import html
import requests

page = requests.get('http://gremlin.bakerlab.org/pfam.php?id=PF12569&year=2013')
tree = html.fromstring(page.text)

consensus = tree.xpath('//*[@id="textarea"]/text()')

print "found sequence", consensus

with open("baker.fasta",'w') as fout:
    print >>fout, ">baker"
    print >>fout, consensus[0]
