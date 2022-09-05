# Grab example data and create a skeleton solution

from bs4 import BeautifulSoup
from string import Template
import requests
import sys

code = sys.argv[1].rstrip()

link = "https://rosalind.info/problems/" + code
page = requests.get(link)
soup = BeautifulSoup(page.content, "html.parser")
problem, solution = [
    txt.get_text() for txt in soup.find_all("div", class_="codehilite")
]
title = soup.find("h1").get_text().splitlines()[0]

# Create test data
print(f"Creating sample dataset file: rosalind_{code}.txt ...")
soup.find(lambda tag: tag.name == "a" and "List View" in tag.text)
with open(f"rosalind_{code}.txt", "w") as f:
    f.write(problem)

# Create expected solution
print(f"Creating sample output file: {code}.txt ...")
with open(f"{code}.txt", "w") as f:
    f.write(solution)

# Create template code
print(f"Creating skeleton solution: {code}.py ...")
s = Template(open("template.txt").read())
o = s.substitute(code=code, title=title)
with open(f"{code}.py", "w") as f:
    f.write(o)
