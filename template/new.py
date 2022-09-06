# Grab example data and create a skeleton solution
# Input should be the path to the python script we want.

from bs4 import BeautifulSoup
from string import Template
import requests
import sys
from pathlib import Path

code = Path(sys.argv[1]).stem
location = Path(sys.argv[1]).parent

link = "https://rosalind.info/problems/" + code
page = requests.get(link)
soup = BeautifulSoup(page.content, "html.parser")
problem, solution = [
    txt.get_text() for txt in soup.find_all("div", class_="codehilite")
]
title = soup.find("h1").get_text().splitlines()[0]

# Create test data
print(f"Creating sample dataset file: {location}/test/rosalind_{code}.txt ...")
soup.find(lambda tag: tag.name == "a" and "List View" in tag.text)
with open(f"{location}/test/rosalind_{code}.txt", "w") as f:
    f.write(problem)

# Create expected solution
print(f"Creating sample output file: tests/snapshots/{location}/{code}.txt ...")
with open(f"tests/snapshots/{location}/{code}.txt", "w") as f:
    f.write(solution)

# Create template code
print(f"Creating skeleton solution: {location}/{code}.py ...")
s = Template(open(str(Path(__file__).parent) + "/template.txt").read())
o = s.substitute(code=code, title=title)
with open(f"{location}/{code}.py", "w") as f:
    f.write(o)
