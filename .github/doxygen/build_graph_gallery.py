"""Embed Doxygen's generated SVG diagrams in the documentation's graph gallery."""

import html
import re
import sys
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.links = {}
        self.graphs = []
        self.link = None
        self.text = []
        self.feed(path.read_text(encoding="utf-8"))

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and attrs.get("href"):
            self.link = attrs["href"]
            self.text = []
        if tag in ("iframe", "object"):
            source = attrs.get("src", attrs.get("data", ""))
            if source.endswith(".svg"):
                self.graphs.append(attrs | {"source": source})

    def handle_data(self, text):
        if self.link is not None:
            self.text.append(text)

    def handle_endtag(self, tag):
        if tag == "a" and self.link is not None:
            self.links["".join(self.text).strip()] = self.link
            self.link = None


def diagram(root, graph, title, reference=None):
    source = graph["source"]
    if not (root / source).is_file():
        raise RuntimeError(f"Missing graph: {source}")
    # Small Doxygen graphs have fixed SVG dimensions; fit them to mobile frames.
    svg = ET.parse(root / source)
    element = svg.getroot()
    if element.get("viewBox"):
        element.set("width", "100%")
        element.set("height", "100%")
        element.set("preserveAspectRatio", "xMidYMid meet")
        source = "gallery_" + Path(source).name
        svg.write(root / source, encoding="utf-8", xml_declaration=True)
    width_match = re.fullmatch(r"([0-9.]+)(?:px)?", graph.get("width", ""))
    height_match = re.fullmatch(r"([0-9.]+)(?:px)?", graph.get("height", ""))
    width = float(width_match[1]) if width_match else 850
    height = float(height_match[1]) if height_match else 400
    height = round(max(150, min(600, height * min(1, 850 / width) + 24)))
    links = f'<a href="{html.escape(source)}" target="_blank" rel="noopener">Open full-size diagram</a>'
    if reference:
        links += f' &middot; <a href="{html.escape(reference)}">Open reference</a>'
    return (
        f"<h3>{html.escape(title)}</h3>\n"
        '<figure class="as3-diagram">\n'
        f'<iframe src="{html.escape(source)}" title="{html.escape(title)}" '
        f'loading="lazy" height="{height}"></iframe>\n'
        f"<figcaption>{links}</figcaption>\n</figure>\n"
    )


def main(root):
    page = root / "diagrams.html"
    source = page.read_text(encoding="utf-8")
    class_links = Page(root / "annotated.html").links
    file_links = Page(root / "files.html").links
    sections = {}
    inheritance = []
    families = {
        "ISolver": "Zone solvers", "ITemporal": "Time integration",
        "IBoundary": "Boundary conditions", "IInterface": "Zone interfaces",
        "IRiemannSolver": "Riemann solvers", "ITensorProduct": "Tensor products",
        "IInitialCondition": "Initial conditions", "IFileVTK": "Visualization output",
        "IBoundaryParamMarker": "Boundary marker parameters",
    }
    ns = {"svg": "http://www.w3.org/2000/svg"}
    for graph in Page(root / "inherits.html").graphs:
        svg = ET.parse(root / graph["source"])
        nodes = svg.findall('.//svg:g[@class="node"]', ns)
        if len(nodes) < 2:
            continue
        labels = [" ".join(t.itertext()) for node in nodes for t in node.findall(".//svg:text", ns)]
        title = next((title for name, title in families.items() if name in labels), " / ".join(labels))
        inheritance.append((title, diagram(root, graph, title, "inherits.html")))
    if not inheritance:
        raise RuntimeError("No inheritance relationships were generated")
    sections["INHERITANCE"] = "".join(markup for _, markup in sorted(inheritance))

    collaboration = []
    for name in ("CDriver", "CEESolver", "CMultizoneGeometry"):
        reference = class_links[name]
        graphs = [g for g in Page(root / reference).graphs if g["source"].endswith("__coll__graph.svg")]
        if not graphs:
            raise RuntimeError(f"No class dependency graph was generated for {name}")
        collaboration.append(diagram(root, graphs[0], name, reference))
    sections["COLLABORATION"] = "".join(collaboration)

    includes = []
    for name, suffix, title in (
        ("driver_structure.hpp", "__incl.svg", "Driver: included headers"),
        ("geometry_structure.hpp", "__dep__incl.svg", "Geometry: headers and sources that include it"),
        ("solver_structure.hpp", "__incl.svg", "Solver: included headers"),
    ):
        reference = file_links[name]
        graphs = [g for g in Page(root / reference).graphs if g["source"].endswith(suffix) and (suffix != "__incl.svg" or "__dep__incl.svg" not in g["source"])]
        if not graphs:
            raise RuntimeError(f"No header dependency graph was generated for {name}")
        includes.append(diagram(root, graphs[0], title, reference))
    sections["INCLUDE"] = "".join(includes)

    for name, markup in sections.items():
        pattern = rf"(<!-- AS3_{name}_START -->).*?(<!-- AS3_{name}_END -->)"
        source, count = re.subn(pattern, lambda match: match[1] + "\n" + markup + match[2], source, flags=re.S)
        if count != 1:
            raise RuntimeError(f"Expected one {name.lower()} gallery placeholder in {page}")
    page.write_text(source, encoding="utf-8")
    print(f"Graph gallery: {len(inheritance)} inheritance families, {len(collaboration)} class dependency graphs, {len(includes)} header dependency graphs")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "build/docs/html"))
